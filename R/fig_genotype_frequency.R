# ###############################################################
# COVARIATES AND MAPS
# ###############################################################
# PART A   covariates, run FIRST, writes data_raw/kigali_covariates.csv
# PART B   maps, run AFTER the simulation
#
# The order is always: Part A, simulation, Part B. The simulation reads the CSV
# Part A writes, and Part B reads what the simulation writes.
#
# Which package does what:
#   terra      reads the .tif files, and disagg / crop / mask for the grids.
#              raster cannot read GDAL band descriptions, so "nets_2010" would
#              never match, and it has no touches argument on mask().
#   sf         the boundary, the district lookup, the outlines on the figures
#   dplyr,tidyr  tables, joins, reshaping
#   plain R    the Gaussian smoother: outer() and one matrix multiply.
#              sf.kde gave a density, which is not bounded by the values it is
#              built from, so it pulled the frequencies towards zero. Dropped,
#              and spatialEco with it.
#
# Run the whole file with source(), never block by block.

library(terra)
library(sf)
library(geodata)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)


# ###############################################################
# PART A  -  COVARIATES
# ###############################################################
covariates_insec<- 
#############################################################
### A0- Configuration
#############################################################
YEAR       <- 2010
STUDY_AREA <- "kigali"

# The simulation reads this exact path. Change it here and you change it there.
OUT_CSV <- "data_raw/kigali_covariates.csv"

# annual cubes give `prefix`, the band is prefix_YEAR
# static layers give `bands`, and every band is summed
CUBES <- list(
  itn  = list(path = "data_raw/raster_clean/net_use_scaled_cube.tif",      prefix = "nets"),
  irs  = list(path = "data_raw/raster_clean/irs_coverage_scaled_cube.tif", prefix = "irs"),
  crop = list(path = "data_raw/raster_clean/crop_scaled.tif",              bands  = "ALL"),
  pop  = list(path = "data_raw/raster_clean/pop_scaled_cube.tif",          prefix = "pop")
)

LOG_FIRST <- c("pop")     # right-skewed: without the log one urban cell takes
# the whole upper half of the scale

# Real ITN over Kigali spans 0.60 to 0.73 and real IRS is flat, so neither can
# identify betamat. Both are replaced. They must not look alike or their
# coefficients trade off: the districts run northeast to southwest, so a
# north-south IRS gave cor = 0.79. West-east runs across that instead.
#   "real" | "district" | "gradient_we" | "gradient_ns"
OVERRIDE <- list(itn = "district", irs = "gradient_we", crop = "real", pop = "real")

SMOOTH_BW       <- 0.055    # degrees, MUST exceed the cell width of 0.0417
MODEL_FACTOR    <- 1        # 1 keeps the native grid. 2 quadruples the pixels.
TEMPLATE_FACTOR <- 10       # display grid, figures only
SCALE_RANGE     <- c(0.02, 0.98)
KEEP_EDGE       <- TRUE     # keep cells the border only clips, so the thin
# northern tip of Kigali gets a pixel

dir.create("data_raw/gadm", recursive = TRUE, showWarnings = FALSE)
dir.create("dataoutput",    recursive = TRUE, showWarnings = FALSE)


#############################################################
### A1- Boundary
#############################################################
# GADM writes the province in Kinyarwanda, "Umujyi wa Kigali", so an exact match
# on "Kigali City" returns zero rows and fails much later as an NA bbox.
ADMIN2_SHP <- "data_raw/rwanda_admin2.shp"

if (!file.exists(ADMIN2_SHP)) {
  gadm("RWA", level = 2, path = "data_raw/gadm") %>%
    st_as_sf() %>% st_transform(4326) %>%
    st_write(ADMIN2_SHP, delete_dsn = TRUE, quiet = TRUE)
}

admin_sel <- st_read(ADMIN2_SHP, quiet = TRUE) %>%
  st_transform(4326) %>%
  filter(grepl("Kigali", NAME_1, ignore.case = TRUE))

stopifnot(nrow(admin_sel) > 0)

kigali_shp  <- st_union(admin_sel) %>% st_sf()
kigali_vect <- vect(kigali_shp)

bb <- st_bbox(kigali_shp)

outside <- st_difference(
  st_as_sfc(st_bbox(c(xmin = bb[["xmin"]] - 0.2, ymin = bb[["ymin"]] - 0.2,
                      xmax = bb[["xmax"]] + 0.2, ymax = bb[["ymax"]] + 0.2),
                    crs = 4326)),
  st_union(st_geometry(kigali_shp))
)

# the white polygon is bigger than Kigali, so hold the view on Kigali itself
map_coord <- function() {
  coord_sf(xlim = c(bb[["xmin"]], bb[["xmax"]]),
           ylim = c(bb[["ymin"]], bb[["ymax"]]), expand = FALSE)
}
#############################################################
### A2- Load the layers
#############################################################
load_layer <- function(spec, year, label) {
  r <- rast(spec$path) %>% crop(kigali_vect) %>% mask(kigali_vect, touches = KEEP_EDGE)
  out <- if (is.null(spec$bands)) r[[paste0(spec$prefix, "_", year)]]
  else mask(sum(r, na.rm = TRUE), kigali_vect, touches = KEEP_EDGE)
  names(out) <- label
  out
}

layers_raw <- lapply(names(CUBES), function(k) load_layer(CUBES[[k]], YEAR, k))
names(layers_raw) <- names(CUBES)

# every layer must sit on the same grid, or the smoothing mixes geometries
ref <- layers_raw[[1]]
for (nm in names(layers_raw)) {
  if (!compareGeom(layers_raw[[nm]], ref, stopOnError = FALSE)) {
    layers_raw[[nm]] <- resample(layers_raw[[nm]], ref, method = "bilinear")
  }
}


#############################################################
### A3- Impose the designed surfaces
#############################################################
# a point inside a polygon is at distance 0, so nearest_feature gives the
# district it sits in, and the closest one for pixels on the border
surface_district <- function(xy) {
  pts <- st_as_sf(as.data.frame(xy), coords = c("x", "y"), crs = 4326)
  d   <- as.character(admin_sel$NAME_2[st_nearest_feature(pts, admin_sel)])
  setNames(seq(0.05, 0.95, length.out = length(unique(d))), sort(unique(d)))[d]
}
surface_we <- function(xy) 0.02 + 0.96 * (xy[, 1] - min(xy[, 1])) / diff(range(xy[, 1]))
surface_ns <- function(xy) 0.02 + 0.96 * (max(xy[, 2]) - xy[, 2]) / diff(range(xy[, 2]))

impose <- function(r, mode, label) {
  cells <- which(!is.na(values(r)))
  xy    <- xyFromCell(r, cells)
  new   <- switch(mode,
                  district    = surface_district(xy),
                  gradient_we = surface_we(xy),
                  gradient_ns = surface_ns(xy),
                  stop("unknown mode: ", mode))
  out <- r
  values(out)[cells] <- unname(new)
  names(out) <- label
  out
}

for (nm in names(OVERRIDE)) {
  if (OVERRIDE[[nm]] != "real") {
    layers_raw[[paste0(nm, "_real")]] <- layers_raw[[nm]]   # kept for the map
    layers_raw[[nm]] <- impose(layers_raw[[nm]], OVERRIDE[[nm]], nm)
  }
}


#############################################################
### A4- Rescale
#############################################################
rescale01 <- function(r, log_first = FALSE) {
  v   <- values(r)
  if (log_first) v <- log1p(v * 1000)     # values are tiny, lift before the log
  fin <- is.finite(v)
  rng <- range(v[fin])
  if (diff(rng) < 1e-12) return(r)
  v[fin] <- SCALE_RANGE[1] + diff(SCALE_RANGE) * (v[fin] - rng[1]) / diff(rng)
  values(r) <- v
  r
}

log_this <- function(n) sub("_real$", "", n) %in% LOG_FIRST

layers_scaled <- lapply(names(layers_raw), function(n)
  rescale01(layers_raw[[n]], log_first = log_this(n)))
names(layers_scaled) <- names(layers_raw)
cov_cols <- names(layers_scaled)


#############################################################
### A5- Smooth and write the CSV
#############################################################
# A weighted mean of the source cells, weights falling off with distance. A mean
# cannot leave the range of its inputs, so a frequency stays a frequency. Below
# the cell width the kernel never reaches a neighbour and nothing is smoothed.
stopifnot(SMOOTH_BW > res(ref)[1])

gauss_smooth <- function(target, src_lon, src_lat, src_val, bw) {
  d2 <- outer(target$longitude, src_lon, "-")^2 +
    outer(target$latitude,  src_lat,  "-")^2
  w  <- exp(-d2 / (2 * bw^2))
  as.vector(w %*% src_val) / rowSums(w)
}

smooth_to <- function(r, target, bw = SMOOTH_BW, rescale = FALSE) {
  src <- as.data.frame(r, xy = TRUE)
  names(src) <- c("longitude", "latitude", "v")
  src <- src[is.finite(src$v), ]
  out <- gauss_smooth(target, src$longitude, src$latitude, src$v, bw)
  # smoothing pulled the extremes in, so the imposed surfaces are stretched
  # back. The real ones are not: that would invent contrast the data lack.
  if (rescale) {
    rng <- range(out, na.rm = TRUE)
    if (diff(rng) > 1e-12)
      out <- SCALE_RANGE[1] + diff(SCALE_RANGE) * (out - rng[1]) / diff(rng)
  }
  out
}

needs_rescale <- function(n) !grepl("_real$", n) && OVERRIDE[[sub("_real$", "", n)]] != "real"

# disagg cuts the existing cells into smaller ones, so the fine cells sit
# exactly inside the coarse ones. rast(vect, ...) would build a fresh grid from
# the polygon bounding box, and its cells land off the source grid: that is what
# makes the colour sit half a cell away from the outline.
# build_grid <- function(fact) {
#   g <- disagg(ref, fact) %>% crop(kigali_vect)
#   values(g) <- 1
#   as.data.frame(mask(g, kigali_vect, touches = KEEP_EDGE),
#                 xy = TRUE, na.rm = TRUE)[, 1:2] %>%
#     setNames(c("longitude", "latitude"))
# }

build_grid <- function(fact, touch = KEEP_EDGE) {
  g <- disagg(ref, fact) %>% crop(kigali_vect)
  values(g) <- 1
  as.data.frame(mask(g, kigali_vect, touches = touch),
                xy = TRUE, na.rm = TRUE)[, 1:2] %>%
    setNames(c("longitude", "latitude"))
}

tgt      <- build_grid(MODEL_FACTOR)
smoothed <- lapply(cov_cols, function(n)
  smooth_to(layers_scaled[[n]], tgt, rescale = needs_rescale(n)))
names(smoothed) <- cov_cols

cov_tbl <- bind_cols(tgt, as.data.frame(smoothed))

# district name, for reporting only
cov_tbl$district <- as.character(admin_sel$NAME_2[st_nearest_feature(
  st_as_sf(cov_tbl, coords = c("longitude", "latitude"), crs = 4326), admin_sel)])
cov_tbl$village <- paste0("px_", STUDY_AREA, "_", seq_len(nrow(cov_tbl)))

write.csv(cov_tbl, "C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/data_raw/kigali_covariates.csv", row.names = FALSE)
cat("written:", OUT_CSV, nrow(cov_tbl), "rows\n")

# if two covariates look alike the model cannot tell their effects apart
model_cols <- grep("_real$", cov_cols, value = TRUE, invert = TRUE)
print(round(cor(cov_tbl[, model_cols]), 2))


#############################################################
### A6- Covariate map
#############################################################
tmpl_pts <- build_grid(TEMPLATE_FACTOR, touch = FALSE)
disp     <- lapply(cov_cols, function(n)
  smooth_to(layers_scaled[[n]], tmpl_pts, rescale = needs_rescale(n)))
names(disp) <- cov_cols
disp_tbl <- bind_cols(tmpl_pts, as.data.frame(disp))

map_theme <- function() {
  theme_minimal(base_size = 13) +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 10),
          legend.position = "bottom")
}

cov_labels <- c(itn = "ITN (imposed)", itn_real = "ITN (real)",
                irs = "IRS (imposed)", irs_real = "IRS (real)",
                crop = "Cropland",     pop = "Population")

p_cov <- disp_tbl %>%
  pivot_longer(all_of(cov_cols), names_to = "covariate", values_to = "value") %>%
  mutate(covariate = factor(covariate, levels = cov_cols,
                            labels = cov_labels[cov_cols])) %>%
  ggplot(aes(longitude, latitude, fill = value)) +
  geom_raster() +
  geom_sf(data = admin_sel, fill = NA, colour = "grey20", inherit.aes = FALSE) +
  scale_fill_gradientn(name = "Scaled value", colours = brewer.pal(9, "Purples"),
                       limits = c(0, 1), na.value = NA,
                       guide = guide_colorbar(barwidth = 16, barheight = 0.7,
                                              title.position = "top")) +
  facet_wrap(~covariate, nrow = 2) +
  coord_sf(expand = FALSE) +
  map_theme()

print(p_cov)
ggsave("C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/data_raw/kigali_covariates.png", p_cov, width = 10, height = 8, dpi = 300)


# ###############################################################
# PART B  -  MAPS
# ###############################################################
# Everything below needs ref, kigali_vect, admin_sel, gauss_smooth and cov_tbl
# from Part A, so run Part A in the same session, not just its CSV.

#############################################################
### B0- Setup
#############################################################
BW_MAP <- 0.05 # separate from SMOOTH_BW: Part A needs a bandwidth above the
# cell width, the maps do not. Smaller keeps more contrast.

base_rast  <- ref
template   <- disagg(ref, TEMPLATE_FACTOR) %>% crop(kigali_vect)
gen_labels <- c("1" = "Generation 1", "2" = "Generation 2")

template_pts <- as.data.frame(template, xy = TRUE, na.rm = FALSE)[, 1:2] %>%
  setNames(c("longitude", "latitude"))

to_long_df <- function(r) {
  as.data.frame(r, xy = TRUE, na.rm = FALSE) %>%
    pivot_longer(-c(x, y), names_to = "panel", values_to = "value")
}


#############################################################
### B1- Read the simulated data
#############################################################
mcmc_data_all_sim <- read.csv("C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/dataoutput/data_to_use.csv") %>%
  mutate(longitude = round(longitude, 4), latitude = round(latitude, 4))

allele_freq_data <- read.csv("C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/dataoutput/allele_frequency_data.csv") %>%
  left_join(distinct(mcmc_data_all_sim, village, longitude, latitude),
            by = "village") %>%
  mutate(longitude = round(longitude, 4), latitude = round(latitude, 4))

allele_loci   <- unique(allele_freq_data$locus)
marker_labels <- setNames(paste0("Marker ", seq_along(allele_loci)), allele_loci)


#############################################################
### B2- One raster per generation and per locus
#############################################################
panels <- expand_grid(timepoint = as.numeric(names(gen_labels)),
                      locus     = allele_loci) %>%
  mutate(label = paste0("Gen", timepoint, "_", marker_labels[locus]))

make_raw <- function(tp, loc) {
  d <- allele_freq_data %>% filter(locus == loc, timepoint == tp)
  r <- base_rast
  values(r) <- NA
  r[cellFromXY(r, as.matrix(d[, c("longitude", "latitude")]))] <- d$allele_frequency
  mask(r, kigali_vect, touches = FALSE)
}

make_smooth <- function(tp, loc) {
  d <- allele_freq_data %>% filter(locus == loc, timepoint == tp)
  r <- template
  values(r) <- gauss_smooth(template_pts, d$longitude, d$latitude,
                            d$allele_frequency, BW_MAP)
  mask(r, kigali_vect, touches = FALSE)
}

allele_raw      <- rast(setNames(Map(make_raw,    panels$timepoint, panels$locus), panels$label))
allele_smoothed <- rast(setNames(Map(make_smooth, panels$timepoint, panels$locus), panels$label))

allele_raw_long      <- to_long_df(allele_raw)
allele_smoothed_long <- to_long_df(allele_smoothed)


#############################################################
### B3- Allele frequency maps
#############################################################
# admin_sel, not kigali_shp: kigali_shp is the union, so it has only the outer
# edge. admin_sel keeps the districts, which gives the internal borders too.
# Limits from each figure's own values, so the full palette is used. Raw and
# smoothed therefore sit on different scales: raw has one village per pixel and
# reaches the extremes, smoothed averages neighbours and cannot.
# to_long_df <- function(r) {
#   as.data.frame(r, xy = TRUE, na.rm = FALSE) %>%
#     pivot_longer(-c(x, y), names_to = "panel", values_to = "value") %>%
#     filter(!is.na(value))
# }
# 
# plot_allele <- function(long_df, title_text) {
#   ggplot(long_df, aes(x, y, fill = value)) +
#     geom_raster() +
#     geom_sf(data = admin_sel, fill = NA, colour = "grey20", inherit.aes = FALSE) +
#     scale_fill_gradientn(name = "Allele frequency",
#                          colours = brewer.pal(9, "Blues"),
#                          limits = range(long_df$value, na.rm = TRUE),
#                          na.value = NA,
#                          guide = guide_colorbar(barwidth = 18, barheight = 1,
#                                                 title.position = "top")) +
#     facet_wrap(~panel, ncol = length(allele_loci)) +
#     coord_sf(expand = FALSE) +
#     labs(title = title_text) +
#     map_theme()
# }
# 
# p_allele_raw      <- plot_allele(allele_raw_long, "Allele frequency, raw")
# p_allele_smoothed <- plot_allele(allele_smoothed_long,
#                                  paste0("Allele frequency, smoothed, bw = ", BW_MAP))
# 
# print(p_allele_raw)
# print(p_allele_smoothed)
to_long_df <- function(r) {
  as.data.frame(r, xy = TRUE, na.rm = FALSE) %>%
    pivot_longer(-c(x, y), names_to = "panel", values_to = "value") %>%
    filter(!is.na(value))
}

allele_raw_long      <- to_long_df(allele_raw)
allele_smoothed_long <- to_long_df(allele_smoothed)

sum(is.na(allele_raw_long$value))
sum(is.na(allele_smoothed_long$value))

plot_allele <- function(long_df, title_text) {
  ggplot(long_df, aes(x, y, fill = value)) +
    geom_raster() +
    geom_sf(data = admin_sel, fill = NA, colour = "grey20", inherit.aes = FALSE) +
    scale_fill_gradientn(name = "Allele frequency",
                         colours = brewer.pal(9, "Blues"),
                         limits = range(long_df$value, na.rm = TRUE),
                         na.value = NA,
                         guide = guide_colorbar(barwidth = 18, barheight = 1,
                                                title.position = "top")) +
    facet_wrap(~panel, ncol = length(allele_loci)) +
    coord_sf(expand = FALSE) +
    labs(title = title_text) +
    map_theme()
}

p_allele_raw      <- plot_allele(allele_raw_long, "Allele frequency, raw")
p_allele_smoothed <- plot_allele(allele_smoothed_long,
                                 paste0("Allele frequency, smoothed, bw = ", BW_MAP))

print(p_allele_raw)
print(p_allele_smoothed)

ggsave("dataoutput/kigali_allele_raw.png",      p_allele_raw,      width = 10, height = 7, dpi = 300)
ggsave("dataoutput/kigali_allele_smoothed.png", p_allele_smoothed, width = 10, height = 7, dpi = 300)
#############################################################
### B4- Change from Generation 1 to Generation 2
#############################################################
a_gen1 <- grep("^Gen1_", names(allele_smoothed), value = TRUE)
a_gen2 <- grep("^Gen2_", names(allele_smoothed), value = TRUE)

allele_diff_long <- to_long_df(rast(setNames(
  Map(function(g2, g1) allele_smoothed[[g2]] - allele_smoothed[[g1]], a_gen2, a_gen1),
  sub("^Gen1_", "", a_gen1))))

# Purples here, Blues above, so a change map cannot be mistaken for a frequency
# map. Light to dark is honest only while every value is positive: check first.
print(range(allele_diff_long$value, na.rm = TRUE))

p_allele_diff <- ggplot(allele_diff_long, aes(x, y, fill = value)) +
  geom_raster() +
  geom_sf(data = admin_sel, fill = NA, colour = "grey20", inherit.aes = FALSE) +
  scale_fill_gradientn(name = "Change (Gen 2 - Gen 1)",
                       colours = brewer.pal(9, "Purples"), na.value = NA,
                       guide = guide_colorbar(barwidth = 18, barheight = 1,
                                              title.position = "top")) +
  facet_wrap(~panel, ncol = length(allele_loci)) +
  coord_sf(expand = FALSE) +
  labs(title = "Change from Generation 1 to Generation 2") +
  map_theme()

print(p_allele_diff)
ggsave("dataoutput/kigali_allele_change.png", p_allele_diff,
       width = 10, height = 6, dpi = 300)


#############################################################
### B5- Selection pressure
#############################################################
# s depends only on betamat and the covariates, which have no time index, so it
# is the same in both generations: one panel per locus.

setdiff(paste(s_tbl$longitude, s_tbl$latitude),
        paste(round(cov_tbl$longitude, 4), round(cov_tbl$latitude, 4)))
s_tbl <- read.csv("C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/dataoutput/selection_pressure.csv") %>%
  mutate(marker = marker_labels[locus],
         longitude = round(longitude, 4), latitude = round(latitude, 4)) %>%
  pivot_wider(id_cols = c(longitude, latitude), names_from = marker, values_from = s)

s_cols <- unname(marker_labels[allele_loci])
stopifnot(all(s_cols %in% names(s_tbl)))

make_s <- function(col_name) {
  r <- template
  values(r) <- gauss_smooth(template_pts, s_tbl$longitude, s_tbl$latitude,
                            s_tbl[[col_name]], BW_MAP)
  mask(r, kigali_vect, touches = TRUE)
}


s_long <- to_long_df(rast(setNames(lapply(s_cols, make_s), s_cols)))
print(range(s_long$value, na.rm = TRUE))

p_selection <- ggplot(s_long, aes(x, y, fill = value)) +
  geom_raster() +
  geom_sf(data = admin_sel, fill = NA, colour = "grey20", inherit.aes = FALSE) +
  # scale_fill_gradientn(name = "Selection pressure s",
  #                      colours = brewer.pal(9, "Reds"), na.value = NA,
  #                      guide = guide_colorbar(barwidth = 18, barheight = 0.7,
  #                                             title.position = "top")) +
  scale_fill_gradientn(name = "Selection pressure s",
                       colours = brewer.pal(9, "Reds"),
                       trans = "log10",
                       na.value = NA,
                       labels = scales::label_number(accuracy = 0.01),
                       guide = guide_colorbar(barwidth = 18, barheight = 1,
                                              title.position = "top")) +
  facet_wrap(~panel, ncol = length(allele_loci)) +
  coord_sf(expand = FALSE) +
  labs(title = "Selection pressure by locus") +
  map_theme()

print(p_selection)
ggsave("dataoutput/selection.png", p_selection,
       width = 10, height = 6, dpi = 300)

p_selection <- ggplot(s_long, aes(x, y, fill = value)) +
  geom_raster() +
  geom_sf(data = admin_sel, fill = NA, colour = "grey20", inherit.aes = FALSE) +
  
  scale_fill_gradient(name = "Selection pressure s",
                      low = "white", high = "#A5675C",
                      trans = "log2",
                      breaks = c(0.1, 0.25, 0.5, 1, 1.5),
                      limits = range(s_long$value, na.rm = TRUE),
                      na.value = NA,
                      labels = scales::label_number(accuracy = 0.01),
                      guide = guide_colorbar(barwidth = 18, barheight = 1,
                                             title.position = "top")) +
  # scale_fill_gradient2(name = "Selection pressure s",
  #                      low = "white", high = "#67000D",
  #                      # low = "#FFF5F0", mid = "#FB6A4A", high = "brown",
  #                      midpoint = log10(median(s_long$value, na.rm = TRUE)),
  #                      trans = "log10",
  #                      na.value = NA,
  #                      labels = scales::label_number(accuracy = 0.01),
  #                      guide = guide_colorbar(barwidth = 18, barheight = 0.8,
  #                                             title.position = "top")) +
  facet_wrap(~panel, ncol = length(allele_loci)) +
  coord_sf(expand = FALSE) +
  labs(title = "Selection pressure by locus") +
  map_theme()

print(p_selection)


p_selection <- ggplot(s_long, aes(x, y, fill = value)) +
  geom_raster() +
  geom_sf(data = admin_sel, fill = NA, colour = "grey20", inherit.aes = FALSE) +
  scale_fill_gradient(name = "Selection pressure s",
                      low = "white", high = "#67000D",
                      limits = range(s_long$value, na.rm = TRUE),
                      na.value = NA,
                      labels = scales::label_number(accuracy = 0.03 ),
                      guide = guide_colorbar(barwidth = 18, barheight = 1,
                                             title.position = "top")) +
  
  facet_wrap(~panel, ncol = length(allele_loci)) +
  coord_sf(expand = TRUE) +
  labs(title = "Selection pressure by locus") +
  map_theme()

print(p_selection)
ggsave("dataoutput/kigali_selection.png", p_selection, width = 10, height = 6, dpi = 300)


p_selection <- ggplot(s_long, aes(x, y, fill = value)) +
  geom_raster() +
  
  geom_sf(
    data = admin_sel,
    fill = NA,
    colour = "grey20",
    inherit.aes = FALSE
  ) +
  
  scale_fill_gradientn(
    name = "Selection pressure s",
    colours = c(
      "white",
      "#FFF7BC",
      "#FEC44F",
      "#FC4E2A",
      "#E31A1C",
      "#800026"
    ),
    limits = range(s_long$value, na.rm = TRUE),
    na.value = NA,
    labels = scales::label_number(accuracy = 0.03),
    guide = guide_colorbar(
      barwidth = 18,
      barheight = 1,
      title.position = "top"
    )
  ) +
  
  facet_wrap(
    ~panel,
    ncol = length(allele_loci)
  ) +
  
  coord_sf(expand = TRUE) +
  
  labs(
    title = "Selection pressure by locus"
  ) +
  
  map_theme()

print(p_selection)

ggsave(
  "dataoutput/kigali_selection.png",
  p_selection,
  width = 10,
  height = 6,
  dpi = 300
)
#############################################################
### B6- Does the covariate drive the change?
#############################################################
# The change should rise with coverage. Not a clean line: it is proportional to
# s * p * (1 - p), so a pixel starting near 0 or 1 barely moves however strong
# the selection. Hence the colour by starting frequency.
change_vs_cov <- allele_freq_data %>%
  dplyr::select(locus, timepoint, allele_frequency, longitude, latitude) %>%
  pivot_wider(names_from = timepoint, values_from = allele_frequency,
              names_prefix = "gen") %>%
  mutate(delta = gen2 - gen1) %>%
  left_join(cov_tbl %>%
              mutate(longitude = round(longitude, 4), latitude = round(latitude, 4)) %>%
              dplyr::select(longitude, latitude, itn),
            by = c("longitude", "latitude"))

# if this fails, data_to_use.csv came from a different grid than the covariates
stopifnot(!any(is.na(change_vs_cov$itn)))

p_cov_check <- ggplot(change_vs_cov, aes(itn, delta, colour = gen1)) +
  geom_point(size = 2) +
  scale_colour_gradientn(colours = brewer.pal(9, "Blues"),
                         name = "Starting\nfrequency") +
  facet_wrap(~locus) +
  labs(x = "ITN coverage", y = "Change in allele frequency",
       title = "Does the covariate drive the change?") +
  theme_minimal(base_size = 13)

print(p_cov_check)
ggsave("dataoutput/kigali_change_vs_covariate.png", p_cov_check,
       width = 9, height = 5, dpi = 300)



p_selection <- Reduce(`|`, lapply(s_cols, one_marker)) +
  plot_annotation(
    title = "Selection pressure by locus",
    subtitle = "Each panel uses an independent colour scale to show the spatial gradient within each locus.",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "grey35")
    )
  )

print(p_selection)







#############################################################
### B5a- Selection pressure, raw
#############################################################
make_s_raw <- function(col_name) {
  r <- base_rast
  values(r) <- NA
  r[cellFromXY(r, as.matrix(s_tbl[, c("longitude", "latitude")]))] <- s_tbl[[col_name]]
  mask(r, kigali_vect, touches = FALSE)
}

s_raw_long <- to_long_df(rast(setNames(lapply(s_cols, make_s_raw), s_cols)))
print(range(s_raw_long$value, na.rm = TRUE))

plot_selection <- function(long_df, title_text) {
  ggplot(long_df, aes(x, y, fill = value)) +
    geom_raster() +
    geom_sf(data = admin_sel, fill = NA, colour = "grey20", inherit.aes = FALSE) +
    scale_fill_gradientn(name = "Selection pressure s",
                         colours = brewer.pal(9, "Reds"),
                         limits = range(long_df$value, na.rm = TRUE),
                         na.value = NA,
                         labels = scales::label_number(accuracy = 0.001),
                         guide = guide_colorbar(barwidth = 18, barheight = 1,
                                                title.position = "top")) +
    facet_wrap(~panel, ncol = length(allele_loci)) +
    coord_sf(expand = FALSE) +
    labs(title = title_text) +
    map_theme()
}

p_selection_raw      <- plot_selection(s_raw_long, "Selection pressure, raw")
p_selection_smoothed <- plot_selection(s_long, paste0("Selection pressure, smoothed, bw = ", BW_MAP))

print(p_selection_raw)
print(p_selection_smoothed)

ggsave("dataoutput/kigali_selection_raw.png",      p_selection_raw,      width = 10, height = 6, dpi = 300)
ggsave("dataoutput/kigali_selection_smoothed.png", p_selection_smoothed, width = 10, height = 6, dpi = 300)

