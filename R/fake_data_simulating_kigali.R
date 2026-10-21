#############################################################
### 0- Packages and settings
#############################################################
source("R/Packages.R")               # greta, greta.dynamics, terra, tidyverse
source("R/create_dummy_matrices.R")  # L and R matrices
source("R/function_geno_pheno.R")    # probability_genotype_rows_greta()
set.seed(123)

locus_names <- c("marker_1", "marker_2")   # marker_1 = vgsc 995, marker_2 = ace1 280
n_loci      <- 2L

# the years: 4 years = 4 states of the allele frequency
years   <- 2018:2021
n_trans <- length(years) - 1L      # 3 generations: 2018->2019->2020->2021
Tmax    <- n_trans + 1L            # 4 states

# number of mosquitoes tested in each row
n_tested_per_row <- 10L

#############################################################
### 1- The design: where, when and how many
#############################################################
# 1a) the cells = pixels of the mask inside Kigali
mask   <- rast("C:/Users/Sylviane/Desktop/training_perth_2024/ir_cube/data/clean/raster_mask.tif")
kigali <- vect("data_raw/kigali_shp.shp")
kigali <- project(kigali, crs(mask))
kigali_raster <- crop(mask, kigali, mask = TRUE)

# one row per cell with its coordinates
# na.rm = TRUE removes the pixels outside the Kigali boundary
cells <- as.data.frame(kigali_raster, cells = TRUE, xy = TRUE, na.rm = TRUE) %>%
  select(cell, longitude = x, latitude = y)

unique_cells   <- cells$cell
n_cells        <- length(unique_cells)
n_unique_cells <- n_cells
n_cells

# 1b) like in the VA, for marker_1 some cells only report the FREQUENCY
# (no genotype counts), here we choose 10 of them at random
freq_cells <- sample(unique_cells, 10)

# 1c) the design: one row = one cell x one year x one locus
design <- expand_grid(
  cell  = unique_cells,
  year  = years,
  locus = locus_names
) %>%
  left_join(cells, by = "cell") %>%
  mutate(
    n_tested  = n_tested_per_row,
    data_type = if_else(locus == "marker_1" & cell %in% freq_cells,
                        "frequency", "genotype")
  )

# check
#table(design$locus, design$data_type)

######COVARIATES#############################
# extract the values at our cells, long on year, wide on the variable
# the covariates of year t move p from t to t+1, so we need 2018, 2019, 2020
nets_cube <- rast("data_raw/raster_clean/net_use_cube.tif")
irs_cube  <- rast("data_raw/raster_clean/irs_coverage_scaled_cube.tif")
pop_cube  <- rast("data_raw/raster_clean/pop_scaled_cube.tif")

# irs and pop do not reach the last years, so we repeat the last year
# irs_cube <- post_pad_cube(irs_cube, max(years))
# pop_cube <- post_pad_cube(pop_cube, max(years))

# put the cubes on the Kigali grid (same as the main script)
align_to_grid <- function(x, grid) {
  if (!identical(crs(x), crs(grid))) x <- project(x, crs(grid))
  resample(crop(x, ext(grid)), grid, method = "bilinear")
}
nets_k <- align_to_grid(nets_cube, kigali_raster)
irs_k  <- align_to_grid(irs_cube,  kigali_raster)
pop_k  <- align_to_grid(pop_cube,  kigali_raster)

# extract the values at our cells, long on year, wide on the variable
# the covariates of year t move p from t to t+1, so we need 2018, 2019, 2020
all_extract <- bind_cols(
  terra::extract(nets_k, unique_cells),
  terra::extract(irs_k,  unique_cells),
  terra::extract(pop_k,  unique_cells)
) %>%
  mutate(cell = unique_cells, .before = everything()) %>%
  pivot_longer(cols = -cell, names_sep = "_",
               names_to = c("variable", "year"), values_to = "value") %>%
  pivot_wider(names_from = "variable", values_from = "value") %>%
  mutate(year = as.numeric(year)) %>%
  filter(year %in% years[-Tmax]) %>%
  mutate(cell_id = match(cell, unique_cells)) %>%
  arrange(cell_id, year)                 # cell first, year inside: my model needs this order
# covariate matrix: (n_cells * n_trans) x K
x_cell_years <- all_extract %>% select(nets, irs, pop) %>% as.matrix()
dim(x_cell_years)
# number of covariates
K <- ncol(x_cell_years)

stopifnot(nrow(x_cell_years) == n_cells * n_trans,
          !anyNA(x_cell_years))

# loading the model
# possible combination of the genotype
mats <- create_dummy_matrices(n_loci)
L    <- mats$left
R    <- mats$right
SS_mask <- (L == 1L) & (R == 1L)
RR_mask <- (L == 0L) & (R == 0L)
SR_mask <- (L != R)
storage.mode(SS_mask) <- "double"
storage.mode(RR_mask) <- "double"
storage.mode(SR_mask) <- "double"
stopifnot(nrow(SS_mask) == 9)

dim(L)
# 3b) the priors
h          <- uniform(0, 1, dim = n_loci)                 # dominance of R, one per locus
dim(h)
rho_z      <- uniform(0.01, 0.6)                          # overdispersion of genotype counts
alpha0_z   <- (1 - rho_z^2) / rho_z^2
phi_allele <- lognormal(log(20), 0.5)                     # overdispersion of frequency data
p_cells    <- uniform(0, 1, dim = c(n_unique_cells, n_loci))  # p in 2018
# effect of covariates on selection (mean = -2, overdispersion = 0.7)
# exp(-2) = 0.13, we are meaning that the valeu is centered around -2
betamat    <- normal(-2, 0.7, dim = c(K, n_loci))   

## 3c) selection and fitness for each cell-year
# exp(betamat) first so the effects are always positive, then multiply
effect_type          <- exp(betamat)
selection_cell_years <- x_cell_years %*% effect_type      # selection pressure s
fitness_cell_years   <- 1 + selection_cell_years          # relative fitness w of RR

# 3d) fitness of each of the 9 genotypes (log trick: product = sum of logs)
# RR gets w, SR gets h*w + (1-h), SS gets 1, then multiply across loci
ones_G <- ones(nrow(SS_mask), nrow(SS_mask))
h_mat  <- ones(n_cells * n_trans, 1) %*% t(h)

resistant_probability<- log(fitness_cell_years) %*% as_data(t(RR_mask)) 

heterozygous <- log(h_mat * fitness_cell_years + (1 - h_mat)) %*% as_data(t(SR_mask))

Gw <- resistant_probability + heterozygous

r_vec <- exp(Gw)
dim(r_vec)
 dim(r_vec) <- c(n_trans, n_cells, nrow(SS_mask), 1)       # time first for greta.dynamics

 
 # 3e) one generation: HWE -> selection -> back to allele frequency
 genotype_next <- function(state, iter, r_vec, ones_G) {
   p <- state
   dim(p) <- c(nrow(state), ncol(state))     # n_loci x n_cells
   p <- t(p)                                 # n_cells x n_loci
   r <- r_vec
   dim(r) <- c(nrow(p), ncol(ones_G))        # n_cells x G
   z <- probability_genotype_rows_greta(p)   # HWE genotype frequencies
   genotype_post <- z * r                    # selection
   genotype_post <- genotype_post / (genotype_post %*% ones_G)
   p_next <- 0.5 * (genotype_post %*% as_data(SR_mask + 2 * RR_mask))
   out <- t(p_next)
   dim(out) <- dim(state)
   out
 }

 
 
 # 3f) run the dynamics from 2018 for 3 generations
 initial_state <- t(p_cells)
 dim(initial_state) <- c(n_loci, n_cells, 1)
 
 dynamics <- iterate_dynamic_function(
   transition_function       = genotype_next,
   initial_state             = initial_state,
   niter                     = n_trans,
   tol                       = -1,                # never stop early
   r_vec                     = r_vec,
   ones_G                    = ones_G,
   parameter_is_time_varying = "r_vec",
   state_limits              = c(1e-6, 1 - 1e-6)
 ) 

 
 #############################################################
 ### All the years in one matrix (KEEP, this is yours)
 #############################################################
 get_year <- function(k) {
   p <- dynamics$all_states[, , k]     # loci x cells for jump k
   dim(p) <- c(n_loci, n_cells)
   t(p)                                # cells x loci
 }
 
 allele_freq_matrix <- do.call(rbind, c(list(p_cells), lapply(seq_len(n_trans), get_year)))
 # rows = all cells of 2018, then all cells of 2019, ... (YEAR FIRST)
 
 #############################################################
 ### My likelihood as random variables (no data yet)
 #############################################################
 # the population: if p = 0.2 then SS = 0.64, RS = 0.32, RR = 0.04 (HWE)
 Z_matrix     <- probability_genotype_rows_greta(allele_freq_matrix)   # 9 genotypes
 alpha_matrix <- Z_matrix * alpha0_z + 1e-8
 
 n_rows      <- n_cells * Tmax
 geno_counts <- dirichlet_multinomial(size  = rep(n_tested_per_row, n_rows),
                                      alpha = alpha_matrix)
 
 allele_size   <- matrix(2 * n_tested_per_row, n_rows, n_loci)
 allele_counts <- beta_binomial(size  = allele_size,
                                alpha = allele_freq_matrix * phi_allele,
                                beta  = (1 - allele_freq_matrix) * phi_allele)
 
 #############################################################
 ### calculate(): ONE draw = THE TRUTH + THE FAKE COUNTS
 #############################################################
 sim <- calculate(betamat, h, p_cells, rho_z, phi_allele,
                  allele_freq_matrix, Z_matrix, geno_counts, allele_counts,
                  nsim = 1)
 
 # the truth (taken AFTER calculate, from the same draw)
 true_betamat    <- sim$betamat[1, , ]
 true_h          <- sim$h[1, , 1]
 true_p_cells    <- sim$p_cells[1, , ]
 true_rho_z      <- sim$rho_z[1, 1, 1]
 true_phi_allele <- sim$phi_allele[1, 1, 1]
 true_p          <- sim$allele_freq_matrix[1, , ]   # (n_cells * Tmax) x 2
 true_Z          <- sim$Z_matrix[1, , ]             # (n_cells * Tmax) x 9
 
 dimnames(true_betamat) <- list(colnames(x_cell_years), locus_names)
 true_betamat
 true_h
 
 # the fake counts
 fake_geno   <- sim$geno_counts[1, , ]     # (n_cells * Tmax) x 9
 fake_allele <- sim$allele_counts[1, , ]   # (n_cells * Tmax) x 2
 
 #############################################################
 ### Tables: one row per year x cell x locus
 #############################################################
 # the rows of the matrix: YEAR FIRST, cell inside
 rows <- expand_grid(year = years, cell = unique_cells)
 
 # true p
 p_true <- bind_rows(lapply(seq_along(locus_names), function(l) {
   rows %>% mutate(locus = locus_names[l], p_true = true_p[, l])
 }))
 
 # true genotype frequencies per locus, from the model's Z_matrix
 # (e.g. z_SS at marker_1 = z_SS_SS + z_SS_SR + z_SS_RR)
 z_true <- bind_rows(lapply(seq_along(locus_names), function(l) {
   rows %>% mutate(locus = locus_names[l],
                   z_SS  = as.vector(true_Z %*% SS_mask[, l]),
                   z_RS  = as.vector(true_Z %*% SR_mask[, l]),
                   z_RR  = as.vector(true_Z %*% RR_mask[, l]))
 }))
 
 # fake counts per locus (same trick with the masks)
 per_locus <- bind_rows(lapply(seq_along(locus_names), function(l) {
   rows %>% mutate(locus      = locus_names[l],
                   n_SS       = as.vector(fake_geno %*% SS_mask[, l]),
                   n_RS       = as.vector(fake_geno %*% SR_mask[, l]),
                   n_RR       = as.vector(fake_geno %*% RR_mask[, l]),
                   n_R_allele = fake_allele[, l])
 }))
 
 # check: z_true is the HWE of p_true, and each row sums to 1
 z_check <- z_true %>% left_join(p_true, by = c("year", "cell", "locus"))
 stopifnot(sum(is.na(z_check$p_true)) == 0,
           all(abs(z_check$z_SS - (1 - z_check$p_true)^2) < 1e-4),
           all(abs(z_check$z_SS + z_check$z_RS + z_check$z_RR - 1) < 1e-6))
 
 #############################################################
 ### Fill the design: truth + counts
 #############################################################
 # back to the empty design (only the original columns)
 design <- design %>%
   select(cell, year, locus, longitude, latitude, n_tested, data_type)
 
 names(design)   # should be: cell year locus longitude latitude n_tested data_type
 unique(design$locus);  unique(p_true$locus)    # same names?
 range(design$year);    range(p_true$year)      # same years?
 head(design$cell);     head(p_true$cell)       # same cells?
 
 # how many design rows find no match in p_true (should be 0)
 p_true    <- p_true    %>% mutate(locus = tolower(locus))
 z_true    <- z_true    %>% mutate(locus = tolower(locus))
 per_locus <- per_locus %>% mutate(locus = tolower(locus))
 
 unique(p_true$locus); unique(z_true$locus); unique(per_locus$locus)   # all "marker_1" "marker_2"
 
 design_full <- design %>%
   left_join(p_true,    by = c("cell", "year", "locus")) %>%
   left_join(z_true,    by = c("cell", "year", "locus")) %>%
   left_join(per_locus, by = c("cell", "year", "locus")) %>%
   mutate(across(c(n_SS, n_RS, n_RR), ~ if_else(data_type == "frequency", NA, .x)),
          n_R_allele = if_else(data_type == "genotype", NA, n_R_allele))
 
 stopifnot(!anyNA(design_full$p_true), !anyNA(design_full$z_SS))
 stopifnot(!anyNA(design_full$p_true), !anyNA(design_full$z_SS))
 head(design_full)
 names(design)      # if you see n_SS, p_true, n_SS.x ... -> design is not empty anymore
 names(per_locus)   # should contain n_SS, n_RS, n_RR, n_R_allele
 #############################################################
 ### The true parameters and selection
 #############################################################
 true_params <- tibble(
   parameter = c(paste0("betamat[", rep(1:K, n_loci), ",", rep(1:n_loci, each = K), "]"),
                 paste0("h[", 1:n_loci, ",1]"),
                 "rho_z", "phi_allele"),
   covariate = c(rep(colnames(x_cell_years), n_loci), rep(NA, n_loci), NA, NA),
   locus     = c(rep(locus_names, each = K), locus_names, NA, NA),
   true      = c(as.vector(true_betamat), true_h, true_rho_z, true_phi_allele)
 )
 true_params
 
 true_s <- all_extract %>% select(cell, year)
 for (l in seq_along(locus_names)) {
   true_s[[paste0("s_", locus_names[l])]] <- as.vector(x_cell_years %*% exp(true_betamat[, l]))
 }
 
 #############################################################
 ### Save: the fake data and the truth SEPARATELY
 #############################################################
 dir.create("dataoutput", showWarnings = FALSE)
 
 # full table (truth + counts): for you to check, never given to the model
 write.csv(design_full, "dataoutput/design_with_truth_kigali_2loci.csv", row.names = FALSE)
 fake_data <- design_full  
 # the FAKE DATA: only what a field team would see -> the fit reads this
 write.csv(fake_data, "dataoutput/fake_data_kigali_2loci.csv", row.names = FALSE)
 names(fake_data)
 # the TRUTH
 write.csv(true_params, "dataoutput/truth_params_kigali_2loci.csv",    row.names = FALSE)
 write.csv(p_true,      "dataoutput/truth_p_kigali_2loci.csv",         row.names = FALSE)
 write.csv(z_true,      "dataoutput/truth_z_kigali_2loci.csv",         row.names = FALSE)
 write.csv(true_s,      "dataoutput/truth_selection_kigali_2loci.csv", row.names = FALSE)
 saveRDS(list(params = true_params, p_true = p_true, z_true = z_true,
              p_matrix = true_p, selection = true_s),
         "dataoutput/truth_kigali_2loci.rds")
 