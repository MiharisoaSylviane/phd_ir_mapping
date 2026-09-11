# ###############################################################
# This script is executing those following in the order:
#   1) call the fix data like dummy matrices data L/R, villages, and number of loci
#   2) Preparing the covariates and read the data
#   3) Preparing the data
#   3) calling the function of the phenotype of each genotype
#   4) Fake data generating with the prior and likelihood
#   5) Plotting the DAG to see the nodes, plot the parameters by using mcmc_draws
#   6) Estimate the parameters (betamat, h, rho-z and p) by using mcmc package and the fake data
#   7) Interpreting the results of the mcmc 
#   8) 
#
# Loading packages
source("R/Packages.R")
set.seed(123)
#############################################################
### 1- Possible genotype combination
#############################################################
# Dummy matrices 
source("R/create_dummy_matrices.R")
source("R/function_geno_pheno.R")

# number of loci
n_loci <- 2L

# creation of the matrix genotype in row and loci in the column
mats   <- create_dummy_matrices(n_loci)
L      <- mats$left
R      <- mats$right
# G      <- nrow(L)
# B      <- ncol(L) # B = number of loci 

# check the dimension of the object L
dim(L)

# naming the locus for the sake of the dataframe later
locus_names <- c("Marker_1", "Marker_2")
colnames(L) <- locus_names
colnames(R) <- locus_names

# Locus transformation to match the Greta syntax
# version greta as the matrix creation won't be identified by Greta
# because Greta can't keep the matrix with qualitative data part of the plain R we have then 
# here we are transforming the L==1L to boleen TRUE OR FALSE
genotype_label_df <- map_dfc(locus_names, function(loc) {
  tibble(!!loc := case_when(
    L[, loc] == 1 & R[, loc] == 1 ~ "SS",
    L[, loc] == 0 & R[, loc] == 0 ~ "RR",
    TRUE                          ~ "SR"
  ))
})

genotype_lookup <- genotype_label_df %>%
  mutate(genotype_id = row_number(), .before = 1)
### defining the combination SS
SS_mask <- (L == 1L) & (R == 1L)
RR_mask <- (L == 0L) & (R == 0L)
SR_mask <- (L != R)
# here we are saving them as 0 or 1
storage.mode(SS_mask) <- "double"

storage.mode(RR_mask) <- "double"
storage.mode(SR_mask) <- "double"
dup <- 1 + L - R

# defining those as variables to be recognized by Greta
# as_data() is saying, take the value of those to be the data used in our model
L_g       <- as_data(L) 
R_g       <- as_data(R)
SS_mask_g <- as_data(SS_mask) 
RR_mask_g <- as_data(RR_mask)
SR_mask_g <- as_data(SR_mask)
# number of generations
n_trans <- 30L
Tmax    <- n_trans + 1L            

# number of mosquitoes tested
M_z <- 10

###########################################################################
## 2 - Preparing the covariates
#############################################################
# set the start of the timeseries considered in modelling (the start of
# non-negligible levels of resistance) - assume it's before the mass-rollout of
# nets
baseline_year <- 1995

# set the final year of data (insufficient and spatially biased data for 2025)
final_data_year <- 2024

nets_cube <- rast("data_raw/raster_clean/net_use_cube.tif")

# IRS coverage
irs_cube <- rast("data_raw/raster_clean/irs_coverage_scaled_cube.tif")

# human population
pop_cube <- rast("data_raw/raster_clean/pop_scaled_cube.tif")

# for each one pad back to the baseline year, repeating the first
nets_cube <- pre_pad_cube(nets_cube, baseline_year)
irs_cube <- pre_pad_cube(irs_cube, baseline_year)
pop_cube <- pre_pad_cube(pop_cube, baseline_year)

# if necessary post-pad
# if necessary, pad forward to the final data year, repeating the last
nets_cube <- post_pad_cube(nets_cube, final_data_year)
irs_cube <- post_pad_cube(irs_cube, final_data_year)
pop_cube <- post_pad_cube(pop_cube, final_data_year)


# load the non-temporal crop covariate layers

# collated total yields of crop types
crops_group <- rast("data_raw/raster_clean/crop_group_scaled.tif")

# yields of individual crops
crops_all <- rast("data_raw/raster_clean/crop_scaled.tif")

# Pull out crop types implicated in risk for IR. refer to this review, crop type
# section:
# https://malariajournal.biomedcentral.com/articles/10.1186/s12936-016-1162-4
crops_implicated <- c(
  # "increased resistance at cotton growing sites, a finding subsequently
  # supported in eight other papers from five different African countries", "the
  # cash crop with the highest intensity insecticide use of any crop"
  crops_all$cotton,
  # "In eight studies, vegetable cultivation strongly related to
  # insecticide-resistant field collections", "Vegetable production requires
  # significantly higher quantities and/or more frequent application of
  # pesticides than other food crops"
  crops_all$vegetables,
  # "Seven of the studies reviewed here examined the insecticide susceptibility
  # of vector populations at rice-growing sites, and found low-to-moderate
  # resistance levels in these mosquito populations."
  crops_all$rice)

# combine all temporally-static covariates
covs_flat <- c(crops_group, crops_implicated)

# read the mask and reduce to kigali
mask <-  rast("C:/Users/Sylviane/Desktop/training_perth_2024/ir_cube/data/clean/raster_mask.tif")
# plot(mask)
kigali <- vect("data_raw/kigali_shp.shp")

# checking the projection
kigali <- project(kigali, crs(mask))
kigali_raster <- crop(mask, kigali, mask = TRUE)
raster_xy <- as.data.frame(
  xyFromCell(kigali_raster, 1:ncell(kigali_raster))
)

# align by checking the projection, the extent and cropping the raster
# resample to match the kigali shapefile itself
align_to_grid <- function(x, grid) {
  if (!identical(crs(x), crs(grid))) x <- project(x, crs(grid))
  resample(crop(x, ext(grid)), grid, method = "bilinear")
}

nets_k  <- align_to_grid(nets_cube, kigali_raster)
irs_k   <- align_to_grid(irs_cube,  kigali_raster)
pop_k   <- align_to_grid(pop_cube,  kigali_raster)
flat_k  <- align_to_grid(covs_flat, kigali_raster)
plot(nets_k)
stopifnot(
  all(dim(nets_k)[1:2] == dim(kigali_raster)[1:2]),
  all(dim(flat_k)[1:2] == dim(kigali_raster)[1:2])
)

###############################################
# 3- Read the data to fit the model
#############################################33
allele_data <- read.csv("dataoutput/allele_frequency_data.csv")
genotype    <- read.csv("dataoutput/data_to_use.csv")

# the grid, and the cells, come from the genotype file since only it has coords
# kigali_raster <- crop(mask_r, kigali, mask = TRUE)
coords <- genotype %>%
  distinct(village_id, longitude, latitude) %>%
  mutate(cell = cellFromXY(mask, as.matrix(cbind(longitude, latitude)))) %>%
  arrange(village_id)

unique_cells <- coords$cell
n_cells      <- nrow(coords)
stopifnot(!any(is.na(unique_cells)), !any(duplicated(unique_cells)))

# Because the genotype data I generated wasn't matching the allele data
# so I needed to arrange it a bit, also making it look like VA data (see obs_wide)

# creating a new colmun to store Marker_1 and 2
geno_levels <- genotype %>%
  distinct(genotype_id, Marker_1, Marker_2) %>%
  arrange(genotype_id) %>%
  mutate(geno_col = paste0("Marker1", Marker_1, "_Marker2", Marker_2))
# the idea here is to get a dataframe with the structure
# marker1_RR_SS and marker2_RR_SS,.... as column, the row value will be the observed
# or the allele frequency
genotype_wide <- genotype %>%
  left_join(geno_levels, by = c("genotype_id", "Marker_1", "Marker_2")) %>%
  mutate(geno_col = factor(geno_col, levels = geno_levels$geno_col)) %>%
  select(-genotype_id, -Marker_1, -Marker_2, -genotype) %>%
  pivot_wider(
    names_from   = geno_col,
    values_from  = c(n_observed,z_true),
    values_fill  = 0,
    names_expand = TRUE
  ) %>%
  arrange(village_id, timepoint)          

# expected dataframe structure is Marker_1_allele and Marker_2_allele
# the row value is the number and the frequency of the resistant allele
allele_wide <- allele_data %>%
  mutate(locus = factor(locus, levels = locus_names)) %>%
  # select(village_id, timepoint, locus, allele_count_observed) %>%
  pivot_wider(
    names_from   = locus,
    values_from  = c(allele_count_observed,allele_frequency),
    values_fill  = 0,
    names_expand = TRUE
  ) %>%
  arrange(village_id, timepoint)

# this is the dataframe that joining both allele and genotype frequency dataframe
obs_wide <- genotype_wide %>%
  left_join(allele_wide, by = c("village_id", "timepoint"),
            suffix = c("_geno", "_allele")) %>%
  arrange(village_id, timepoint)

# this is a check to see if the row data is matching
stopifnot(
  nrow(obs_wide) == nrow(genotype_wide),      
  nrow(obs_wide) == n_cells * length(unique(genotype$timepoint)),
  !any(is.na(obs_wide)))                     

# Matching the dataframe and the raster coordinates
obs_wide <- obs_wide %>%
  mutate(cell = cellFromXY(kigali_raster,
                           as.matrix(cbind(longitude, latitude))))

stopifnot(!any(is.na(obs_wide$cell)))

# create some indices 
unique_cells <- unique(obs_wide$cell)
n_cells      <- length(unique_cells)
years <- baseline_year - 1 + sort(unique(obs_wide$timepoint))

# create a new cell_id 
obs_wide <- obs_wide %>%
  mutate(cell_id = match(cell, unique_cells)) %>%
  arrange(cell_id, timepoint)

stopifnot(
  n_cells == length(unique(obs_wide$village_id)),   # one pixel per village
  !any(duplicated(obs_wide[, c("cell_id", "timepoint")]))
)

#Pulling out the two matrices
geno_cols <- geno_levels$geno_col

observed_counts_geno   <- as.matrix(obs_wide[, paste0("n_observed_", geno_cols)])
z_true_matrix          <- as.matrix(obs_wide[, paste0("z_true_", geno_cols)])
observed_counts_allele <- as.matrix(obs_wide[, paste0("allele_count_observed_", locus_names)])
allele_true_matrix     <- as.matrix(obs_wide[, paste0("allele_frequency_", locus_names)])

storage.mode(observed_counts_geno)   <- "double"
storage.mode(z_true_matrix)          <- "double"
storage.mode(observed_counts_allele) <- "double"
storage.mode(allele_true_matrix)     <- "double"

stopifnot(
  all(dim(observed_counts_geno)   == c(98, nrow(SS_mask))),
  all(dim(observed_counts_allele) == c(98, n_loci)),
  all(rowSums(observed_counts_geno) == M_z),
  all(abs(rowSums(z_true_matrix) - 1) < 1e-8),
  all(allele_true_matrix >= 0 & allele_true_matrix <= 1)
)

obs_rows <- (obs_wide$village_id - 1L) * Tmax + obs_wide$timepoint
stopifnot(!any(duplicated(obs_rows)), max(obs_rows) <= n_cells * Tmax)
# pull out temporally-static covariates for all cells
flat_extract <- covs_flat %>%
  extract(unique_cells) %>%
  mutate(
    cell = unique_cells,
    .before = everything()
  )

# extract spatiotemporal covariates from the cube
all_extract <- bind_cols(
  terra::extract(nets_k, unique_cells),
  terra::extract(irs_k, unique_cells),
  terra::extract(pop_k, unique_cells)
) %>%
  mutate(
    cell = unique_cells,
    .before = everything()
  ) %>%
  # this stacks all the different cubes in long format, but we want wide on the
  # variable but long on year, so pivot_wider immediately after
  pivot_longer(
    cols = -one_of("cell"),
    names_sep = "_",
    names_to = c("variable", "year"),
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = "variable",
    values_from = "value"
  ) %>%
  mutate(
    year = as.numeric(year)
  ) %>%
  left_join(
    flat_extract,
    by = "cell"
  ) %>%
  mutate(
    cell_id = match(cell, unique_cells),
    year_id = year - baseline_year + 1,
    .before = everything()
  ) %>%
  filter(
    year >= baseline_year
  ) %>%
  select(
    -cell,
    -year,
  )
all_extract <- all_extract %>% 
  select( cell_id, 
          year_id,
          nets,
          irs,
          pop)
# pull out index to cells and years
cell_years_index <- all_extract %>%
  select(cell_id, year_id)

# get covariates for these cell-years as a matrix
x_cell_years <- all_extract %>%
  select(-cell_id,
         -year_id) %>%
  as.matrix()
# dimensions of things in the fitting stage
K      <- ncol(x_cell_years)
n_obs          <- nrow(obs_wide)            
n_unique_cells <- length(unique_cells)      
n_times        <- max(obs_wide$timepoint) - min(obs_wide$timepoint) + 1


##################################################################
### 2- Defining Priors
##################################################################
# betamat the effect covariate on the selection, h: dominance, 
# rho-z : overdispersion that will be used in the dirichlet distribution
# p_village : initial allele frequency
betamat  <- normal(0, 5, dim = c(n_loci, K))
# betamat <- calculate(betamat, nsim= 1)
effect_type <- exp(betamat)
#### dominance of the resistant allele in the heterozygous population
#h         <- beta(2, 2, dim = n_loci)  # this is saying that the doominance
# h is spinning around 0.5
h      <- uniform(0, 1, dim = n_loci)

#### overdispersion
rho_z     <- uniform(0.01, 0.6)  # to avoid 0 and extreme value like 100 in n_observed                              
alpha0_z    <- (1 - rho_z^2) / rho_z^2
#### p_village:  initial allele frequency that we will use to estimate the selection effect
# p_village <- uniform(0, 1, dim = c(n_villages, n_loci)) 
p_cells <- uniform(0, 1, dim = c(n_unique_cells, n_loci)) 

### overdispersion for the phenotype/bioassay data
phi_add  <- lognormal(log(20), 0.5) 

### overdisperion for the allele frequency (hyperparameter)
phi_allele <- lognormal(log(20), 0.5)

### theta is the locus-specific, insecticide-specific effect\
# meaning the contribution of a particular genetic locus to resistance against 
#  a particular insecticide
theta <- exponential(0.5, dim = c(n_unique_cells * Tmax, n_loci))
# theta <- uniform(0, 1, dim = c(n_villages * Tmax, n_loci))

# 3- calling the function of the probability of having one genotype and getting one genotype
source("R/function_geno_pheno.R")

# Build a list to put the list of rho_z by village 
# this equation is needed for the Dirichlet distribution
alpha0_z <- (1 - rho_z^2) / rho_z^2


# creating a vector where we could store a list of dimension village X time
# it was equivalent to the step where we are creating vector or matrix for data storage
# alpha_rows  <- vector("list", n_unique_cells * Tmax)
# # so this will be null or empty then we will fill that with the 
# Z_rows      <- vector("list", n_unique_cells * Tmax) # list of genotypes frequency
# allele_freq_rows <- vector("list", n_unique_cells * Tmax) # list of allele frequency by village and time
# s_rows <- vector("list", n_unique_cells) # list of the selection pressure by pixel
# w_rows <- vector("list", n_unique_cells) # list of the relative fitness by pixel
# # village_vec <- integer(n_unique_cells * Tmax)
# time_vec    <- integer(n_unique_cells * Tmax)
# n_unique_cells <- length(unique_cells)

# pdied_mult_rows <- vector("list", n_villages * Tmax)  
### Betamat represent the effect of covariates on the selection
betamat   <- normal(-2, 0.7, dim = c(K, n_loci))

# convert beta to positive effect sizes, them matrix multiply with indicator
# variables to constrain to positive effects of the covariates on selection
# coefficient/relative fitness
effect_type <- exp(betamat)

# # compute selection coefficients for cell-years, convert to relative fitness,
# and extract at data locations
# selection pressure
selection_cell_years <- x_cell_years %*% effect_type
dim(x_cell_years)
# check
# dim(selection_cell_years)
# dim(x_cell_years)
# dim(effect_type)

# relative fitness
fitness_cell_years <- 1 + selection_cell_years

dim(fitness_cell_years)

# reformat this in to a 3D array with dimensions:
#   n_times x n_unique_cells x n_types x 1
# to solve dynamics with time-varying fitness (time must be first, then other
# two must match state variable, which has a trailing dimension of size 1)
fitness_array <- fitness_cell_years
dim(fitness_array) <- c(n_trans, n_cells, n_loci, 1)
dim(fitness_array)
# this is how I used the loop to compute w and s
# w_rows <- vector("list", n_villages)
# 
# for (i in seq_len(n_villages)) {
#   w_rows[[i]] <- t(compute_w_greta(betamat, X_villages[i, ]))
# }
# w_matrix <- do.call(greta::abind, c(w_rows, list(along = 1)))   # n_villages x n_loci
# 
# stopifnot(inherits(w_matrix, "greta_array"))
# starting point 
# dim(SR_mask) # 9 x 2
# dim(p_cells) # 49 x 2 meaning p_cells x SR_mask
# is impossible so that's why we need the transpose and it will give 
# 49 x 2 X 2 x 9
# the log here is needed as we are trying to use the product function in
# in the mosquitoes, product won't be easy to compute
# so we are using the log(a x b) = log(a) + log(b)
# the product is telling biologcally, what is the probability
# of having SSSS at the loci 1 and 2 or RRRR or RSRS or any combination
# calling the function that computing the probability of having 
# a gentoype X
# dimensions should be n_cells x G
Z_init     <- probability_genotype_rows_greta(p_cells)  
# here we transpose to be G x n_cells
state_init <- t(Z_init)                                  

# chk <- calculate(Z_init, nsim = 1)
# # One generation, for greta.dynamics::iterate_dynamic_function
# # state   : n_cells x n_loci x 1 allele frequency
# # iter    : generation number, unused (fitness arrives already sliced)
# # fitness : n_cells x n_loci x 1 relative fitness for this generation
# fit_1 <- fitness_array[, , 1]    # adjust indexing to your array's rank
# dim(initial_state)      # what the state must look like every step
# dim(fitness_array)      # last dim must equal n_times - 1
# n_times - 1 

# diploid_next <- function(state, iter, fitness, h, L_g, R_g,
#                          SS_mask, RR_mask, SR_mask) {
#   
#   n_cells <- dim(state)[1]
#   n_loci  <- dim(state)[2]
#   
#   p <- state
#   dim(p) <- c(n_cells, n_loci)
#   
#   w <- fitness
#   dim(w) <- c(n_cells, n_loci)
#   
#   z      <- probability_genotype_fast_greta(p, L_g, R_g)
#   z_post <- polygenic_multilocus_next_step_greta(z, w, h,
#                                                  SS_mask, RR_mask, SR_mask)
#   p_next <- allele_frequency_from_genotype_greta(z_post, L, R)
#   
#   dim(p_next) <- c(n_cells, n_loci, n_times, 1)
#   p_next
# }

# tol is negative so the iteration never stops early.


# initial state, needs the trailing dimension of 1
initial_state <- t(Z_init)
dim(initial_state) <- c(nrow(SS_mask), n_cells, 1)
# this step will fix the data as a matrix as nrow(SS_mask = 9)
# so we have ones_G 9 X 9
ones_G <- as_data(matrix(1, nrow(SS_mask), nrow(SS_mask)))

# Gw and r_vec, built once, fitness does not vary in time 
# dominance should be a matrix with n_cells * n_loci 
# reminder ones(m,)
# dimension  =  n_cells x n_loci
h_mat <- ones(n_cells * n_trans, 1) %*% t(h)                       
# dimensions for the polygenic score n_cells x G
# this is to estimate the probablity of having one genotype X
Gw <- log(fitness_cell_years) %*% as_data(t(RR_mask)) +
  log(h_mat * fitness_cell_years + (1 - h_mat)) %*% as_data(t(SR_mask))
dim(Gw)
# dimensions should be  G x n_cells
r_vec <- t(exp(Gw))                                     

# checks before iterating
# the idea is that the dim should be nearly equal to allow the multiplication

dim(initial_state)      
dim(r_vec)              
dim(ones_G)     
dim(selection_cell_years)
dim(fitness_cell_years)

# transition, one generation
r_vec <- exp(Gw)
dim(r_vec) <- c(n_trans, n_cells, nrow(SS_mask), 1)      # 30 x 49 x 9 x 1

genotype_next <- function(state, iter, r_vec, ones_G) {
  
  p <- state
  dim(p) <- c(nrow(state), ncol(state))     # n_loci x n_cells
  p <- t(p)                                 # n_cells x n_loci
  
  r <- r_vec
  dim(r) <- c(nrow(p), ncol(ones_G))        # n_cells x G
  
  z <- probability_genotype_rows_greta(p)
  
  genotype_post <- z * r
  genotype_post <- genotype_post / (genotype_post %*% ones_G)
  
  p_next <- 0.5 * (genotype_post %*% as_data(SR_mask + 2 * RR_mask))
  
  out <- t(p_next)
  dim(out) <- dim(state)
  out
}
initial_state <- t(p_cells)
dim(initial_state) <- c(n_loci, n_cells, 1)
dim(initial_state)
class(initial_state)
dynamics <- iterate_dynamic_function(
  transition_function = genotype_next,
  initial_state       = initial_state,
  niter               = n_trans,
  tol                 = -1,
  r_vec               = r_vec,
  ones_G              = ones_G,
  parameter_is_time_varying = "r_vec",
  state_limits        = c(1e-6, 1 - 1e-6)
)

dim(dynamics$all_states)       # expect n_loci x n_cells x n_trans

# check
chk <- calculate(dynamics$all_states, nsim = 1)
st  <- chk$all_states[1, , , , ]
sim_result <- calculate(betamat, h, rho_z, p_cells, selection_cell_years,Z_init, fitness_cell_years, nsim = 1)
# observational model for allele frequency
allele_freq_matrix <- do.call(rbind, c(
  list(p_cells),
  lapply(seq_len(n_trans), function(k) {
    slice <- dynamics$all_states[, , k]
    dim(slice) <- c(n_loci, n_cells)
    t(slice)
  })
))[as.vector(t(matrix(seq_len(n_cells * Tmax), nrow = n_cells, ncol = Tmax))), ]


# likelihood 
Z_matrix     <- probability_genotype_rows_greta(allele_freq_matrix)
alpha_matrix <- Z_matrix * alpha0_z + 1e-8

alpha_matrix <- Z_matrix * alpha0_z + 1e-8
size_vector  <- rep(M_z, n_cells * Tmax)
obs_geno   <- as_data(observed_counts_geno)     # 98 x 9
obs_allele <- as_data(observed_counts_allele)   # 98 x 2

eps <- 1e-6

# genotype counts
distribution(obs_geno) <- dirichlet_multinomial(
  size  = rep(M_z, nrow(observed_counts_geno)),
  alpha = alpha_matrix[obs_rows, ]
)



# allele counts
allele_freq_safe <- eps + (1 - 2 * eps) * allele_freq_matrix[obs_rows, ]

obs_allele <- as_data(observed_counts_allele)


distribution(obs_allele) <- beta_binomial(
  size  = 2 * M_z,
  alpha = allele_freq_safe * phi_allele,
  beta  = (1 - allele_freq_safe) * phi_allele
)

# build the model
m <- model(
  # dominance
  h,
  # effect of covariates
  betamat,
  p_cells
)



# draws 
draws <- greta::mcmc(
  model     = m,
  n_samples = 1000,
  warmup    = 1000,
  chains    = 4
)

# plotting the chains iteration
library(bayesplot)
mcmc_trace(draws, regex_pars = c("h"))

mcmc_trace(draws, regex_pars = "p_cells\\[(1|2|3|4|5|6|7|8|9|10),")

mcmc_trace(draws, regex_pars = c("betamat"))


 # check convergence
rhats <- coda::gelman.diag(draws,
                           autoburnin = FALSE,
                           multivariate = FALSE)
summary(rhats$psrf)


#these have to match the arguments to model(), above, and neet to be greta
# variable nodes (not operation nodes)
posts <- calculate(
  h_mat,
  betamat,
  p_cells,
  values = draws,
  nsim = 100
)
posts <- calculate(allele_freq_matrix, selection_cell_years, values = draws, nsim = 100)
post   <- do.call(rbind, lapply(draws, as.matrix))
p1_hat <- matrix(colMeans(post[, grep("^p_cells", colnames(post))]),
                 n_cells, n_loci)

r <- kigali_raster
values(r) <- NA
r[unique_cells] <- p1_hat[, 1]              # locus 1

plot(r, range = c(0, 1), col = hcl.colors(20, "Blues", rev = TRUE))
plot(kigali, add = TRUE, border = "grey30")

dim(Gw)
dim(r_vec)
dim(h_mat)

calculate(fitness_cell_years, values = draws, nsim = 2)
calculate(h_mat, values = draws, nsim = 2)
calculate(Gw, values = draws, nsim = 2)

Gw <- log(fitness_cell_years) %*% as_data(t(RR_mask)) +
  log(h_mat * fitness_cell_years + (1 - h_mat)) %*% as_data(t(SR_mask))

calculate(Gw, values = draws, nsim = 2)
calculate(log(fitness_cell_years) %*% as_data(t(RR_mask)), values = draws, nsim = 2)
calculate(log(h_mat * fitness_cell_years + (1 - h_mat)) %*% as_data(t(SR_mask)), values = draws, nsim = 2)

traj  <- calculate(allele_freq_matrix, values = draws, nsim = 200)$allele_freq_matrix
p_hat <- apply(traj, c(2, 3), mean)                    # (n_cells * Tmax) x n_loci
p_cube <- array(p_hat, dim = c(Tmax, n_cells, n_loci)) # rows are time-fastest

r <- kigali_raster
values(r) <- NA
r[unique_cells] <- p_cube[Tmax, , 1]        # last generation, locus 1

plot(r, range = c(0, 1), col = hcl.colors(20, "Blues", rev = TRUE))
plot(kigali, add = TRUE, border = "grey30")


dim(betamat)
dim(effect_type)
dim(selection_cell_years)
dim(fitness_cell_years)
dim(Gw)
# # function to simulate
# for(i in seq_len(n_villages)) {
#   
#   w_i <- compute_w_greta(betamat, X_villages[i, ])
#   # transpose is used here because we have (p1,p2,p3) which is 1 X 3
#   # however those should be in the column like
#   # [p1]
#   # [p2]
#   # [p3] which give us (3 X 1)
#   s_i <- compute_s_greta(betamat, X_villages[i, ])
#   s_rows[[i]] <- t(s_i)
#   w_rows[[i]] <- t(w_i)
#     
# 
#   #  
#   Z_list_i <- simulate_genotype_timecourse_greta(
#     p = p_i, w = w_i, h = h, Tmax = Tmax,
#     SS_mask = SS_mask_g, RR_mask = RR_mask_g, SR_mask = SR_mask_g,
#     L = L_g, R = R_g
#   )
#   
#   for (t in seq_len(Tmax)) {
#     # genotype frequency (1 X G)
#     Z_rows[[row_id]]     <- t(Z_list_i[[t]])
#     alpha_rows[[row_id]] <- Z_rows[[row_id]] * alpha0_z + 1e-8   
#     village_vec[row_id]  <- i
#     time_vec[row_id]    <- t
#     theta_it            <- t(theta[row_id, ])
#     Ugc <- compute_Ugc(theta_it, h, SS_mask_g, RR_mask_g, SR_mask_g)
#   
#     # allele
#     allele_freq_rows[[row_id]] <- t(allele_frequency_from_genotype_greta(Z_list_i[[t]], L_g, R_g))
#     
#     U_add_i  <- compute_Ustar_additive(Ugc_i, theta)
#     U_add_i <- compute_Ustar_additive(Ugc, theta_it)
#     # U_mult_i <- compute_Ustar_multiplicative(Ugc, theta_it)
#     
#     pdied_add_rows[[i]]  <- t(compute_p_died(U_add_i))
#     pdied_add_rows[[row_id]] <- t(compute_p_died_additive(U_add_i))
#     # pdied_mult_rows[[row_id]] <- t(compute_p_died(U_mult_i))
#     row_id              <- row_id + 1
#     
#   }
# }
# 
# # (n_villages*Tmax) x G, true frequency
# Z_matrix     <- do.call(greta::abind, c(Z_rows, list(along = 1)))  
# # (n_villages*Tmax) x G
# alpha_matrix <- do.call(greta::abind, c(alpha_rows, list(along = 1)))  
# size_vector  <- rep(M_z, length(alpha_rows))  
# # alleles
# allele_freq_matrix   <- do.call(greta::abind, c(allele_freq_rows, list(along = 1)))
# # n_villages x G
# # p_died_add_matrix  <- do.call(greta::abind, c(pdied_add_rows,  list(along = 1)))
# # p_died_mult_matrix <- do.call(greta::abind, c(pdied_mult_rows, list(along = 1)))
# s_matrix <- do.call(greta::abind, c(s_rows, list(along = 1)))   # n_villages x n_loci
# w_matrix <- do.call(greta::abind, c(w_rows, list(along = 1)))
# 
# class(s_matrix)
# dim(s_matrix)




# to check bug or any errors
# dim(alpha_matrix)        # 70 9
# dim(allele_freq_matrix)  # 70 2
# dim(p_died_add_matrix)   # 70 9
# length(size_vector)      # 70
# 
# all(village_vec == rep(seq_len(n_villages), each = Tmax))
# all(time_vec    == rep(seq_len(Tmax), times = n_villages))
#########################################################
### 4- Fake data generating with the prior and likelihood
#########################################################
#sim_result <- calculate(alpha_matrix, theta, Z_matrix, betamat, h, rho_z, p_village, p_died_mult_matrix, nsim = 1)
# sim_result <- calculate(alpha_matrix, theta, Z_matrix, betamat, h, rho_z, p_died_add_matrix, p_village, phi_add, allele_freq_matrix, phi_allele, s_matrix, w_matrix, nsim = 1)
# #sim_result <- calculate(alpha_matrix, Z_matrix, betamat, h, rho_z, p_village, allele_freq_matrix, s_matrix, w_matrix, nsim = 1)
# # str(sim_result)
# # range(true_p_village)
# # range(true_allele_freq)
# # range(true_Z_matrix)
# # (n_villages*Tmax) x G
# 
# # range(rowSums(true_Z_matrix))
# 
# alpha_numeric  <- sim_result$alpha_matrix[1, , ] 
# # (n_villages*Tmax) x G
# true_Z_matrix  <- sim_result$Z_matrix[1, , ]       
# true_betamat   <- sim_result$betamat[1, , ]
# true_h         <- sim_result$h[1, , ]
# true_rho_z     <- as.numeric(sim_result$rho_z)[1]
# true_p_village <- sim_result$p_village[1, , ]
# 
# true_theta          <- sim_result$theta[1, , ]
# true_p_died_add      <- sim_result$p_died_add_matrix[1, , ]
# # true_p_died_mult      <- sim_result$p_died_mult_matrix[1, , ]
# true_allele_freq <- sim_result$allele_freq_matrix[1, , ]
# true_phi_allele  <- as.numeric(sim_result$phi_allele)[1]
# true_w <- sim_result$w_matrix[1, , ]
# true_s <- sim_result$s_matrix[1, , ]
# 
# # allele sanity check
# if (any(true_allele_freq < 0 | true_allele_freq > 1)) {
#   stop("true_allele_freq out of [0,1]. range = ", paste(round(range(true_allele_freq), 6), collapse = " to "))
# }
# # Gneotype count
# # Likelihood, a here is one row of the alpha_numeric at a time
# fake_counts_matrix <- t(apply(alpha_numeric, 1, function(a) {
#   z_disp <- as.numeric(MCMCpack::rdirichlet(1, a))
#   as.vector(rmultinom(1, M_z, z_disp / sum(z_disp)))
#   
# }))
# 
# n_tested_vec <- rowSums(fake_counts_matrix) 
# 
# # Phenotype count
# 
# # fake_dead_add  <- matrix(rbinom(length(true_p_died_add),  M_z, true_p_died_add),
# #                          nrow = n_villages * Tmax)
# 
# # fake_dead_mult <- matrix(rbinom(length(true_p_died_mult), M_z, true_p_died_mult),
# #                          nrow = n_villages * Tmax)
# 
# 
# # tibble of the fake data we tried before: row =  village x génotype x timepoint
# ## the option (scipen =999) allow us to avoid the e-10 that would make the dataframe
# # strange
# options(scipen = 999)
# # mcmc_data_all_sim <- map_dfr(seq_len(nrow(fake_counts_matrix)), function(r) {
# #   genotype_lookup %>%
# #     mutate(
# #       n_observed = fake_counts_matrix[r, ],
# #       z_true     = true_Z_matrix[r, ],
# #       n_tested   = n_tested_vec[r],
# #       timepoint  = time_vec[r],
# #       village    = villages$village[village_vec[r]],
# #       latitude   = villages$latitude[village_vec[r]],
# #       longitude  = villages$longitude[village_vec[r]],
# #       # p_died_add_true     = true_p_died_add[r, ],     
# #       # p_died_mult_true    = true_p_died_mult[r, ],     
# #       # dead_add_observed   = fake_dead_add[r, ],        
# #       # dead_mult_observed  = fake_dead_mult[r, ] 
# #     )
# # })
# # mcmc_data_all_sim <- mcmc_data_all_sim %>%
# #   mutate(across(where(is.numeric), ~ format(., scientific = FALSE)))
# # locus_cols <- setdiff(names(genotype_lookup), "genotype_id")
# # mcmc_data_all_sim <- map_dfr(seq_len(nrow(fake_counts_matrix)), function(r) {
# #   genotype_lookup %>%
# #     mutate(
# #       n_observed = fake_counts_matrix[r, ],
# #       z_true     = true_Z_matrix[r, ],
# #       n_tested   = n_tested_vec[r],
# #       timepoint  = time_vec[r],
# #       village    = villages$village[village_vec[r]],
# #       latitude   = villages$latitude[village_vec[r]],
# #       longitude  = villages$longitude[village_vec[r]]
# #     )
# # }) %>%
# #   mutate(genotype   = do.call(paste0, across(all_of(locus_cols))),
# #     # genotype   = paste0(L1014F, L1014S),
# #          village_id = match(village, villages$village))
# # 
# # 
# # summary(mcmc_data_all_sim)
# # true_s
# # 
# # true_w[1:5, ]
# # true_s[1:5, ]
# # write.csv(mcmc_data_all_sim, "dataoutput/data_to_use.csv", row.names = FALSE)
# # # writeRaster(mcmc_data_all_sim, "dataoutput/data_to_use.tif")
# # 
# # # # 
# # # mcmc_data_all_sim <- mcmc_data_all_sim %>%
# # #   mutate(village_id = match(village, villages$village))
# # 
# # mcmc_data_all_sim <- mcmc_data_all_sim %>%
# #   mutate(genotype   = do.call(paste0, across(all_of(locus_cols))),
# #     # genotype = paste0(L1014F, L1014S),
# #          village_id = match(village, villages$village))
# # 
# # 
# # # allele-frequency table, long format: one row per pixel x timepoint x locus
# # allele_freq_long <- as.data.frame(true_allele_freq) %>%
# #   setNames(locus_names) %>%
# #   mutate(
# #     village   = villages$village[village_vec],
# #     timepoint = time_vec,
# #     village_id = village_vec
# #   ) %>%
# #   pivot_longer(cols = all_of(locus_names), names_to = "locus", values_to = "allele_frequency")
# # 
# # # we transform the data to matrix for the model to be able to read it
# # 
# # # fake allele counts — 2*M_z "trials" per row/locus, since alleles are diploid
# # eps <- 1e-6
# # # here we are multiplying the value from true_allele_freq to avoid it will fall to 0
# true_allele_freq_safe <- eps + (1 - 2 * eps) * true_allele_freq
# # generate a fake allele fequency, here because one loci contains 2 allele
# # that's why we say that the size will be 2* Mz
# fake_allele_count <- matrix(
#   rbinom(length(true_allele_freq_safe), size = 2 * M_z, prob = true_allele_freq_safe),
#   nrow = n_villages * Tmax
# )
# # matching the column names to match the loci
# colnames(fake_allele_count) <- locus_names
# colnames(fake_allele_count)
# 
# # modify the name to village_id and timepoint
# allele_count_long <- as.data.frame(fake_allele_count) %>%
#   mutate(village_id = village_vec, timepoint = time_vec) %>%
#   pivot_longer(cols = all_of(locus_names), names_to = "locus", values_to = "allele_count_observed")
# allele_freq_long <- allele_freq_long %>%
#   left_join(allele_count_long, by = c("village_id", "timepoint", "locus"))
# write.csv(allele_freq_long, "dataoutput/allele_frequency_data.csv", row.names = FALSE)
# 
# # fake slection pressure and relative fitness
# ##################################################################
# ### Extract true s and w through calculate()
# ##################################################################
# # sim_result <- calculate(
# #   alpha_matrix, theta, Z_matrix, betamat, h, rho_z, p_village,
# #   p_died_add_matrix, phi_add, allele_freq_matrix, phi_allele,
# #   s_matrix, w_matrix,
# #   nsim = 1
# # )
# 
# # true_s <- sim_result$s_matrix[1, , ]   # n_villages x n_loci
# # true_w <- sim_result$w_matrix[1, , ]   # n_villages x n_loci
# # true_betamat     <- sim_result$betamat[1, , ]
# # because we will do the same thing for the fake phenotype drawn from the additive
# # effect and multiplicative effect so we are using a function
# # pivot_genotype_matrix <- function(df, value_col) {
# #   df %>%
# #     mutate(across(any_of(c("village_id", "timepoint", "genotype_id")), as.numeric)) %>%
# #     arrange(village_id, timepoint, genotype_id) %>%
# #     pivot_wider(
# #       id_cols     = c(village_id, timepoint),
# #       names_from  = genotype_id,
# #       names_sort  = TRUE,
# #       values_from = {{ value_col }}
# #     ) %>%
# #     arrange(village_id, timepoint) %>%
# #     dplyr::select(-village_id, -timepoint) %>%
# #     as.matrix() %>%
# #     unname()
# # }
# 
# # fake_counts_matrix_pivoted <- pivot_genotype_matrix(mcmc_data_all_sim, n_observed)
# # #dead_add_matrix_pivoted    <- pivot_genotype_matrix(mcmc_data_all_sim, dead_add_observed)
# # # selection pressure
# # ##################################################################
# # selection_pressure_tbl <- as.data.frame(true_w) %>%
# #   setNames(locus_names) %>%
# #   mutate(village = villages$village, latitude = villages$latitude, longitude = villages$longitude) %>%
# #   tidyr::pivot_longer(cols = all_of(locus_names), names_to = "locus", values_to = "w") %>%
# #   left_join(
# #     as.data.frame(true_s) %>%
# #       setNames(locus_names) %>%
# #       mutate(village = villages$village) %>%
# #       tidyr::pivot_longer(cols = all_of(locus_names), names_to = "locus", values_to = "s"),
# #     by = c("village", "locus")
# #   )
# # 
# # selection_pressure_tbl
# # write.csv(selection_pressure_tbl, "dataoutput/selection_pressure.csv", row.names = FALSE)
# # 
# # # fake_counts_matrix_pivoted <- mcmc_data_all_sim %>%
# # #   arrange(village, timepoint, genotype_id) %>%
# # #   pivot_wider(
# # #     id_cols     = c(village, timepoint),
# # #     names_from  = genotype_id,
# # #     values_from = n_observed
# # #   ) %>%
# # #   arrange(village, timepoint) %>%
# # #   dplyr::select(-village, -timepoint) %>%
# # #   as.matrix()
# # fake_counts_matrix_pivoted  <- pivot_genotype_matrix(mcmc_data_all_sim, n_observed)
# # dead_add_matrix_pivoted     <- pivot_genotype_matrix(mcmc_data_all_sim, dead_add_observed)
# # # view(fake_counts_matrix_pivoted)
# # # dead_add_matrix_pivoted <- mcmc_data_all_sim %>%
# # #   mutate(across(c(village_id, timepoint, genotype_id), as.numeric)) %>%
# # #   arrange(village_id, timepoint, genotype_id) %>%
# # #   pivot_wider(
# # #     id_cols     = c(village_id, timepoint),
# # #     names_from  = genotype_id,
# # #     names_sort  = TRUE,
# # #     values_from = dead_add_observed
# # #   ) %>%
# # #   arrange(village_id, timepoint) %>%
# # #   dplyr::select(-village_id, -timepoint) %>%
# # #   as.matrix()
# # 
# # allele_count_matrix_pivoted <- allele_count_long %>%
# #   arrange(village_id, timepoint, locus) %>%
# #   pivot_wider(id_cols = c(village_id, timepoint), names_from = locus, values_from = allele_count_observed) %>%
# #   arrange(village_id, timepoint) %>%
# #   dplyr::select(-village_id, -timepoint) %>% 
# #   as.matrix()
# # 
# # stopifnot(all.equal(unname(allele_count_matrix_pivoted), unname(fake_allele_count)))
# # mcmc_data_all_sim %>% distinct(village, village_id) %>% arrange(village_id) %>% head(5)
# # ## slection pressure
# # # relative_fitness_long <- as.data.frame(true_w) %>%
# # #   setNames(locus_names) %>%
# # #   mutate(village = villages$village, latitude = villages$latitude, longitude = villages$longitude) %>%
# # #   pivot_longer(cols = all_of(locus_names), names_to = "locus", values_to = "w") %>%
# # #   left_join(
# # #     as.data.frame(true_s) %>%
# # #       setNames(locus_names) %>%
# # #       mutate(village = villages$village) %>%
# # #       pivot_longer(cols = all_of(locus_names), names_to = "locus", values_to = "s"),
# # #     by = c("village", "locus")
# # #   )
# # summary(relative_fitness_long$w)
# # head(relative_fitness_long )
# # write.csv(relative_fitness_long, "dataoutput/relative_fitness_long.csv", row.names = FALSE)
# # ## multiplicative
# # # dead_mult_matrix_pivoted <- mcmc_data_all_sim %>%
# # #   arrange(village_id, timepoint, genotype_id) %>%
# # #   pivot_wider(
# # #     id_cols     = c(village_id, timepoint),
# # #     names_from  = genotype_id,
# # #     values_from = dead_mult_observed
# # #   ) %>%
# # #   arrange(village_id, timepoint) %>%
# # #   dplyr::select(-village_id, -timepoint) %>%
# # #   as.matrix()
# # 
# # # we have to verify if it has the same format that data that our model is giving
# # #  here we named it fake_counts_matrix
# # 
# # fake_counts_matrix_pivoted <- apply(fake_counts_matrix_pivoted, 2, function(x) as.numeric(trimws(x)))
# # #dead_add_matrix_pivoted    <- apply(dead_add_matrix_pivoted,    2, function(x) as.numeric(trimws(x)))
# # 
# # # re-check
# # 
# # #stopifnot(all.equal(unname(dead_add_matrix_pivoted),  unname(fake_dead_add)))
# # # stopifnot(all.equal(unname(dead_mult_matrix_pivoted), unname(fake_dead_mult)))
# # stopifnot(inherits(allele_freq_matrix, "greta_array"))
# # 
# # stopifnot(all.equal(unname(fake_counts_matrix_pivoted), unname(fake_counts_matrix)))
# # 
# # observed_counts             <- as_data(fake_counts_matrix_pivoted)
# # #observed_counts_pheno_add   <- as_data(dead_add_matrix_pivoted)
# # observed_counts_allele     <- as_data(allele_count_matrix_pivoted)
# # # observed_counts_pheno_mult  <- as_data(dead_mult_matrix_pivoted)
# # # fitting the data by using his likelihood
# # ## genotype
# # distribution(observed_counts)            <- dirichlet_multinomial(size = size_vector, alpha = alpha_matrix)
# # 
# # 
# # # To avoid a new node and to fix the bug that was created here
# # # we tell R in advance that observed_counts_pheno_add is a data
# # #bserved_counts_pheno_add <- as_data(dead_add_matrix_pivoted)
# # 
# # # and then we could apply the likelihood by telling in advance how should 
# # # the alpha and beta in the betabinomial likelihood should be
# # # betabin(n, alpha, beta)
# # eps <- 1e-6
# # # this part is the mechanistic part correlating genotype and phenotype
# # # p_died_matrix comes from the theta and dominance
# # # multiplying it with the phi-add meaning the overdispersion of the values
# # # we consider that beta = 1 - mean, and alpha = mean * precision (phi-add)
# # # here the eps = 0.000001,
# # #p_died_add_safe <- eps + (1 - 2 * eps) * p_died_add_matrix
# # #numerically this equation above allow us to not falling to 0 for p and 
# # # will break the mcmc so if p_died_matrix = 0, we have 0.000001 + (1- 2*0.000001)*0
# # # = 0.0000001
# # alpha_beta_add <- p_died_add_safe * phi_add
# # beta_beta_add  <- (1 - p_died_add_safe) * phi_add
# # 
# # distribution(observed_counts_pheno_add) <- beta_binomial(size = M_z, alpha = alpha_beta_add, beta = beta_beta_add)
# # 
# # # allele frequency data likelihood
# # true_allele_freq_safe <- eps + (1 - 2 * eps) * true_allele_freq
# # fake_allele_count <- matrix(
# #   rbinom(length(true_allele_freq_safe), size = 2 * M_z, prob = true_allele_freq_safe),
# #   nrow = n_villages * Tmax
# # )
# # 
# # allele_freq_safe  <- eps + (1 - 2 * eps) * allele_freq_matrix
# # alpha_beta_allele <- allele_freq_safe * phi_allele
# # beta_beta_allele  <- (1 - allele_freq_safe) * phi_allele
# # 
# # distribution(observed_counts_allele) <- beta_binomial(size = 2 * M_z, alpha = alpha_beta_allele, beta = beta_beta_allele)
# # 
# # # alpha_beta_mult <- p_died_mult_matrix * phi_mult
# # # beta_beta_mult  <- (1 - p_died_mult_matrix) * phi_mult
# # # p_sample_mult   <- beta(alpha_beta_mult, beta_beta_mult)
# #  
# # ## looking on available distribution in greta package we have
# # #beta_binomial(size, alpha, beta, dim = NULL)
# # #distribution(observed_counts_pheno_add)  <- beta_binomial(size = M_z, alpha = alpha_beta_add, beta = beta_beta_add)
# # # distribution(observed_counts_pheno_add)  <- beta_binomial(size = M_z, prob = p_sample_add)
# # #distribution(observed_counts_pheno_mult)  <- beta_binomial(size = M_z, alpha = alpha_beta_mult, beta = beta_beta_mult)
# # # distribution(observed_counts_pheno_mult) <- binomial(size = M_z, prob = p_sample_mult)
# # # distribution(observed_counts_pheno_add)  <- binomial(size = M_z, prob = p_died_add_matrix)
# # # distribution(observed_counts_pheno_mult) <- binomial(size = M_z, prob = p_died_mult_matrix)
# # 
# # # estimation of the parameters byb using the model f unction of greta
# # geno_model <- model(betamat, h, rho_z, p_village, theta,  phi_add, phi_allele)
# # # geno_model <- model(betamat, h, rho_z, p_village, s_matrix)
# #  # 5- Plotting the DAG to see the nodes, plot the parameters by using mcmc_draws
# # # this code was trying to get the png of the dag but it didn't work
# # # png(filename = "almost_model.png", 
# # #     width = 280, height = 100, units = "mm", res = 200)
# # # dag <- plot(geno_model)
# # # print(dag)
# # # dev.off()
# #  
# # # here is am alternative
# # # library(DiagrammeR)
# # # library(DiagrammeRsvg)  
# # # library(rsvg)
# # dag <- plot(geno_model)
# # print(dag)
# # # this is the code to get the best version of the dag
# # svg_code <- export_svg(dag)
# # rsvg_png(charToRaw(svg_code), file = "dataoutput/2IRattempt_dag.png", width = 3000, height = 1200)
# # 
# # dev.list() # this is to check because here our code were stuck at the dag graph
# # dev.off() # this is to remove all images
# # 
# # # by using the mcmc, here we are trying to recover the parameters
# # # warmup named also burnin are the draw that would be not considered as they
# # # could consider as test for the draw
# # # chain is the chaine de valeur produite par le mcmc
# # # we think like we are investigated an area and identified a breeding site but
# # # didn't take coordinates (that's dumb), so we are sending others entomologists
# # # to investigate around the village to find it
# # # we could send 1 entomologist to do the task but we will be more confident
# # # if we send 4 entomologists that they will find the areas
# # draws <- greta::mcmc(
# #   model     = geno_model,
# #   n_samples = 200,
# #   warmup    = 200,
# #   chains    = 4
# # )
# # ####################################################
# # ### 7) Interpreting the results of the mcmc 
# # ####################################################
# # # investigating that the model given by mcmc is really giving what we are expecting
# # geno_model$dag$node_list |> length()
# # draws <- greta::mcmc(geno_model, n_samples = 50, warmup = 50, chains = 1)
# # colnames(as.matrix(draws))
# # 
# # library(bayesplot)
# # mcmc_trace(draws, regex_pars = c("rho_z"))
# # 
# # mcmc_trace(draws, regex_pars = "p_village\\[(1|2|3|4|5|6|7|8|9|10),")
# # 
# # mcmc_trace(draws, regex_pars = c("p_village"))
# # 
# # # compare the values of true priors with the posterior priors
# # true_p_village
# # true_h
# # true_betamat
# # 
# # # let's save our posterior in an object
# # posterior_draws <- as.matrix(draws)
# # 
# # 
# # # to get the summary statistic
# # statistic_summary <- summary(draws)$statistics
# # 
# # ####################################################################
# # ### 8) Comparing TRUE parameter values (used to simulate fake data)
# # ###    against the MCMC POSTERIOR draws — values + plots together
# # ####################################################################
# # ### 8.1 - Build a lookup: param name -> true value
# # # betamat[i,k]  (loci x covariates)
# # # reminder to myself that expand_grid is giving all possible combinaiton of
# # # a list, sp expand_grid (vect1, vect2)
# # # paste0 = create a text value
# # # so here we are saying for each possible combination, 
# # true_betamat <- matrix(as.numeric(true_betamat), nrow = n_loci, ncol = K)
# # betamat_truth <- expand_grid(row = seq_len(n_loci), col = seq_len(K)) %>%
# #   mutate(
# #     param = paste0("betamat[", row, ",", col, "]"),
# #     true_value = true_betamat[cbind(row, col)], # select the value of 
# #     # based on the row and col values, and put a matrix with
# #     family= "betamat",
# #     label = paste0("betamat,[", locus_names[row], ", cov", col, "]" )
# #   )
# # 
# # # because tibble is more for scalar and vectors
# # # and our h is more a vector, as it is changing based on the loci 
# # # dominance will also changed based on the type of insecticide
# # #  but as here we don't consider insecticide type yet then we
# # # we are assuming that dominance change based on gene
# # h_truth <- tibble(
# #   row= seq_len(n_loci),
# #   param = paste0("h[", row, ",1]"),
# #   true_value = as.numeric(true_h),
# #   family = "h",
# #   label = paste0("h[", locus_names, "]")
# # )
# # 
# # 
# # # rho_z is a scalar
# # rho_z_truth <- tibble(
# #   param = "rho_z",
# #   true_value = as.numeric(true_rho_z),
# #   family = "rho_z",
# #   label = "rho_z"
# # )
# # 
# # # p_village: which depend on the loci and the village
# # # next step should be in depend of the type of insecticide 
# # # so it is a matrix with row = n_villages and col = n_loci
# # p_truth <- expand.grid(row = seq_len(n_villages), col = seq_len(n_loci)) %>%
# #   mutate(
# #     param  = paste0("p_village[", row, ",", col, "]"),
# #     true_value = true_p_village[cbind(row, col)],
# #     family = "p_village",
# #     label  = paste0("p[", villages$village[row], ", ", locus_names[col])
# #   )
# # 
# # # here we are trying to get a dataframe to show the value of the 4 parameters
# # # parameters_combined <- bind_rows(betamat_truth, h_truth, rho_z_truth, p_truth)
# # parameters_combined<- bind_rows(
# #   betamat_truth %>% dplyr::select(param, true_value, family, label),
# #   h_truth        %>% dplyr::select(param, true_value, family, label),
# #   rho_z_truth      %>% dplyr::select(param, true_value, family, label),
# #   p_truth        %>% dplyr::select(param, true_value, family, label)
# # )
# # 
# # setdiff(parameters_combined$param, colnames(posterior_draws))
# # 
# # ##############################################  
# # ### 8.2 - Long-format posterior draws + attach truth
# # ##############################################
# # posterior_long <- as_tibble(posterior_draws) %>%
# #   mutate(draw = row_number()) %>%
# #   pivot_longer(-draw, names_to = "param", values_to = "posterior_value") %>%
# #   inner_join(parameters_combined, by = "param")
# 
# 
# rhat <- coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)
##############################################
### 8.3 - Numeric summary table: posterior vs truth
##############################################
# 
# recovery_summary <- posterior_long %>%
#   group_by(family, label, param, true_value) %>%
#   summarise(
#     post_mean   = mean(posterior_value),
#     post_median = median(posterior_value),
#     post_sd     = sd(posterior_value),
#     # the ci_low is telling us the 2.5% of the posterior value is low
#     # than the ci_low
#     ci_low      = quantile(posterior_value, 0.025),
#     # ci_high is telling us that that 97.5 % of the posterior value is 
#     # low than ci_high
#     ci_high     = quantile(posterior_value, 0.975),
#     rmse        = sqrt(mean((posterior_value - true_value)^2)),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     bias          = post_mean - true_value,
#     # so here we have the credible intervals that telling us that
#     # there is 95% that the true parameters values lies within that range
#     # so it means true_value is greater than the low_ci and less than ci_high
#     covered_95    = true_value >= ci_low & true_value <= ci_high,
#   ) %>%
#   arrange(family, label)
# 
# #write_csv(recovery_summary, "dataoutput/recovery_summary1.csv")  
# 
# # we read the recovery_summary here
# recovery_summary_weak <- read.csv("dataoutput/recovery_summary2.csv")
# recovery_summary_wide <- recovery_summary
# 
# # changing the name of the dataframe
# recovery_summary_tight <- recovery_summary_weak%>% mutate(prior_setup = "tight")
# recovery_summary_wide  <- recovery_summary_wide  %>% mutate(prior_setup = "wide")
# recovery_compare <- bind_rows(recovery_summary_tight, recovery_summary_wide)
# 
# # to get both data together
# recovery_diff <- recovery_summary_tight %>%
#   dplyr::select(family, label, param, true_value,
#          post_mean_tight = post_mean, ci_low_tight = ci_low, ci_high_tight = ci_high,
#          rmse_tight = rmse, bias_tight = bias, covered_95_tight = covered_95) %>%
#   inner_join(
#     recovery_summary_wide %>%
#       dplyr::select(param,
#              post_mean_wide = post_mean, ci_low_wide = ci_low, ci_high_wide = ci_high,
#              rmse_wide = rmse, bias_wide = bias, covered_95_wide = covered_95),
#     by = "param"
#   ) %>%
#   mutate(
#     ci_width_tight  = ci_high_tight - ci_low_tight,
#     ci_width_wide   = ci_high_wide  - ci_low_wide,
#     rmse_diff       = rmse_wide - rmse_tight,        
#     ci_width_diff   = ci_width_wide - ci_width_tight, 
#     rmse_pct_change = 100 * (rmse_wide - rmse_tight) / rmse_tight
#   ) %>%
#   arrange(family, label)
# 
# write.csv(recovery_diff, "dataoutput/comparison_param.csv")
# ## plotting rmse
# p_rmse_compare <- recovery_compare %>%
#   mutate(label = fct_reorder(label, rmse)) %>%
#   ggplot(aes(x = label, y = rmse, fill = prior_setup)) +
#   geom_col(position = position_dodge(width = 0.7), width = 0.6) +
#   coord_flip() +
#   facet_wrap(~family, scales = "free", ncol = 1) +
#   scale_fill_manual(values = c(tight = "steelblue", wide = "darkorange")) +
#   labs(title = "RMSE by parameter: tight vs wide prior", x = NULL, y = "RMSE", fill = "Prior setup") +
#   theme_minimal(base_size = 11)
# 
# print(p_rmse_compare)
# ##############################################
### 9- Plot : posterior density per parameter,
###       true value drawn as a vertical line
##############################################
# the function stat_halfeye is a function that allow us to plot the posterior density
# .width : sets which intervals to draw below the density
# about the fucntion stat_halfeye from ggdist https://r-statistics.co/ggdist-Package-in-R.html
plot_family_halfeye <- function(fam_name) {
  df <- posterior_long %>% filter(family == fam_name)
  truths <- parameters_combined %>% filter(family == fam_name)
  
  ggplot(df, aes(x = posterior_value, y = 0)) +
    stat_halfeye(
      fill = "#E8C55A",           # gold/yellow, matching the reference image
      color = "black",
      point_color = "black",
      interval_color = "black",
      .width = c(0.5, 0.95),       # thick bar = 50% CI, thin bar = 95% CI
      point_size = 2.2,
      slab_alpha = 1
    ) +
    geom_vline(
      data = truths,
      aes(xintercept = true_value),
      colour = "firebrick", linewidth = 0.9
    ) +
    facet_wrap(~label, scales = "free_x") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 12)
    )
}

p_betamat <- plot_family_halfeye("betamat")
p_h       <- plot_family_halfeye("h")
p_rho     <- plot_family_halfeye("rho_z")
p_p_village <- plot_family_halfeye("p_village")

keep <- sample(unique(parameters_combined$label[parameters_combined$family == "p_village"]), 2)
print(p_betamat)
posterior_long_bk      <- posterior_long
parameters_combined_bk <- parameters_combined

posterior_long      <- posterior_long      %>% filter(family != "p_village" | label %in% keep)
parameters_combined <- parameters_combined %>% filter(family != "p_village" | label %in% keep)

p_p_village <- plot_family_halfeye("p_village")
print(p_p_village)

## Posterior plot
betamat_truth <- expand_grid(row = seq_len(n_loci), col = seq_len(K)) %>%
  mutate(
    param      = paste0("betamat[", row, ",", col, "]"),
    true_value = true_betamat[cbind(row, col)],
    family     = "betamat",
    label      = paste0("betamat[Marker ", row, ", cov", col, "]")
  )

h_truth <- tibble(
  row        = seq_len(n_loci),
  param      = paste0("h[", row, ",1]"),
  true_value = as.numeric(true_h),
  family     = "h",
  label      = paste0("h[Marker ", row, "]")
)

rho_z_truth <- tibble(
  param      = "rho_z",
  true_value = as.numeric(true_rho_z),
  family     = "rho_z",
  label      = "rho_z"
)

# the closing bracket was missing in the original, which is why the panel
# titles read "p[px_kigali_24, L1014S" with nothing at the end
p_truth <- expand_grid(row = seq_len(n_villages), col = seq_len(n_loci)) %>%
  mutate(
    param      = paste0("p_village[", row, ",", col, "]"),
    true_value = true_p_village[cbind(row, col)],
    family     = "p_village",
    label      = paste0("p[", villages$village[row], ", Marker ", col, "]")
  )

#############################################################
### 9-a Selection pressure as a family
#############################################################
# s = exp(betamat %*% X) is deterministic, so I rebuild the draws from the
# betamat columns that are already in posterior_draws. No calculate(), no
# re-running the mcmc.
s_mat <- do.call(cbind, lapply(seq_len(n_loci), function(l) {
  b <- posterior_draws[, paste0("betamat[", l, ",", seq_len(K), "]"), drop = FALSE]
  m <- exp(b %*% t(X_villages))                      # n_draws x n_villages
  colnames(m) <- paste0("s[", seq_len(n_villages), ",", l, "]")
  m
}))

# the truth, rebuilt the same arithmetic way from true_betamat. This also
# sidesteps the bug where the second calculate(nsim = 1) overwrote true_s
# with a fresh draw from the prior.
true_s <- t(exp(true_betamat %*% t(X_villages)))     # n_villages x n_loci

s_truth <- expand_grid(row = seq_len(n_villages), col = seq_len(n_loci)) %>%
  mutate(
    param      = paste0("s[", row, ",", col, "]"),
    true_value = true_s[cbind(row, col)],
    family     = "s",
    label      = paste0("s[", villages$village[row], ", Marker ", col, "]")
  )

#############################################################
### 9-b Combine, then long format draws
#############################################################
parameters_combined <- bind_rows(
  betamat_truth %>% dplyr::select(param, true_value, family, label),
  h_truth       %>% dplyr::select(param, true_value, family, label),
  rho_z_truth   %>% dplyr::select(param, true_value, family, label),
  p_truth       %>% dplyr::select(param, true_value, family, label),
  s_truth       %>% dplyr::select(param, true_value, family, label)
)

# monitored parameters come straight from the draws
posterior_long <- as_tibble(posterior_draws) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "param", values_to = "posterior_value") %>%
  inner_join(parameters_combined, by = "param")

# s is not in posterior_draws, so it is appended separately
s_long <- as_tibble(s_mat[, s_truth$param, drop = FALSE]) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "param", values_to = "posterior_value") %>%
  inner_join(s_truth %>% dplyr::select(param, true_value, family, label), by = "param")

posterior_long <- bind_rows(posterior_long, s_long)

# every family should be listed here, including s
posterior_long %>% count(family) %>% print()

# nothing should be missing on the parameter side either
setdiff(parameters_combined$param, c(colnames(posterior_draws), colnames(s_mat)))
keep_labels <- c("p[px_kigali_15, Marker 1]",
                 "p[px_kigali_2, Marker 2]")

keep_labels_s <- sub("^p\\[", "s[", keep_labels)
keep_labels_s
# to plot for the selection
df_s     <- posterior_long      %>% filter(family == "s", label %in% keep_labels_s)
truths_s <- parameters_combined %>% filter(family == "s", label %in% keep_labels_s)


# to plot for the initial allele frequency
df_keep     <- posterior_long      %>% filter(family == "p_village", label %in% keep_labels)
truths_keep <- parameters_combined %>% filter(family == "p_village", label %in% keep_labels)

stopifnot(nrow(df_keep) > 0, nrow(truths_keep) == 2)

halfeye_panel <- function(df, truths, title = NULL) {
  ggplot(df, aes(x = posterior_value, y = 0)) +
    stat_halfeye(fill = "#E8C55A", color = "black",
                 point_color = "black", interval_color = "black",
                 .width = c(0.5, 0.95), point_size = 4, slab_alpha = 1) +
    geom_vline(data = truths, aes(xintercept = true_value),
               colour = "firebrick", linewidth = 1.4) +
    facet_wrap(~label, scales = "free_x", ncol = 2) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 20) +
    theme(strip.text   = element_text(face = "bold", size = 18),
          axis.text.x  = element_text(size = 18),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid   = element_blank(),
          plot.title   = element_text(face = "bold", size = 19))
}

#p_p_village <- halfeye_panel(df_p, truths_p, "Initial allele frequency")
p_s_village <- halfeye_panel(df_s, truths_s, "Selection pressure")

p_s_village
ggsave("dataoutput/p_village_28_37.png", p_p_village, width = 9, height = 4, dpi = 300)
ggsave("dataoutput/s_village_28_37.png", p_s_village, width = 9, height = 4, dpi = 300)

## ---- stacked, for the slide 
p_both <- p_p_village / p_s_village
ggsave("dataoutput/p_and_s_28_37.png", p_both, width = 9, height = 8, dpi = 300)

# Printing the plot
p_betamat   <- plot_family_halfeye("betamat")
p_h         <- plot_family_halfeye("h")
p_rho       <- plot_family_halfeye("rho_z")
p_p_village <- plot_family_halfeye("p_village", n_sample = 2, seed =10)

p_s         <- plot_family_halfeye("s",         n_sample = 4)

print(p_betamat)
print(p_h)
print(p_rho)
print(p_p_village)
print(p_s)

# save our works
ggsave("dataoutput/betamat_halfeye.png",     p_betamat,   width = 10, height = 6, dpi = 300)
ggsave("dataoutput/dominance.png",           p_h,         width = 8,  height = 5, dpi = 300)
ggsave("dataoutput/rho_z_halfeye.png",       p_rho,       width = 6,  height = 4, dpi = 300)
ggsave("dataoutput/p_village_4_sample.png",  p_p_village, width = 12, height = 8, dpi = 300)
ggsave("dataoutput/s_4_sample.png",          p_s,         width = 12, height = 8, dpi = 300)

#############################################################
### 10- Recovery numbers for the write up
#############################################################
recovery_summary <- posterior_long %>%
  group_by(family, label, param, true_value) %>%
  summarise(
    post_mean   = mean(posterior_value),
    post_median = median(posterior_value),
    post_sd     = sd(posterior_value),
    ci_low      = quantile(posterior_value, 0.025),
    ci_high     = quantile(posterior_value, 0.975),
    rmse        = sqrt(mean((posterior_value -  true_value)^2)),
    .groups = "drop"
  ) %>%
  mutate(
    bias       = post_mean - true_value,
    ci_width   = ci_high - ci_low,
    covered_95 = true_value >= ci_low & true_value <= ci_high
  ) %>%
  arrange(family, label)

# the per family line is what goes in the simulation study table. Four sampled
# panels can look fine while the other 192 do not, so read this before the plots
recovery_summary %>%
  group_by(family) %>%
  summarise(
    coverage_95   = mean(covered_95),
    mean_bias     = mean(bias),
    mean_rmse     = mean(rmse),
    mean_ci_width = mean(ci_width),
    .groups = "drop"
  ) %>%
  print()

write.csv(recovery_summary, "dataoutput/recovery_summary.csv", row.names = FALSE)

# n_sample = NULL keeps every panel, which is what I want for betamat, h and
# rho_z. p_village and s have 196 labels each, so those get sampled.
# Same seed across families means the same villages appear in both figures.
# plot_family_halfeye <- function(fam_name, n_sample = NULL, seed = 42) {
#   
#   df     <- posterior_long      %>% filter(family == fam_name)
#   truths <- parameters_combined %>% filter(family == fam_name)
#   
#   stopifnot(nrow(df) > 0)
#   
#   if (!is.null(n_sample)) {
#     set.seed(seed)
#     keep   <- sample(unique(truths$label), size = min(n_sample, n_distinct(truths$label)))
#     df     <- df     %>% filter(label %in% keep)
#     truths <- truths %>% filter(label %in% keep)
#   }
#   
#   ggplot(df, aes(x = posterior_value, y = 0)) +
#     stat_halfeye(
#       fill = "#E8C55A",            # gold/yellow, matching the reference image
#       color = "black",
#       point_color = "black",
#       interval_color = "black",
#       .width = c(0.5, 0.95),       # thick bar = 50% CI, thin bar = 95% CI
#       point_size = 2.2,
#       slab_alpha = 1
#     ) +
#     geom_vline(
#       data = truths,
#       aes(xintercept = true_value),
#       colour = "firebrick", linewidth = 0.9
#     ) +
#     facet_wrap(~label, scales = "free_x", ncol = 2) +
#     labs(x = NULL, y = NULL) +
#     theme_minimal(base_size = 13) +
#     theme(
#       axis.text.y  = element_blank(),
#       axis.ticks.y = element_blank(),
#       panel.grid   = element_blank(),
#       strip.text   = element_text(face = "bold", size = 12)
#     )
# }
# 
# 
# plot_family_halfeye <- function(fam_name, n_sample = NULL, seed = 42, exclude = NULL) {
#   
#   df     <- posterior_long      %>% filter(family == fam_name)
#   truths <- parameters_combined %>% filter(family == fam_name)
#   
#   # drop unwanted labels before anything else, so sampling cannot pick them
#   if (!is.null(exclude)) {
#     drop   <- grepl(paste(exclude, collapse = "|"), truths$label)
#     truths <- truths[!drop, ]
#     df     <- df %>% filter(label %in% truths$label)
#   }
#   
#   if (!is.null(n_sample)) {
#     set.seed(seed)
#     keep   <- sample(unique(truths$label), size = min(n_sample, n_distinct(truths$label)))
#     df     <- df     %>% filter(label %in% keep)
#     truths <- truths %>% filter(label %in% keep)
#   }
#   
#   ggplot(df, aes(x = posterior_value, y = 0)) +
#     stat_halfeye(
#       fill = "#E8C55A", color = "black",
#       point_color = "black", interval_color = "black",
#       .width = c(0.5, 0.95), point_size = 2.2, slab_alpha = 1
#     ) +
#     geom_vline(data = truths, aes(xintercept = true_value),
#                colour = "firebrick", linewidth = 0.9) +
#     facet_wrap(~label, scales = "free_x", ncol = 2) +
#     labs(x = NULL, y = NULL) +
#     theme_minimal(base_size = 13) +
#     theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#           panel.grid = element_blank(),
#           strip.text = element_text(face = "bold", size = 12))
# }
# 
# p_p_village <- plot_family_halfeye("p_village", n_sample = 4,
#                                    exclude = "px_kigali_13,")



class(x_cell_years)
class(fitness_cell_years)
class(effect_type)
