# ###############################################################
# 1- BUILDING THE FAKE DATA
# a) preparing a table to store the data
# b) and store the data
# 2- FITTING TO ESTIMATE THE PARAMETERS
# This script is executing those following in the order:
#   1) call the fix data like dummy matrices data L/R, villages, and number of loci
#   2) Preparing the covariates and read the data
#   3) Preparing the data
#   3) calling the function of the phenotype of each genotype
#   4) Fake data generating with the prior and likelihood
#   5) Plotting the DAG to see the nodes, plot the parameters by using mcmc_draws
#   6) Estimate the parameters (betamat, h, rho-z and p) by using mcmc package and the fake data
#   7) Interpreting the results of the mcmc 
source("R/Packages.R") 
source("R/create_dummy_matrices.R")
source("R/function_geno_pheno.R")
# Dummy matrices 


# number of generations
Tmax <- 2

# number of mosquitoes tested
M_z <- 10
years <- 2018:2021
length(years)
# this is the jump from one year to the next
n_trans <- length(years) - 1L
# the year we have data about the allele frequency
Tmax    <- n_trans + 1L            
stopifnot(is.numeric(years)) 


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
names(nets_cube)
# IRS coverage
irs_cube <- rast("data_raw/raster_clean/irs_coverage_scaled_cube.tif")

# human population
pop_cube <- rast("data_raw/raster_clean/pop_scaled_cube.tif")
names(pop_cube)
# for each one pad back to the baseline year, repeating the first
nets_cube <- pre_pad_cube(nets_cube, baseline_year)
irs_cube <- pre_pad_cube(irs_cube, baseline_year)
pop_cube <- pre_pad_cube(pop_cube, baseline_year)

# if necessary post-pad
# if necessary, pad forward to the final data year, repeating the last
# nets_cube <- post_pad_cube(nets_cube, final_data_year)
irs_cube <- post_pad_cube(irs_cube, final_data_year)
pop_cube <- post_pad_cube(pop_cube, final_data_year)


# load the non-temporal crop cov ariate layers

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

### The cells: all the pixels of Kigali
#############################################################
cells <- as.data.frame(kigali_raster, cells = TRUE, xy = TRUE, na.rm = TRUE) %>%
  select(cell, longitude = x, latitude = y)
unique_cells   <- cells$cell
n_cells        <- length(unique_cells)
n_unique_cells <- n_cells


#############################################################
### Covariates for each cell and year -> x_cell_years
#############################################################
# the covariates of year t make the p of year t+1 -> 2018, 2019, 2020
# extract spatiotemporal covariates from the cube
all_extract <- bind_cols(
  terra::extract(nets_k, unique_cells),
  terra::extract(irs_k,  unique_cells),
  terra::extract(pop_k,  unique_cells)
) %>%
  mutate(cell = unique_cells, .before = everything()) %>%
  pivot_longer(cols = -one_of("cell"), names_sep = "_",
               names_to = c("variable", "year"), values_to = "value") %>%
  pivot_wider(names_from = "variable", values_from = "value") %>%
  mutate(year = as.numeric(year)) %>%
  filter(year %in% years[-Tmax]) %>%
  mutate(cell_id = match(cell, unique_cells)) %>%
  arrange(cell_id, year)

head(all_extract) 


x_cell_years <- all_extract %>% select(nets, irs, pop) %>% as.matrix()
K <- ncol(x_cell_years)

stopifnot(nrow(x_cell_years) == n_cells * n_trans,
          !anyNA(x_cell_years))


#############################################################
### 1- Possible genotype combination
#############################################################
# number of loci
n_loci <- 2L

# creation of the matrix genotype in row and loci in the column
mats   <- create_dummy_matrices(n_loci)
L      <- mats$left
R      <- mats$right

# check the dimension of the object L
dim(L)

# naming the locus for the sake of the dataframe later
locus_names <- c("marker_1", "marker_2")
colnames(L) <- locus_names
colnames(R) <- locus_names

# Locus transformation to match the Greta syntax
# version greta as the matrix creation won't be identified by Greta
# because Greta can't keep the matrix with qualitative data part of the plain R we have then 
# here we are transforming the L==1L to boleen TRUE OR FALSE
### defining the combination SS
SS_mask <- (L == 1L) & (R == 1L)
RR_mask <- (L == 0L) & (R == 0L) 
SR_mask <- (L != R)
# here we are saving them as 0 or 1
storage.mode(SS_mask) <- "double"
storage.mode(RR_mask) <- "double"
storage.mode(SR_mask) <- "double"

# defining those as variables to be recognized by Greta
# as_data() is saying, take the value of those to be the data used in our model
L_g       <- as_data(L) 
R_g       <- as_data(R)
SS_mask_g <- as_data(SS_mask) 
RR_mask_g <- as_data(RR_mask)
SR_mask_g <- as_data(SR_mask)

# c) a 9 x 9 matrix of ones, to sum each row
ones_G <- ones(nrow(SS_mask), nrow(SS_mask))
##################################################################
### 2- Defining Priors
##################################################################
# betamat the effect covariate on the selection, h: dominance, 
# rho-z : overdispersion that will be used in the dirichlet distribution
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
# Build a list to put the list of rho_z by village 
# this equation is needed for the Dirichlet distribution
alpha0_z <- (1 - rho_z^2) / rho_z^2
### Betamat represent the effect of covariates on the selection
betamat   <- normal(-2, 0.7, dim = c(K, n_loci))

# # ### overdispersion for the phenotype/bioassay data
# # phi_add  <- lognormal(log(20), 0.5) 
# # 
# # ### overdisperion for the allele frequency (hyperparameter)
phi_allele <- lognormal(log(20), 0.5)
# 
# ### theta is the locus-specific, insecticide-specific effect\
# # meaning the contribution of a particular genetic locus to resistance against 
# #  a particular insecticide
# theta <- exponential(0.5, dim = c(n_unique_cells * Tmax, n_loci))
# # theta <- uniform(0, 1, dim = c(n_villages * Tmax, n_loci))

#####Selection and fitness
# convert beta to positive effect sizes, them matrix multiply with indicator
# variables to constrain to positive effects of the covariates on selection
# coefficient/relative fitness
effect_type <- exp(betamat)

# # compute selection coefficients for cell-years, convert to relative fitness,
# and extract at data locations
# selection pressure
selection_cell_years <- x_cell_years %*% effect_type
dim(selection_cell_years)
# check
# dim(selection_cell_years)
# dim(x_cell_years)
# dim(effect_type)
# relative fitness
fitness_cell_years <- 1 + selection_cell_years
dim(fitness_cell_years)
calculate(head(fitness_cell_years), nsim =1)
# reformat this in to a 3D array with dimensions:
#   n_times x n_unique_cells x n_types x 1
# to solve dynamics with time-varying fitness (time must be first, then other
# two must match state variable, which has a trailing dimension of size 1)
fitness_array <- fitness_cell_years
dim(fitness_array) <- c(n_trans, n_cells, n_loci, 1)
dim(fitness_array)

# fitness of each genotype
h_mat <- ones(n_cells * n_trans, 1) %*% t(h)  # Gw and r_vec, built once, fitness does not vary in time 
# dominance should be a matrix with n_cells * n_loci 
# reminder ones(m,)
# dimension  =  n_cells x n_loci
# 3- calling the function of the probability of having one genotype and getting one genotype

# dimensions for the polygenic score n_cells x G
# this is to estimate the probablity of having one genotype X
resistant <- log(fitness_cell_years) %*% as_data(t(RR_mask))
heterozygous <- log(h_mat * fitness_cell_years + (1 - h_mat)) %*% as_data(t(SR_mask))
Gw <-  resistant + heterozygous
dim(Gw)

# dimensions should be  G x n_cells
r_vec <- exp(Gw)
dim(r_vec) <- c(n_trans, n_cells, nrow(SS_mask), 1)      # 30 x 49 x 9 x 1


# dimensions should be n_cells x G
Z_init     <- probability_genotype_rows_greta(p_cells)  
# defining the function for the change
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
# dim(initial_state)
# class(initial_state)
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


# allele frequency
n_dims_states <- length(dim(dynamics$all_states))

# p of jump k, as n_cells x n_loci
get_year <- function(k) {
  if (n_dims_states == 4) {
    p <- dynamics$all_states[, , 1, k]   # n_loci x n_cells x 1 x n_trans
  } else {
    p <- dynamics$all_states[, , k]      # n_loci x n_cells x n_trans
  }
  dim(p) <- c(n_loci, n_cells)
  t(p)                                   # n_cells x n_loci
}

# 2018 (p_cells), then 2019, 2020, 2021 stacked under each other
allele_freq_matrix <- do.call(rbind, c(list(p_cells), lapply(seq_len(n_trans), get_year)))
# check
# chk <- calculate(dynamics$all_states, nsim = 1)
# st  <- chk$all_states[1, , , , ]
# sim_result <- calculate(betamat, h, rho_z, p_cells, selection_cell_years,Z_init, fitness_cell_years, nsim = 1)
#############################################################
### 4- Read the fake data
#############################################################
fake_data <- read.csv("dataoutput/fake_data_kigali_2loci.csv")
names(fake_data)
# keep only what a field team would see: drop the truth columns if they are there
fake_data <- fake_data %>%
  select(-any_of(c("p_true", "z_SS", "z_RS", "z_RR")))

# checks: the data match the model (same cells, loci and years)
stopifnot(all(fake_data$cell  %in% unique_cells),
          all(fake_data$locus %in% locus_names),
          all(fake_data$year  %in% years))

#############################################################
### 5- Find the p of each data row in allele_freq_matrix
#############################################################
# the list of rows of allele_freq_matrix, in the SAME order (year first, cell inside)
matrix_rows <- expand_grid(year = years, cell = unique_cells) %>%
  mutate(row = row_number())

head(matrix_rows)

# each data row gets its row number by joining on year and cell
fake_data <- fake_data %>%
  left_join(matrix_rows, by = c("year", "cell"))

stopifnot(!anyNA(fake_data$row))   # every data row found its row


# split the data: one table per locus and data type
geno_m1 <- fake_data %>% filter(data_type == "genotype",  locus == locus_names[1])
geno_m2 <- fake_data %>% filter(data_type == "genotype",  locus == locus_names[2])
freq_m1 <- fake_data %>% filter(data_type == "frequency", locus == locus_names[1])

c(geno_m1 = nrow(geno_m1), geno_m2 = nrow(geno_m2), freq_m1 = nrow(freq_m1))

# likelihood 
# a) the p of the model at each data row: the row from the join, the column = the locus
p_geno_m1 <- allele_freq_matrix[geno_m1$row, 1]
p_geno_m2 <- allele_freq_matrix[geno_m2$row, 2]
p_freq_m1 <- allele_freq_matrix[freq_m1$row, 1]

# b) genotype counts, marker_1
z_m1   <- cbind((1 - p_geno_m1)^2, 2 * p_geno_m1 * (1 - p_geno_m1), p_geno_m1^2)
obs_m1 <- as_data(as.matrix(geno_m1[, c("n_SS", "n_RS", "n_RR")]))
distribution(obs_m1) <- dirichlet_multinomial(size  = geno_m1$n_tested,
                                              alpha = z_m1 * alpha0_z + 1e-8)

# c) genotype counts, marker_2
z_m2   <- cbind((1 - p_geno_m2)^2, 2 * p_geno_m2 * (1 - p_geno_m2), p_geno_m2^2)
obs_m2 <- as_data(as.matrix(geno_m2[, c("n_SS", "n_RS", "n_RR")]))
distribution(obs_m2) <- dirichlet_multinomial(size  = geno_m2$n_tested,
                                              alpha = z_m2 * alpha0_z + 1e-8)
# allele counts
# d) frequency rows, marker_1: R alleles among 2 x n_tested
obs_f1 <- as_data(freq_m1$n_R_allele)
distribution(obs_f1) <- beta_binomial(size  = 2 * freq_m1$n_tested,
                                      alpha = p_freq_m1 * phi_allele,
                                      beta  = (1 - p_freq_m1) * phi_allele)

# build the model
m <- model(
  # dominance
  h,
  # effect of covariates
  betamat,
  p_cells,
  rho_z, phi_allele
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

post   <- do.call(rbind, lapply(draws, as.matrix))
p1_hat <- matrix(colMeans(post[, grep("^p_cells", colnames(post))]),
                 n_cells, n_loci)

# truth vs posterior
p_truth <- read.csv("dataoutput/fake_data_kigali_2loci.csv") %>%   # <- the data file
  mutate(locus = tolower(trimws(locus))) %>%
  select(year, cell, locus, p_true)

# p_cells = the p of the FIRST year (2018)
# greta names it p_cells[i,j]: i = position of the cell in unique_cells, j = locus
truth_pcells <- p_truth %>%
  filter(year == min(years)) %>%
  mutate(i         = match(cell, unique_cells),
         j         = match(locus, locus_names),
         parameter = paste0("p_cells[", i, ",", j, "]")) %>%
  select(parameter, cell, locus, true = p_true)

setdiff(truth_pcells$parameter, colnames(post))   # must be EMPTY

set.seed(1)
some <- sample(truth_pcells$parameter, 8)     # 8 random cell-loci

post_pc <- as_tibble(post[, some]) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value")

ggplot(post_pc, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +   # the prior: uniform(0, 1)
  geom_vline(data = filter(truth_pcells, parameter %in% some),
             aes(xintercept = true), colour = "red") +
  facet_wrap(~ parameter, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "p_cells: posterior (blue), prior (dashed), truth (red)", x = "p in 2018", y = NULL)

# 3) the model's estimate for 2018 = posterior mean of p_cells[i, j]
p_hat_2018 <- expand_grid(i = seq_len(n_cells), j = seq_len(n_loci)) %>%
  mutate(cell      = unique_cells[i],
         locus     = locus_names[j],
         year      = min(years),
         parameter = paste0("p_cells[", i, ",", j, "]"),
         p_hat     = colMeans(post[, parameter])) %>%
  select(year, cell, locus, p_hat)


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
