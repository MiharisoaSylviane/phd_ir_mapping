# ###############################################################
# Simulation Estimation Study : Multi-locus Model (using Greta)
# ###############################################################
# This script is executing those following in the order:
#   1) call the fix data like dummy matrices data L/R, villages, and number of loci
#   2) define priors (betamat (effect of intervention), dominance (h), overdispersion (rho_z),
#   initial allele_frequency (p))
#   3) calling the function of the probability of having one genotype and getting one genotype
#   4) Fake data generating with the prior and likelihood
#   5) Plotting the DAG to see the nodes, plot the parameters by using mcmc_draws
#   6) Estimate the parameters (betamat, h, rho-z and p) by using mcmc package and the fake data
#   7) Interpreting the results of the mcmc 

### Loading Library
library(tidyverse)
library(MCMCpack)
library(coda)
library(greta)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

set.seed(123)
#############################################################
### 1- Defining all fix data, we consider as known
#############################################################
# Dummy matrices 
source("R/create_dummy_matrices.R")

# because we want to get a good DAG then we will reduce to 2 loci first
n_loci <- 2

# creation of the matrix genotype in row and loci in the column
mats   <- create_dummy_matrices(n_loci)
L      <- mats$left
R      <- mats$right
G      <- nrow(L)
B      <- ncol(L) # B = number of loci 

# check the dimension of the object L
dim(L)

# naming the locus for the sake of the dataframe later
locus_names <- c("L1014F", "L1014S")
colnames(L) <- locus_names
colnames(R) <- locus_names

# locus_names <- c("L1014F", "L1014S", "Ace1")

# defining our village 
villages <- tibble(
  village   = c("Nkolondom", "Nkolbikon"),
  latitude  = c( 3.948,  5.602),
  longitude = c(11.505, 13.675)
)
n_villages <- nrow(villages)

# number of generations
Tmax <- 2

# number of mosquitoes tested
M_z <- 100


## Covariates per village
# K : number of covariates, X_village = the intensity of use in the village
K <- 2
X_villages <- cbind(
  intercept = 1,
  coverage  = c(0.10, 0.35)
)
# X_villages <- cbind(
#   intercept = 1,
#   coverage  = c(0.10, 0.35, 0.50, 0.65, 0.80, 0.95)
# )

# Locus transformation to match the Greta syntax
# version greta as the matrix creation won't be identified by Greta
# because Greta can't keep the matrix with qualitative data part of the plain R we have then 
# here we are transforming the L==1L to boleen TRUE OR FALSE
# 
genotype_label_df <- map_dfc(locus_names, function(loc) {
  tibble(!!loc := case_when(
    L[, loc] == 1 & R[, loc] == 1 ~ "SS",
    L[, loc] == 0 & R[, loc] == 0 ~ "RR",
    TRUE                          ~ "SR"
  ))
})

genotype_lookup <- genotype_label_df %>%
  mutate(genotype_id = row_number(), .before = 1)

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

##################################################################
### 2- Defining Priors
##################################################################
# betamat the effect covariate on the selection, h: dominance, 
# rho-z : overdispersion that will be used in the dirichlet distribution
# p_village : initial allele frequency
betamat   <- normal(0, 1, dim = c(n_loci, K))             
h         <- beta(2, 2, dim = n_loci)                      
rho_z     <- uniform(0, 1)                                 
p_village <- uniform(0, 1, dim = c(n_villages, n_loci))    

# 3- calling the function of the probability of having one genotype and getting one genotype
source("R/function_geno_pheno.R")

# Build a list to put the list of rho_z by village 
# this equation is needed for the Dirichlet distribution
alpha0_z <- (1 - rho_z^2) / rho_z^2

# creating a vector where we could store a list of dimension village X time
# it was equivalent to the step where we are creating vector or matrix for data storage
alpha_rows  <- vector("list", n_villages * Tmax)
# so this will be null or empty then we will fill that with the 
Z_rows      <- vector("list", n_villages * Tmax)
village_vec <- integer(n_villages * Tmax)
time_vec    <- integer(n_villages * Tmax)
row_id <- 1

# function to simulate

for (i in seq_len(n_villages)) {
  
  w_i <- compute_w_greta(betamat, X_villages[i, ])
  # transpose is used here because we have (p1,p2,p3) which is 1 X 3
  # however those should be in the column like
  # [p1]
  # [p2]
  # [p3] which give us (3 X 1)
  p_i <- t(p_village[i, ])   
  # 
  Z_list_i <- simulate_genotype_timecourse_greta(
    p = p_i, w = w_i, h = h, Tmax = Tmax,
    SS_mask = SS_mask_g, RR_mask = RR_mask_g, SR_mask = SR_mask_g,
    L = L_g, R = R_g
  )
  
  for (t in seq_len(Tmax)) {
    # genotype frequency (1 X G)
    Z_rows[[row_id]]     <- t(Z_list_i[[t]])                    
    alpha_rows[[row_id]] <- Z_rows[[row_id]] * alpha0_z + 1e-8   
    village_vec[row_id]  <- i
    time_vec[row_id]      <- t
    row_id <- row_id + 1
  }
}

# (n_villages*Tmax) x G, true frequency
Z_matrix     <- do.call(greta::abind, c(Z_rows, list(along = 1)))  
# (n_villages*Tmax) x G
alpha_matrix <- do.call(greta::abind, c(alpha_rows, list(along = 1)))  
size_vector  <- rep(M_z, length(alpha_rows))  

#########################################################
### 4- Fake data generating with the prior and likelihood
#########################################################
sim_result <- calculate(alpha_matrix, Z_matrix, betamat, h, rho_z, p_village, nsim = 1)
# (n_villages*Tmax) x G
alpha_numeric  <- sim_result$alpha_matrix[1, , ] 
# (n_villages*Tmax) x G
true_Z_matrix  <- sim_result$Z_matrix[1, , ]       
true_betamat   <- sim_result$betamat[1, , ]
true_h         <- sim_result$h[1, , ]
true_rho_z     <- as.numeric(sim_result$rho_z)[1]
true_p_village <- sim_result$p_village[1, , ]

# Likelihood
fake_counts_matrix <- t(apply(alpha_numeric, 1, function(a) {
  z_disp <- as.numeric(MCMCpack::rdirichlet(1, a))
  as.vector(rmultinom(1, M_z, z_disp / sum(z_disp)))
}))

n_tested_vec <- rowSums(fake_counts_matrix) 

# tibble of the fake data we tried before: row =  village x génotype x timepoint
mcmc_data_all_sim <- map_dfr(seq_len(nrow(fake_counts_matrix)), function(r) {
  genotype_lookup %>%
    mutate(
      n_observed = fake_counts_matrix[r, ],
      z_true     = true_Z_matrix[r, ],
      n_tested   = n_tested_vec[r],
      timepoint  = time_vec[r],
      village    = villages$village[village_vec[r]],
      latitude   = villages$latitude[village_vec[r]],
      longitude  = villages$longitude[village_vec[r]]
    )
})

# 
mcmc_data_all_sim <- mcmc_data_all_sim %>%
  mutate(village_id = match(village, villages$village))
# we transform the data to matrix for the model to be able to read it
fake_counts_matrix_pivoted <- mcmc_data_all_sim %>%
  arrange(village_id, timepoint, genotype_id) %>%
  pivot_wider(
    id_cols     = c(village_id, timepoint),
    names_from  = genotype_id,
    values_from = n_observed
  ) %>%
  arrange(village_id, timepoint) %>%
  dplyr::select(-village_id, -timepoint) %>%
  as.matrix()
# view(fake_counts_matrix_pivoted)

# we have to verify if it has the same format that data that our model is giving
#  here we named it fake_counts_matrix
stopifnot(all.equal(unname(fake_counts_matrix_pivoted), unname(fake_counts_matrix)))

# so we are considering it as the number of mosquitoes tested positive
observed_counts <- as_data(fake_counts_matrix_pivoted)

# then we apply the likelihood again
distribution(observed_counts) <- dirichlet_multinomial(size = size_vector, alpha = alpha_matrix)

geno_model <- model(betamat, h, rho_z, p_village)


# 5- Plotting the DAG to see the nodes, plot the parameters by using mcmc_draws
# this code was trying to get the png of the dag but it didn't work
# png(filename = "almost_model.png", 
#     width = 280, height = 100, units = "mm", res = 200)
# dag <- plot(geno_model)
# print(dag)
# dev.off()

# here is am alternative
# library(DiagrammeR)
# library(DiagrammeRsvg)
# library(rsvg)
dag <- plot(geno_model)
# this is the code to get the best version of the dag
svg_code <- export_svg(dag)
rsvg_png(charToRaw(svg_code), file = "dataoutput/1IRattempt_dag.png", width = 3000, height = 1200)

dev.list() # this is to check because here our code were stuck at the dag graph
dev.off() # this is to remove all images

###################################################################
# 6- Estimate the parameters (betamat, h, rho-z and p) by using mcmc package and the fake data
#################################################################

# by using the mcmc, here we are trying to recover the parameters
# warmup named also burnin are the draw that would be not considered as they
# could consider as test for the draw
# chain is the chaine de valeur produite par le mcmc
# we think like we are investigated an area and identified a breeding site but
# didn't take coordinates (that's dumb), so we are sending others entomologists
# to investigate around the village to find it
# we could send 1 entomologist to do the task but we will be more confident
# if we send 4 entomologists that they will find the areas
draws <- greta::mcmc(
  model     = geno_model,
  n_samples = 1000,
  warmup    = 1000,
  chains    = 4
)
####################################################
### 7) Interpreting the results of the mcmc 
####################################################
# investigating that the model given by mcmc is really giving what we are expecting

library(bayesplot)
bayesplot::mcmc_trace(draws, regex_pars = c("rho_z", "h"))

mcmc_trace(draws, regex_pars = c("p_village"))

estim_summary <- summary(draws)$statistics
<<<<<<< HEAD
=======

>>>>>>> 99f374446a218cd065491aced1d2cc001cac0b68
