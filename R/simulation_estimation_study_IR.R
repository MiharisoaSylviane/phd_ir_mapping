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
### Priors 1
#betamat   <- normal(0, 1, dim = c(n_loci, K))
betamat  <- normal(0, 5, dim = c(n_loci, K))
#### dominance
#h         <- beta(2, 2, dim = n_loci)  # this is saying that the doominance
# h is spinning around 0.5
h      <- uniform(0, 1, dim = n_loci)
#### overdispersion
rho_z     <- uniform(0, 1)                                 
p_village <- uniform(0, 1, dim = c(n_villages, n_loci)) 

### 

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

# fitting the data by using his likelihood
distribution(observed_counts) <- dirichlet_multinomial(size = size_vector, alpha = alpha_matrix)

# estimation of the parameters byb using the model function of greta
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
mcmc_trace(draws, regex_pars = c("rho_z", "h"))

mcmc_trace(draws, regex_pars = c("p_village"))

# compare the values of true priors with the posterior priors
true_p_village
true_h
true_betamat

# let's save our posterior in an object
posterior_draws <- as.matrix(draws)


# to get the summary statistic
statistic_summary <- summary(draws)$statistics

####################################################################
### 8) Comparing TRUE parameter values (used to simulate fake data)
###    against the MCMC POSTERIOR draws — values + plots together
####################################################################
### 8.1 - Build a lookup: param name -> true value
# betamat[i,k]  (loci x covariates)
# reminder to myself that expand_grid is giving all possible combinaiton of
# a list, sp expand_grid (vect1, vect2)
# paste0 = create a text value
# so here we are saying for each possible combination, 
betamat_truth <- expand_grid(row = seq_len(n_loci), col = seq_len(K)) %>%
  mutate(
    param = paste0("betamat[", row, ",", col, "]"),
    true_value = true_betamat[cbind(row, col)], # select the value of 
    # based on the row and col values, and put a matrix with
    family= "betamat",
    label = paste0("betamat,[", locus_names[row], ", cov", col, "]" )
  )


# because tibble is more for scalar and vectors
# and our h is more a vector, as it is changing based on the loci 
# dominance will also changed based on the type of insecticide
#  but as here we don't consider insecticide type yet then we
# we are assuming that dominance change based on gene
h_truth <- tibble(
  row= seq_len(n_loci),
  param = paste0("h[", row, ",1]"),
  true_value = as.numeric(true_h),
  family = "h",
  label = paste0("h[", locus_names, "]")
)


# rho_z is a scalar
rho_z_truth <- tibble(
  param = "rho_z",
  true_value = as.numeric(true_rho_z),
  family = "rho_z",
  label = "rho_z"
)

# p_village: which depend on the loci and the village
# next step should be in depend of the type of insecticide 
# so it is a matrix with row = n_villages and col = n_loci
p_truth <- expand.grid(row = seq_len(n_villages), col = seq_len(n_loci)) %>%
  mutate(
    param  = paste0("p_village[", row, ",", col, "]"),
    true_value = true_p_village[cbind(row, col)],
    family = "p_village",
    label  = paste0("p[", villages$village[row], ", ", locus_names[col])
  )

# here we are trying to get a dataframe to show the value of the 4 parameters
# parameters_combined <- bind_rows(betamat_truth, h_truth, rho_z_truth, p_truth)
parameters_combined<- bind_rows(
  betamat_truth %>% dplyr::select(param, true_value, family, label),
  h_truth        %>% dplyr::select(param, true_value, family, label),
  rho_z_truth      %>% dplyr::select(param, true_value, family, label),
  p_truth        %>% dplyr::select(param, true_value, family, label)
)

setdiff(parameters_combined$param, colnames(posterior_draws))

##############################################
### 8.2 - Long-format posterior draws + attach truth
##############################################
posterior_long <- as_tibble(posterior_draws) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "param", values_to = "posterior_value") %>%
  inner_join(parameters_combined, by = "param")

##############################################
### 8.3 - Numeric summary table: posterior vs truth
##############################################

recovery_summary <- posterior_long %>%
  group_by(family, label, param, true_value) %>%
  summarise(
    post_mean   = mean(posterior_value),
    post_median = median(posterior_value),
    post_sd     = sd(posterior_value),
    # the ci_low is telling us the 2.5% of the posterior value is low
    # than the ci_low
    ci_low      = quantile(posterior_value, 0.025),
    # ci_high is telling us that that 97.5 % of the posterior value is 
    # low than ci_high
    ci_high     = quantile(posterior_value, 0.975),
    rmse        = sqrt(mean((posterior_value - true_value)^2)),
    .groups = "drop"
  ) %>%
  mutate(
    bias          = post_mean - true_value,
    # so here we have the credible intervals that telling us that
    # there is 95% that the true parameters values lies within that range
    # so it means true_value is greater than the low_ci and less than ci_high
    covered_95    = true_value >= ci_low & true_value <= ci_high,
  ) %>%
  arrange(family, label)

#write_csv(recovery_summary, "dataoutput/recovery_summary1.csv")  

# we read the recovery_summary here
recovery_summary_weak <- read.csv("dataoutput/recovery_summary2.csv")
recovery_summary_wide <- recovery_summary

# changing the name of the dataframe
recovery_summary_tight <- recovery_summary_weak%>% mutate(prior_setup = "tight")
recovery_summary_wide  <- recovery_summary_wide  %>% mutate(prior_setup = "wide")
recovery_compare <- bind_rows(recovery_summary_tight, recovery_summary_wide)

# to get both data together
recovery_diff <- recovery_summary_tight %>%
  dplyr::select(family, label, param, true_value,
         post_mean_tight = post_mean, ci_low_tight = ci_low, ci_high_tight = ci_high,
         rmse_tight = rmse, bias_tight = bias, covered_95_tight = covered_95) %>%
  inner_join(
    recovery_summary_wide %>%
      dplyr::select(param,
             post_mean_wide = post_mean, ci_low_wide = ci_low, ci_high_wide = ci_high,
             rmse_wide = rmse, bias_wide = bias, covered_95_wide = covered_95),
    by = "param"
  ) %>%
  mutate(
    ci_width_tight  = ci_high_tight - ci_low_tight,
    ci_width_wide   = ci_high_wide  - ci_low_wide,
    rmse_diff       = rmse_wide - rmse_tight,        
    ci_width_diff   = ci_width_wide - ci_width_tight, 
    rmse_pct_change = 100 * (rmse_wide - rmse_tight) / rmse_tight
  ) %>%
  arrange(family, label)

write.csv(recovery_diff, "dataoutput/comparison_param.csv")
## plotting rmse
p_rmse_compare <- recovery_compare %>%
  mutate(label = fct_reorder(label, rmse)) %>%
  ggplot(aes(x = label, y = rmse, fill = prior_setup)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  coord_flip() +
  facet_wrap(~family, scales = "free", ncol = 1) +
  scale_fill_manual(values = c(tight = "steelblue", wide = "darkorange")) +
  labs(title = "RMSE by parameter: tight vs wide prior", x = NULL, y = "RMSE", fill = "Prior setup") +
  theme_minimal(base_size = 11)

print(p_rmse_compare)
##############################################
### 8.4 - Plot 1: posterior density per parameter,
###       true value drawn as a vertical line
##############################################

plot_family <- function(fam_name) {
  df <- posterior_long %>% filter(family == fam_name)
  truths <- parameters_combined %>% filter(family == fam_name)
  
  ggplot(df, aes(x = posterior_value)) +
    geom_density(fill = "steelblue", alpha = 0.4, colour = "steelblue") +
    geom_vline(
      data = truths,
      aes(xintercept = true_value),
      colour = "firebrick", linewidth = 0.9, linetype = "dashed"
    ) +
    facet_wrap(~label, )
    labs(
      title = paste0("Posterior recovery — ", fam_name),
      subtitle = "Blue = posterior density, red dashed = true (simulated) value",
      x = "value", y = "density"
    ) +
    theme_minimal(base_size = 12)
}

p_betamat   <- plot_family("betamat")
p_h         <- plot_family("h")
p_rho       <- plot_family("rho_z")
p_p_village <- plot_family("p_village")

print(p_betamat)
print(p_h)
print(p_rho)
print(p_p_village)




