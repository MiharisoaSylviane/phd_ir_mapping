
# this code contain the main script for simulate our model

#load library
library(MCMCpack)
library(tidyverse)
library(VGAM) 
library(reshape2)
#library(randomForest)
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# library(remotes)
# run function : function_geno_pheno.R and create_dummy_matrices.R

# source functions
#getwd() # to see the directory
setwd("C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/R")
list.files() # check all files

n_loci <- 5
# call function to use the simulation script
source("create_dummy_matrices.R")
genotype_combination <- create_dummy_matrices(n_loci)
L <- genotype_combination$left
L
R <- genotype_combination$right

# Parameters
Tmax <- 20
M_z <- 200
rho_z <- 0.05
p <- runif(n_loci, 0.3, 0.7)
w <- runif(n_loci, 0.8, 1.4)
h <- runif(n_loci, 0, 1)
set.seed(123)
# Create genotype matrices
G <- nrow(L)
B <- ncol(L)

source("function_geno_pheno.R")
# Initialize genotype frequencies
Z_list <- list()
Z_list[[1]] <- probability_genotype_fast(p, L, R)


# Simulation over time
for (t in 2:Tmax) {
  Z_list[[t]] <- polygenic_multilocus_next_step(Z_list[[t-1]], w, h, L, R)
}

# Sample observed counts
N_list <- lapply(Z_list, function(z) sample_genotype_counts(z, M_z, rho_z))

# Convert to matrices for plotting
Z_mat <- do.call(cbind, Z_list)
N_mat <- do.call(cbind, N_list)


# Plot
matplot(t(Z_mat), type="l", lwd=2, lty=1,
        xlab="Generation", ylab="Genotype frequency",
        main="Evolution of genotype frequencies (true)")

matplot(t(N_mat), type="l", lwd=2, lty=1,
        xlab="Generation", ylab="Observed counts",
        main="Observed counts (Dirichlet-Multinomial)")


# Genotype → Phenotype

theta <- runif(G, 0.5, 1)
epsilon <- matrix(runif(B^2, -0.2, 0.2), B, B)
n_tested <- 100

Ugc <- compute_Ugc(L, R, w, h)

U_add  <- compute_Ustar(Ugc, theta, "additive")
U_mult <- compute_Ustar(Ugc, theta, "multiplicative")
U_epi  <- compute_Ustar(Ugc, theta, "epistatic", epsilon)

p_add  <- compute_p_died(U_add)
p_mult <- compute_p_died(U_mult)
p_epi  <- compute_p_died(U_epi)

alpha <- 0.33; beta <- 0.33; gamma <- 0.34
U_combined <- alpha * U_add + beta * U_mult + gamma * U_epi
p_combined <- compute_p_died(U_combined)


# Bioassay simulation
sim_add <- simulate_beta_binomial(p_add, n_tested)
sim_mult <- simulate_beta_binomial(p_mult, n_tested)
sim_epi <- simulate_beta_binomial(p_epi, n_tested)
sim_combined <- simulate_beta_binomial(p_combined, n_tested)


# saving data

df <- data.frame(
  Genotype = 1:G,
  p_add = p_add,
  p_mult = p_mult,
  p_epi = p_epi,
  p_combined = p_combined,
  dead_add = sim_add$y,
  dead_mult = sim_mult$y,
  dead_epi = sim_epi$y,
  dead_combined = sim_combined$y
)

print(df)


# PHENOTYPE OVER TIME

# Compute genotype-level phenotype (fixed over time)

Ugc <- compute_Ugc(L, R, w, h)
U_star <- compute_Ustar(Ugc, theta, type = "epistatic", epsilon = epsilon)

U_star <- pmin(U_star, 1) # so U_star will belong to 0 to 1

# Compute population phenotype over time
Tmax <- length(Z_list)
pheno_time <- numeric(Tmax)

for (t in 1:Tmax) {
  freq <- Z_list[[t]]              # genotype frequencies at time t
  pheno_time[t] <- sum(freq * U_star)
}

# plot phenotype over time

df_pheno <- data.frame(
  Time = 1:Tmax,
  Phenotype = pheno_time
)

ggplot(df_pheno, aes(x = Time, y = Phenotype)) +
  geom_line(color = "darkgreen", size = 1.2) +
  labs(
    title = "Population phenotype over time",
    x = "Generation",
    y = "Mean phenotype (resistance)"
  ) +
  theme_minimal()
