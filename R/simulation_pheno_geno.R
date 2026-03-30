
# MAIN SIMULATION SCRIPT
# this code contain the main script for simulate our model

#load library
library(tidyverse)
library(MCMCpack)
library(randomForest)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
library(remotes)

# source functions
# source_url("https://raw.githubusercontent.com/MiharisoaSylviane/phd_ir_mapping/main/R/function_geno_pheno.R")
# source_url("https://raw.githubusercontent.com/MiharisoaSylviane/phd_ir_mapping/main/R/create_dummy_matrices.R")

#getwd()
setwd("C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/R")
list.files()
# call function to use the simulation script
source("function_geno_pheno.R")
source("create_dummy_matrices.R")

set.seed(123)

# Parameters
n_loci <- 2
Tmax <- 20
M_z <- 200
rho_z <- 0.05
p <- runif(n_loci, 0.3, 0.7)
w <- runif(n_loci, 0.8, 1.4)
h <- runif(n_loci, 0, 1)

# Create genotype matrices
genotype_combination <- create_dummy_matrices(n_loci)
L <- genotype_combination$left
R <- genotype_combination$right
G <- nrow(L)
B <- ncol(L)


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

# Simulation Script
# Run the simulation over all generations
pheno_res <- compute_survival_over_time(Z_list, L, R, w, h, theta, epsilon, effect_type = "additive")
# Population-level survival dataframe
df_pop <- data.frame(
  Time = 1:length(pheno_res$population),
  Survival = pheno_res$population
)

# Plot population-level survival
ggplot(df_pop, aes(x = Time, y = Survival)) +
  geom_line(color = "steelblue", size = 1) +
  labs(
    title = "Population-level survival probability over time",
    x = "Time (generations)",
    y = "Survival probability"
  ) +
  theme_minimal()

# Per-genotype survival dataframe

geno_df <- as.data.frame(pheno_res$genotype)
colnames(geno_df) <- paste0("G", 1:ncol(geno_df))
geno_df$Time <- 1:nrow(geno_df)

geno_long <- melt(
  geno_df,
  id.vars = "Time",
  variable.name = "Genotype",
  value.name = "Survival"
)

# Plot per-genotype survival
ggplot() +
  geom_line(data = geno_long, aes(x = Time, y = Survival, color = Genotype), alpha = 0.5) +
  geom_line(data = df_pop, aes(x = Time, y = Survival), color = "red", size = 1.2) +
  labs(
    title = "Evolution of mosquito survival probability over time",
    x = "Time (generations)",
    y = "Survival probability"
  ) +
  theme_minimal() +
  scale_color_viridis_d(option = "plasma") +
  theme(legend.position = "right")

