library(MCMCpack)
library(tidyverse)
library(VGAM) 
library(randomForest)

set.seed(123)


n_loci <- 2
source("C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/R/create_dummy_matrices.R")
genotype_combination <- create_dummy_matrices(n_loci)
# defining the matrix name for the hand of genotype
L <- genotype_combination$left
L
R <- genotype_combination$right
L[3, 2]
G <- nrow(L)
B <- ncol(L)
#### GENOTYPE 
# initial allele frequency
p <- runif(B, 0,1)
Ttime <- 5
# ++++++++hardy–Weinberg genotype probabilities+++++++++++

probability_genotype_fast <- function(p, L, R) {
  prob_left  <- sweep(L, 2, 1 - p, "*") + sweep(1 - L, 2, p, "*")
  prob_right <- sweep(R, 2, 1 - p, "*") + sweep(1 - R, 2, p, "*")
  dup <- 1 + L - R # hetero to not having the duplicate RS
  # F <- matrix(NA_real_, nrow = G, ncol=B)
  F <- prob_left * prob_right * dup
  # z <- array(NA_real_, dim = G)
  z <- apply(F, 1, prod)
  
  return(z)
}


# +++++++++fitness model per genotype+++++++++++++++++++
polygenic_multilocus_next_step <- function(z, w, h, L, R) {
  G <- nrow(L)
  B <- ncol(L)
  Gw <- matrix(NA_real_, nrow = G, ncol = B)
  r_vec <- array(NA_real_, dim = G)
  for (lo in seq_len(B)) {
    SS <- (L[, lo] == 1L) & (R[, lo] == 1L)
    RR <- (L[, lo] == 0L) & (R[, lo] == 0L)
    SR <- (L[, lo] == 1L) & (R[, lo] == 0L)
    Gw[, lo] <- 1 * SS +
      w[lo] * RR +
      (h[lo] * w[lo] + (1 - h[lo])) * SR
  }
  r_vec <- apply(Gw, 1, prod)
  unpost <- z * r_vec
  genotype_post <- unpost / sum(unpost)
  return(genotype_post)
}

# +++++++++++++++++parameters+++++++++++++++++++++
# time steps
Tmax <- 500
# mosquitoes sampled per generation
M_z <- 200 
# overdispersion for Dirichlet–Multinomial
rho_z <- 0.05

# convert to 'burstiness' for dirichlet distribution
# https://en.wikipedia.org/wiki/Dirichlet-multinomial_distribution#Matrix_notation
alpha0_z <- (1 - rho_z ^2 ) / (rho_z ^ 2)

# initial allele frequencies
p <- runif(B, 0.3, 0.7) 
w <- runif(B, 1, 1.4) # should be > 1
h <- runif(B, 0, 1)

# +++++++++++++Initialize+++++++++++++
# true genotype frequencies
Z_list <- list() 
Z_list[[1]] <- probability_genotype_fast(p = p, L = L, R = R)

# +++++++++++++simulation over time+++++++++++++++++++++
for (t in 2:Tmax) {
  # Update genotype frequencies deterministically (selection)
  Z_true <- polygenic_multilocus_next_step(
    z = Z_list[[t-1]],
    w = w,
    h = h,
    L = L,
    R = R
  )
  Z_list[[t]] <- Z_true
  
  # # dirichlet–Multinomial likelihood
  # alpha <- pmax(Z_true * (1 - rho_z) / rho_z, 1e-8)
  # Z_disp <- as.numeric(rdirichlet(1, alpha))
  # Z_disp <- Z_disp / sum(Z_disp)
  # N_obs <- as.vector(rmultinom(1, M_z, Z_disp))
  # N_list[[t]] <- N_obs
  
}

# observed counts (Dirichlet–Multinomial)
N_list <- list()


for (t in 1:Tmax) {
  # dirichlet–Multinomial likelihood
  alpha <- Z_list[[t]] * alpha0_z
  # alpha <- pmax(Z_list[[t]] * (1 - rho_z) / rho_z, 1e-8)
  Z_disp <- as.numeric(rdirichlet(1, alpha))
  # Z_disp <- Z_disp / sum(Z_disp)
  N_obs <- as.vector(rmultinom(1, M_z, Z_disp))
  N_list[[t]] <- N_obs
}


# +++++++++++convert lists to matrices for plotting++++++++++
Z_mat <- do.call(cbind, Z_list)
N_mat <- do.call(cbind, N_list)
# ?do.call
# ++++++++plot genotype frequencies over time+++++++++++++++++++++
matplot(t(Z_mat), type = "l", lwd = 2, lty = 1,
        xlab = "Generation", ylab = "Genotype frequency",
        main = "Evolution of multilocus genotype frequencies (true)")

# ++++++++plot observed counts over time++++++++++++++
matplot(t(N_mat), type = "l", lwd = 2, lty = 1,
        xlab = "Generation", ylab = "Observed genotype counts",
        main = "Observed counts (Dirichlet–Multinomial)")



# ----------------------------------------------------
# GENOTYPE TO PHENOTYPE
# Without taking account of the different type of epistasis, here we consider epsilon as the modifier 
# on how gene interact to each other.
# polygenic 
# Compute locus-level genotype susceptibility
compute_Ugc <- function(L, R, w, h) {
  G <- nrow(L)
  B <- ncol(L)
  f <- matrix(NA, nrow = G, ncol = B)
  
  for (l in 1:B) {
    SS <- (L[, l] == 1 & R[, l] == 1)
    RR <- (L[, l] == 0 & R[, l] == 0)
    SR <- (L[, l] == 1 & R[, l] == 0) | (L[, l] == 0 & R[, l] == 1)
    
    f[, l] <- 1 * SS +
      w[l] * RR +
      (h[l] * w[l] + (1 - h[l])) * SR
  }
  
  return(f)  # Ugc
}


compute_Ustar <- function(Ugc, theta, type = "additive", epsilon = NULL) {
  
  # Compute base contributions
  U_add <- rowSums(Ugc)
  U_mult <- apply(Ugc, 1, prod)
  
  if (type == "additive") {
    U_star <- U_add
    
  } else if (type == "multiplicative") {
    U_star <- U_mult
    
  } else if (type == "epistatic") {
    if (is.null(epsilon)) stop("epsilon required")
    
    # Pairwise epistatic interactions
    pairwise_epi <- rowSums(Ugc %*% epsilon * Ugc)
    
    # Combine contributions
    U_star <- U_add + U_mult + pairwise_epi
  }
  
  # Scale by sum of theta
  U_star <- U_star / sum(theta)
  
  return(U_star)
}



# Compute mortality probability for a genotype
compute_p_died <- function(U_star) {
  p <- 1 - U_star
  p[p < 0] <- 0
  p[p > 1] <- 1
  return(p)
}

# Beta-Binomial simulation for bioassay
simulate_beta_binomial <- function(p, n, phi = 20) {
  alpha <- p * phi
  beta  <- (1 - p) * phi
  p_sample <- rbeta(length(p), alpha, beta)
  y <- rbinom(length(p), size = n, prob = p_sample)
  return(list(p_sample = p_sample, y = y))
}

# Parameters
theta <- runif(G, 0.5, 1)
epsilon <- matrix(runif(B^2, -0.2, 0.2), B, B)
n_tested <- 100  # mosquitoes per bioassay

# Compute locus-level susceptibility
Ugc <- compute_Ugc(L, R, w, h)

# Individual U_star types
U_add <- compute_Ustar(Ugc, theta, "additive")
U_mult <- compute_Ustar(Ugc, theta, "multiplicative")
U_epi <- compute_Ustar(Ugc, theta, "epistatic", epsilon)

# Mortality probabilities per genotype
p_add <- compute_p_died(U_add)
p_mult <- compute_p_died(U_mult)
p_epi <- compute_p_died(U_epi)

# Combined U_star with weights alpha, beta, gamma
alpha <- 0.33; beta <- 0.33; gamma <- 0.34
U_combined <- alpha * U_add + beta * U_mult + gamma * U_epi
p_combined <- compute_p_died(U_combined)

# Simulate bioassay deaths
sim_add <- simulate_beta_binomial(p_add, n_tested)
sim_mult <- simulate_beta_binomial(p_mult, n_tested)
sim_epi <- simulate_beta_binomial(p_epi, n_tested)
sim_combined <- simulate_beta_binomial(p_combined, n_tested)

# Save results
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
df


### Random Forest genotype to phenotype
# Using Z_mat data from the genotype
#n <- nrow(Z_mat)
p <- ncol(Z_mat) # number of genotype in depend of Z_mat
colnames(Z_mat) <- paste0("g", 1:p)

# coefficients for each genotype, this will determine how much the genotype
# is contributing to the phenotype
alpha <- runif(p, 0.1, 1)


# convert Z_mat as a dataframe to use for the model
Z_df <- as.data.frame(Z_mat)
genotype_cols <- select(Z_df, starts_with("g"))


# function to compute our phenotype from the genotype
# here we say that alpha will give the importance of each genotype
# and this (pheno - min(pheno)) / (max(pheno) - min(pheno))  to remain phenotype between 0 and 1
pheno <- function(Z_mat, effect_type = "additive", alpha = NULL, epsilon = NULL) {
  p <- ncol(Z_mat)
  if (is.null(alpha)) alpha <- runif(p, 0.1, 1)  # genotype weights
  
  if (effect_type == "additive") {
    pheno <- as.vector(as.matrix(Z_mat) %*% alpha)
  } else if (effect_type == "epistatic") {
    if (is.null(epsilon)) epsilon <- matrix(runif(p^2, -0.2, 0.2), p, p)
    pheno <- as.vector(as.matrix(Z_mat) %*% alpha + rowSums((Z_mat %*% epsilon) * Z_mat))
  }
  
  # normalize 0-1
  pheno <- (pheno - min(pheno)) / (max(pheno) - min(pheno))
  return(pheno)
}

pheno <- simulate_pheno(Z_mat, effect_type = "additive")

pheno <- compute_phenotype(genotype_cols, alpha)

# saving the data for random forest
df <- data.frame(pheno, genotype_cols)

# Random Forest linear regression
# we consider that each genotype are independent so they could follow the equation:
# pheno = alpha * g1 + alpha * g2 + ....... + alpha gn
rf_model <- randomForest(pheno ~ ., data=df,
                         ntree=500,      # tree number
                         mtry=2,         # splitting the data to 2 for example
                         importance=TRUE) # calculate the importance of each variable

# model 
print(rf_model)

# Importance of each genotype
importance(rf_model)

# Predicted phenotype
pred <- predict(rf_model, newdata=df)

# Verify the correaltion between predicted and observed
cor(pred, df$pheno)



#-------------------------------------------------------
# ALLELE FREQUENCY
# function
allele_frequency_next_step <- function(genotype_next, L, R) {
  G <- nrow(L); B <- ncol(L)
  p_next <- numeric(B)
  for (lo in seq_len(B)) {
    a_lo <- 2L - L[, lo] - R[, lo]   # SS=0, SR=1, RR=2
    p_next[lo] <- 0.5 * sum(genotype_next * a_lo)
  }
  p_next
}

probability_genotype_fast <- function(p, L, R) {
  prob_left  <- sweep(L, 2, 1 - p, "*") + sweep(1 - L, 2, p, "*")
  prob_right <- sweep(R, 2, 1 - p, "*") + sweep(1 - R, 2, p, "*")
  dup <- 1 + L - R
  F <- prob_left * prob_right * dup
  z <- apply(F, 1, prod)
  z / sum(z)
}

polygenic_multilocus_next_step <- function(z, w, h, L, R) {
  G <- nrow(L); B <- ncol(L)
  Gw <- matrix(NA_real_, nrow = G, ncol = B)
  for (lo in seq_len(B)) {
    SS <- (L[, lo] == 1L) & (R[, lo] == 1L)
    RR <- (L[, lo] == 0L) & (R[, lo] == 0L)
    SR <- (L[, lo] == 1L) & (R[, lo] == 0L)
    Gw[, lo] <- 1 * SS + w[lo] * RR + (h[lo] * w[lo] + (1 - h[lo])) * SR
  }
  r_vec <- apply(Gw, 1, prod)
  genotype_post <- z * r_vec
  genotype_post / sum(genotype_post)
}

# Sample observed allele counts using Beta-Binomial
sample_allele_counts <- function(p_next, M, rho) {
  B <- length(p_next)
  N_obs <- numeric(B)
  for(lo in seq_len(B)) {
    phi <- (1 - rho)/rho
    alpha <- p_next[lo] * phi
    beta  <- (1 - p_next[lo]) * phi
    N_obs[lo] <- rbetabinom(1, M, alpha, beta)
  }
  N_obs 
}

# Parameters 
set.seed(123)
Tmax <- 50      # generations
M <- 200        # mosquitoes sampled per locus
rho <- 0.05     # overdispersion
n_loci <- 2
p <- runif(B, 0.3, 0.7)
z <- probability_genotype_fast(p, L, R)
w <- runif(B, 0.8, 1)
h <- runif(B, 0, 1)

#  Simulation
allele_freq_list <- list()

allele_count_list <- list()

allele_freq_list[[1]] <- allele_frequency_next_step(z, L, R)
allele_count_list[[1]] <- sample_allele_counts(allele_freq_list[[1]], M, rho)

for (t in 2:Tmax) {
  z <- polygenic_multilocus_next_step(z, w, h, L, R)
  p_next <- allele_frequency_next_step(z, L, R)
  
  allele_freq_list[[t]] <- p_next
  allele_count_list[[t]] <- sample_allele_counts(p_next, M, rho)
}

# Convert to tidy dataframe
allele_mat <- do.call(rbind, allele_freq_list)
allele_count_mat <- do.call(rbind, allele_count_list)

allele_df <- as.data.frame(allele_mat)
allele_count_df <- as.data.frame(allele_count_mat)
colnames(allele_df) <- paste0("Locus_", 1:B)
colnames(allele_count_df) <- paste0("Locus_", 1:B)

allele_df$Generation <- 1:Tmax
allele_count_df$Generation <- 1:Tmax

allele_long <- allele_df %>%
  pivot_longer(cols = starts_with("Locus_"),
               names_to = "Locus",
               values_to = "Allele_Frequency")

head(allele_long, 20)





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##### GENOTYPE FREQUENCY
#++++ getting dataframe+++++++++++
# create genotype labels
# genotype_labels <- function(L, R) {
#   G <- nrow(L)
#   B <- ncol(L)
#   labels <- character(G)
#   
#   for (g in 1:G) {
#     loci <- character(B)
#     for (lo in 1:B) {
#       if (L[g, lo] == 1 && R[g, lo] == 1) {
#         loci[lo] <- "SS"
#       } else if (L[g, lo] == 0 && R[g, lo] == 0) {
#         loci[lo] <- "RR"
#       } else {
#         loci[lo] <- "SR"
#       }
#     }
#     labels[g] <- paste(loci, collapse = "-")
#   }
#   
#   return(labels)
# }
# 
# labels <- genotype_labels(L, R)
# 
# # convert Z_mat (true frequencies) to long dataframe
# Z_df <- as.data.frame(Z_mat)
# colnames(Z_df) <- paste0("Gen_", 1:ncol(Z_df))
# Z_df$Genotype <- labels
# 
# Z_long <- Z_df %>%
#   pivot_longer(cols = starts_with("Gen_"),
#                names_to = "Generation",
#                values_to = "Frequency") %>%
#   mutate(Generation = as.integer(gsub("Gen_", "", Generation)))
# 
# # optionally: convert observed counts N_mat to dataframe
# N_df <- as.data.frame(N_mat)
# colnames(N_df) <- paste0("Gen_", 1:ncol(N_df))
# N_df$Genotype <- labels
# 
# N_long <- N_df %>%
#   pivot_longer(cols = starts_with("Gen_"),
#                names_to = "Generation",
#                values_to = "Count") %>%
#   mutate(Generation = as.integer(gsub("Gen_", "", Generation)))
# 
# # preview
# head(Z_long)
# head(N_long)
# 
# ### genotype per SS RR RS
# ################ get a dataframe
# # parameters 
# # generations
# Tmax <- ncol(Z_mat) 
# # number of loci
# B <- ncol(L) 
# # number of multilocus genotypes
# G <- nrow(L)         
# genotypes <- c("SS", "SR", "RR")
# 
# # function to get per-locus genotype indices 
# get_indices <- function(L, R, locus, gt) {
#   if(gt == "SS") return(which(L[, locus] == 1 & R[, locus] == 1))
#   if(gt == "SR") return(which((L[, locus] == 1 & R[, locus] == 0) | (L[, locus] == 0 & R[, locus] == 1)))
#   if(gt == "RR") return(which(L[, locus] == 0 & R[, locus] == 0))
# }
# 
# # create expand.grid for all combinations 
# df <- expand.grid(
#   Generation = 1:Tmax,
#   Locus = 1:B,
#   Genotype = genotypes
# )
# 
# # compute per-locus genotype frequencies
# df$Frequency <- mapply(function(gen, loc, gt){
#   idx <- get_indices(L, R, loc, gt)
#   sum(Z_mat[idx, gen])
# }, df$Generation, df$Locus, df$Genotype)
# 
# # convert to tibble
# df_tibble <- as_tibble(df)
# 
# # preview
# df_tibble %>% arrange(Locus, Generation)


## We assume each batch under insecticide c has its own true susceptibility θ_{batch,c}.
## Heterogeneity across batches is modeled by θ_{batch,c} ~ Beta(α_c, β_c) with mean p_c
## and overdispersion ρ_c. 
## We then sample deaths within the batch as Binomial(N, θ_{batch,c}).
# rbetabinom1 <- function(N, p, rho) {
#   stopifnot(N >= 0, p > 0, p < 1, rho >= 0, rho < 1)
#   if (rho == 0) return(stats::rbinom(1, size = N, prob = p))
#   ab <- rho_to_ab(p, rho)
#   theta <- stats::rbeta(1, ab$alpha, ab$beta)
#   stats::rbinom(1, size = N, prob = theta)
# }

### the average phenotype
## - L,R encode each genotype’s alleles per locus (0/1). With locus effects θ (by insecticide c),
##   dominance h, and a baseline susceptible mortality v_s[c], we build the per-genotype
##   susceptible probability:
##       U_{g,c} = v_s[c] * ∏_b u_{g b c}
##   where u_{g b c} = 1 (SS), θ_{b c} (RR), or (1−h_b)+h_b*θ_{b c} (heterozygote; here via SR).
## - Z[,t] are genotype frequencies at time t (sum to 1). The overall susceptible phenotype
##   proportion at time t under insecticide c is the mixture:
##       Q*_{t,c} = Σ_g Z_{g,t} * U_{g,c}.
## This function computes U_{g,c} multiplicatively across loci, then averages over genotypes
## with Z to return Q* over time (rows) and insecticides (cols).
# priors
# C <- 2
# theta <- matrix(runif(B * C, 0, 1), nrow = B, ncol =C)
# h <- array(runif(B, 0,1))
# v_s <- array(runif(C, 0,1))
# # function to compute the average phenotype
# compute_Qstar <- function(theta, h, v_s) {
#   stopifnot(is.matrix(L), is.matrix(R), all(dim(L) == dim(R)))
#   G <- nrow(L); B <- ncol(L)
#   stopifnot(is.matrix(theta) && nrow(theta) == B)
#   stopifnot(length(h) == B, length(v_s) == C)
#   # stopifnot(is.array(z) && nrow(z) == G)
#   
#   # per-genotype x insecticide utility U_gc (G x C)
#   SS <- L * R
#   RR <- (1L - L) * (1L - R)
#   SR <- L * (1L - R)
#   
#   theta_3 <- array(theta, dim = c(1, B, C)); 
#   theta_3 <- theta_3[rep(1, G), , , drop = FALSE]
#   h_3     <- array(h,     dim = c(1, B, 1)); 
#   h_3     <- h_3[rep(1, G), , rep(1, C), drop = FALSE]
#   SS_3 <- array(SS, dim = c(G, B, 1))[, , rep(1, C), drop = FALSE]
#   RR_3 <- array(RR, dim = c(G, B, 1))[, , rep(1, C), drop = FALSE]
#   SR_3 <- array(SR, dim = c(G, B, 1))[, , rep(1, C), drop = FALSE]
#   # to avoid that u_3 value will give 0 or negative, we use the .Machine$double.eps
#   eps  <- .Machine$double.eps
#   u_3  <- SS_3 + RR_3 * theta_3 + SR_3 * (h_3 * theta_3 + (1 - h_3))
#   U_gc <- exp(apply(log(pmax(u_3, eps)), c(1, 3), sum))  
#   # for each c, multiply the column by v_s
#   U_gc <- sweep(U_gc, 2, v_s, "*")                       
#   Q_star <- array(as.numeric(crossprod(z, U_gc)))
#   pmin(pmax(Q_star, 0), 1)
#   return(Q_star)
# }
# compute_Qstar(theta, h, v_s)




