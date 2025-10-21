library(MCMCpack) 
library(tidyverse)
library(VGAM) 

set.seed(123)

# +++++++++setup+++++++++++++++++++++++++++++
n_loci <- 4
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
# ++++++++hardy–Weinberg genotype probabilities+++++++++++
probability_genotype_fast <- function(p, L, R) {
  prob_left  <- sweep(L, 2, 1 - p, "*") + sweep(1 - L, 2, p, "*")
  prob_right <- sweep(R, 2, 1 - p, "*") + sweep(1 - R, 2, p, "*")
  dup <- 1 + L - R
  F <- prob_left * prob_right * dup
  z <- apply(F, 1, prod)
  z <- z / sum(z)
  return(z)
}

# +++++++++fitness model per genotype+++++++++++++++++++
polygenic_multilocus_next_step <- function(z, w, h, L, R) {
  Gw <- matrix(NA_real_, nrow = G, ncol = B)
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
Tmax <- 50 
# mosquitoes sampled per generation
M_z <- 200 
# overdispersion for Dirichlet–Multinomial
rho_z <- 0.05   
# initial allele frequencies
p <- runif(B, 0.3, 0.7) 
w <- runif(B, 0.8, 1)
h <- runif(B, 0, 1)

# +++++++++++++Initialize+++++++++++++
# true genotype frequencies
Z_list <- list() 
# observed counts (Dirichlet–Multinomial)
N_list <- list()    
Z_list[[1]] <- probability_genotype_fast(p, L, R)

# +++++++++++++simulation over time+++++++++++++++++++++
for (t in 2:Tmax) {
  # Update genotype frequencies deterministically (selection)
  Z_true <- polygenic_multilocus_next_step(Z_list[[t-1]], w, h, L, R)
  Z_list[[t]] <- Z_true
  
  # dirichlet–Multinomial likelihood
  alpha <- pmax(Z_true * (1 - rho_z) / rho_z, 1e-8)
  Z_disp <- as.numeric(rdirichlet(1, alpha))
  Z_disp <- Z_disp / sum(Z_disp)
  N_obs <- as.vector(rmultinom(1, M_z, Z_disp))
  N_list[[t]] <- N_obs
  
}

# +++++++++++convert lists to matrices for plotting++++++++++
Z_mat <- do.call(cbind, Z_list)
N_mat <- do.call(cbind, N_list)

# ++++++++plot genotype frequencies over time+++++++++++++++++++++
matplot(t(Z_mat), type = "l", lwd = 2, lty = 1,
        xlab = "Generation", ylab = "Genotype frequency",
        main = "Evolution of multilocus genotype frequencies (true)")

# ++++++++plot observed counts over time++++++++++++++
matplot(t(N_mat), type = "l", lwd = 2, lty = 1,
        xlab = "Generation", ylab = "Observed genotype counts",
        main = "Observed counts (Dirichlet–Multinomial)")

#### GENOTYPE FREQUENCY
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++ getting dataframe+++++++++++
# create genotype labels
genotype_labels <- function(L, R) {
  G <- nrow(L)
  B <- ncol(L)
  labels <- character(G)
  
  for (g in 1:G) {
    loci <- character(B)
    for (lo in 1:B) {
      if (L[g, lo] == 1 && R[g, lo] == 1) {
        loci[lo] <- "SS"
      } else if (L[g, lo] == 0 && R[g, lo] == 0) {
        loci[lo] <- "RR"
      } else {
        loci[lo] <- "SR"
      }
    }
    labels[g] <- paste(loci, collapse = "-")
  }
  
  return(labels)
}

labels <- genotype_labels(L, R)

# convert Z_mat (true frequencies) to long dataframe
Z_df <- as.data.frame(Z_mat)
colnames(Z_df) <- paste0("Gen_", 1:ncol(Z_df))
Z_df$Genotype <- labels

Z_long <- Z_df %>%
  pivot_longer(cols = starts_with("Gen_"),
               names_to = "Generation",
               values_to = "Frequency") %>%
  mutate(Generation = as.integer(gsub("Gen_", "", Generation)))

# optionally: convert observed counts N_mat to dataframe
N_df <- as.data.frame(N_mat)
colnames(N_df) <- paste0("Gen_", 1:ncol(N_df))
N_df$Genotype <- labels

N_long <- N_df %>%
  pivot_longer(cols = starts_with("Gen_"),
               names_to = "Generation",
               values_to = "Count") %>%
  mutate(Generation = as.integer(gsub("Gen_", "", Generation)))

# preview
head(Z_long)
head(N_long)

### genotype per SS RR RS
################ get a dataframe
# parameters 
# generations
Tmax <- ncol(Z_mat) 
# number of loci
B <- ncol(L) 
# number of multilocus genotypes
G <- nrow(L)         
genotypes <- c("SS", "SR", "RR")

# function to get per-locus genotype indices 
get_indices <- function(L, R, locus, gt) {
  if(gt == "SS") return(which(L[, locus] == 1 & R[, locus] == 1))
  if(gt == "SR") return(which((L[, locus] == 1 & R[, locus] == 0) | (L[, locus] == 0 & R[, locus] == 1)))
  if(gt == "RR") return(which(L[, locus] == 0 & R[, locus] == 0))
}

# create expand.grid for all combinations 
df <- expand.grid(
  Generation = 1:Tmax,
  Locus = 1:B,
  Genotype = genotypes
)

# compute per-locus genotype frequencies
df$Frequency <- mapply(function(gen, loc, gt){
  idx <- get_indices(L, R, loc, gt)
  sum(Z_mat[idx, gen])
}, df$Generation, df$Locus, df$Genotype)

# convert to tibble
df_tibble <- as_tibble(df)

# preview
df_tibble %>% arrange(Locus, Generation)


#### ALLELE FREQUENCY
# ----------------- Functions -----------------
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

# ----------------- Parameters -----------------
set.seed(123)
Tmax <- 50      # generations
M <- 200        # mosquitoes sampled per locus
rho <- 0.05     # overdispersion
n_loci <- 4
p <- runif(B, 0.3, 0.7)
z <- probability_genotype_fast(p, L, R)
w <- runif(B, 0.8, 1)
h <- runif(B, 0, 1)

# ----------------- Simulation -----------------
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

# ----------------- Convert to tidy dataframe -----------------
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

# ----------------- Preview -----------------
head(allele_long, 20)

rho_to_ab <- function(p, rho) {
  stopifnot(all(p > 0 & p < 1), rho >= 0, rho < 1)
  if (rho == 0) return(list(alpha = Inf, beta = Inf))  
  m <- 1 / rho - 1
  list(alpha = p * m, beta = (1 - p) * m)
}

## We assume each batch under insecticide c has its own true susceptibility θ_{batch,c}.
## Heterogeneity across batches is modeled by θ_{batch,c} ~ Beta(α_c, β_c) with mean p_c
## and overdispersion ρ_c. 
## We then sample deaths within the batch as Binomial(N, θ_{batch,c}).
rbetabinom1 <- function(N, p, rho) {
  stopifnot(N >= 0, p > 0, p < 1, rho >= 0, rho < 1)
  if (rho == 0) return(stats::rbinom(1, size = N, prob = p))
  ab <- rho_to_ab(p, rho)
  theta <- stats::rbeta(1, ab$alpha, ab$beta)
  stats::rbinom(1, size = N, prob = theta)
}


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
compute_Qstar <- function(L, R, Z, theta, h, v_s) {
  stopifnot(is.matrix(L), is.matrix(R), all(dim(L) == dim(R)))
  G <- nrow(L); B <- ncol(L)
  stopifnot(is.matrix(theta) && nrow(theta) == B)
  C <- ncol(theta); stopifnot(length(h) == B, length(v_s) == C)
  stopifnot(is.matrix(Z) && nrow(Z) == G)
  
  # per-genotype x insecticide utility U_gc (G x C)
  SS <- L * R
  RR <- (1L - L) * (1L - R)
  SR <- L * (1L - R)
  
  theta_3 <- array(theta, dim = c(1, B, C)); theta_3 <- theta_3[rep(1, G), , , drop = FALSE]
  h_3     <- array(h,     dim = c(1, B, 1)); h_3     <- h_3[rep(1, G), , rep(1, C), drop = FALSE]
  SS_3 <- array(SS, dim = c(G, B, 1))[, , rep(1, C), drop = FALSE]
  RR_3 <- array(RR, dim = c(G, B, 1))[, , rep(1, C), drop = FALSE]
  SR_3 <- array(SR, dim = c(G, B, 1))[, , rep(1, C), drop = FALSE]
  
  eps  <- .Machine$double.eps
  u_3  <- SS_3 + RR_3 * theta_3 + SR_3 * (h_3 * theta_3 + (1 - h_3))
  U_gc <- exp(apply(log(pmax(u_3, eps)), c(1, 3), sum))  
  U_gc <- sweep(U_gc, 2, v_s, "*")                       
  
  # combine with genotype distribution over time
  Ttime <- ncol(Z)
  Q <- matrix(NA_real_, nrow = Ttime, ncol = C)
  for (t in seq_len(Ttime)) {
    z_t <- Z[, t]                    # length G, sums to 1
    Q[t, ] <- as.numeric(crossprod(z_t, U_gc))
  }
  pmin(pmax(Q, 0), 1)
}



simulate_fake_data <- function(
    G , B , Ttime , C ,
    N_per_cell , rho ,
    v_s_decay       
) {
  set.seed(123)
  G <- as.integer(G)[1]; B <- as.integer(B)[1]
  Ttime <- as.integer(Ttime)[1]; C <- as.integer(C)[1]
  
  # Genotype encodings (G x B)
  L <- matrix(sample(0:1, G * B, TRUE), nrow = G)
  R <- matrix(sample(0:1, G * B, TRUE), nrow = G)
  
  # Genotype frequencies Z (G x T)
  Z <- matrix(NA_real_, nrow = G, ncol = Ttime)
  for (t in seq_len(Ttime)) { v <- stats::runif(G); Z[, t] <- v / sum(v) }
  
  # Parameters
  theta <- matrix(stats::runif(B * C, 0.3, 0.9), nrow = B, ncol = C)
  h     <- stats::runif(B, 0, 1)
  v_s0  <- stats::runif(C, 0.6, 1.0)
  
  # the baseline of percent mortality
  v_s_t <- matrix(v_s0, nrow = Ttime, ncol = length(v_s0), byrow = TRUE)
  
  # Compute Q* at each time with that time's v_s
  Qstar <- matrix(NA_real_, nrow = Ttime, ncol = C)
  for (t in seq_len(Ttime)) {
    Qstar[t, ] <- compute_Qstar(L, R, Z[, t, drop = FALSE], theta, h, v_s_t[t, ])
  }
  
  # Tidy data (Time x Insecticide), plus Beta–Binomial observations
  df <- expand.grid(Time = seq_len(Ttime), Insecticide = seq_len(C), KEEP.OUT.ATTRS = FALSE)
  df$Qstar <- as.vector(Qstar)
  df$InsecticideName <- paste0("Ins_", df$Insecticide)
  df$N <- N_per_cell
  df$Deaths <- mapply(function(n, q) rbetabinom1(n, q, rho), df$N, df$Qstar)
  df$Observed <- df$Deaths / df$N
  
  list(
    data = tibble::as_tibble(df),  # tidy frame
    Qstar = Qstar,                 # T x C
    Z = Z, L = L, R = R,
    theta = theta, h = h, v_s0 = v_s0, v_s_decay = v_s_decay,
    params = list(G = G, B = B, Ttime = Ttime, C = C,
                  N_per_cell = N_per_cell, rho = rho)
  )
}


sim <- simulate_fake_data(
  G = 9, B = 2, Ttime = 8, C = 2,
  N_per_cell = 120, rho = 0.03, 
  v_s_decay = 0.97   
)
df <- sim$data
head(df)





























