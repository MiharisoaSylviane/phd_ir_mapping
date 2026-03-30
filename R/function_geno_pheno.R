#library(randomForest)

# Genotype probability (Hardy-Weinberg)
probability_genotype_fast <- function(p, L, R) {
  prob_left  <- sweep(L, 2, 1 - p, "*") + sweep(1 - L, 2, p, "*")
  prob_right <- sweep(R, 2, 1 - p, "*") + sweep(1 - R, 2, p, "*")
  dup <- 1 + L - R
  F <- prob_left * prob_right * dup
  z <- apply(F, 1, prod)
  z / sum(z)
}

# Multilocus polygenic selection step
polygenic_multilocus_next_step <- function(z, w, h, L, R) {
  G <- nrow(L); B <- ncol(L)
  Gw <- matrix(NA_real_, nrow = G, ncol = B)
  for (lo in seq_len(B)) {
    SS <- (L[, lo] == 1L) & (R[, lo] == 1L)
    RR <- (L[, lo] == 0L) & (R[, lo] == 0L)
    SR <- (L[, lo] != R[, lo])
    Gw[, lo] <- 1 * SS + w[lo] * RR + (h[lo] * w[lo] + (1 - h[lo])) * SR
  }
  r_vec <- apply(Gw, 1, prod)
  genotype_post <- z * r_vec
  genotype_post / sum(genotype_post)
}


# Allele frequency from genotype

allele_frequency_next_step <- function(genotype_next, L, R) {
  G <- nrow(L); B <- ncol(L)
  p_next <- numeric(B)
  for (lo in seq_len(B)) {
    a_lo <- 2L - L[, lo] - R[, lo]  # SS=0, SR=1, RR=2
    p_next[lo] <- 0.5 * sum(genotype_next * a_lo)
  }
  p_next
}


# Dirichlet-Multinomial sampling

sample_genotype_counts <- function(Z_true, M_z, rho_z) {
  alpha0_z <- (1 - rho_z^2)/rho_z^2
  alpha <- Z_true * alpha0_z
  Z_disp <- as.numeric(rdirichlet(1, alpha))
  as.vector(rmultinom(1, M_z, Z_disp))
}




# Phenotype function
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
  
  return(f)
}

compute_Ustar <- function(Ugc, theta, type = "additive", epsilon = NULL) {
  
  U_add <- rowSums(Ugc)
  U_mult <- apply(Ugc, 1, prod)
  
  if (type == "additive") {
    U_star <- U_add
    
  } else if (type == "multiplicative") {
    U_star <- U_mult
    
  } else if (type == "epistatic") {
    if (is.null(epsilon)) stop("epsilon required")
    pairwise_epi <- rowSums(Ugc %*% epsilon * Ugc)
    U_star <- U_add + U_mult + pairwise_epi
  }
  
  U_star <- U_star / sum(theta)
  return(U_star)
}

compute_p_died <- function(U_star) {
  p <- 1 - U_star
  p[p < 0] <- 0
  p[p > 1] <- 1
  return(p)
}

simulate_beta_binomial <- function(p, n, phi = 20) {
  alpha <- p * phi
  beta  <- (1 - p) * phi
  p_sample <- rbeta(length(p), alpha, beta)
  y <- rbinom(length(p), size = n, prob = p_sample)
  return(list(p_sample = p_sample, y = y))
}

pheno <- function(Z_mat, effect_type = "additive", alpha = NULL, epsilon = NULL) {
  p <- ncol(Z_mat)
  
  if (is.null(alpha)) alpha <- runif(p, 0.1, 1)
  
  if (effect_type == "additive") {
    pheno <- as.vector(as.matrix(Z_mat) %*% alpha)
    
  } else if (effect_type == "epistatic") {
    if (is.null(epsilon)) epsilon <- matrix(runif(p^2, -0.2, 0.2), p, p)
    pheno <- as.vector(as.matrix(Z_mat) %*% alpha +
                         rowSums((Z_mat %*% epsilon) * Z_mat))
  }
  
  pheno <- (pheno - min(pheno)) / (max(pheno) - min(pheno))
  return(pheno)
}
