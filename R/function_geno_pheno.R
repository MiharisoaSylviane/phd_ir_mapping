#library(randomForest)
# intervention model
compute_w_greta <- function(betamat, X) {
  s <- exp(betamat %*% X)
  1 + s
}

## to get the selection pressure advantage of resistant allele
compute_s_greta <- function(betamat, X) {
  s <- exp(betamat %*% X)
  return(s)
}

# Genotype probability (Hardy-Weinberg)
# probability_genotype_fast_greta <- function(p, L, R) {
#   prob_left  <- sweep(L, 2, 1 - p, "*") + sweep(1 - L, 2, p, "*")
#   prob_right <- sweep(R, 2, 1 - p, "*") + sweep(1 - R, 2, p, "*")
#   dup <- 1 + L - R
#   F <- prob_left * prob_right * dup
#   z <- apply(F, 1, prod)
#   z / sum(z)
# }
# greta version
probability_genotype_fast_greta <- function(p, L, R) {
  prob_left  <- sweep(L, 2, 1 - p, "*") + sweep(1 - L, 2, p, "*")
  prob_right <- sweep(R, 2, 1 - p, "*") + sweep(1 - R, 2, p, "*")
  dup <- 1 + L - R
  # z <- apply(prob_left * prob_right * dup, 1, "prod")
  z <- exp(log(prob_left * prob_right * dup + 1e-12) %*% rep(1, ncol(L)))
  # z <- exp(rowSums(log(prob_left * prob_right * dup + 1e-12)))
  z / sum(z)
}


# Multilocus polygenic selection step
# polygenic_multilocus_next_step <- function(z, w, h, L, R) {
#   G <- nrow(L); B <- ncol(L)
#   Gw <- matrix(NA_real_, nrow = G, ncol = B)
#   for (lo in seq_len(B)) {
#     SS <- (L[, lo] == 1L) & (R[, lo] == 1L)
#     RR <- (L[, lo] == 0L) & (R[, lo] == 0L)
#     SR <- (L[, lo] != R[, lo])
#     Gw[, lo] <- 1 * SS + w[lo] * RR + (h[lo] * w[lo] + (1 - h[lo])) * SR
#   }
#   r_vec <- apply(Gw, 1, prod)
#   genotype_post <- z * r_vec
#   genotype_post / sum(genotype_post)
# }
# version Greta
polygenic_multilocus_next_step_greta <- function(z, w, h, SS_mask, RR_mask, SR_mask) {
  Gw <- SS_mask +
    sweep(RR_mask, 2, w, "*") +
    sweep(SR_mask, 2, h * w + (1 - h), "*")
  r_vec <- exp(log(Gw + 1e-12) %*% rep(1, ncol(Gw)))
  # r_vec <- exp(rowSums(log(Gw + 1e-12)))
  # r_vec <- apply(Gw, 1, "prod")
  genotype_post <- z * r_vec
  genotype_post / sum(genotype_post)
}


# Dirichlet-Multinomial sampling
sample_genotype_counts <- function(Z_true, M_z, rho_z) {
  alpha0_z <- (1 - rho_z^2)/rho_z^2
  alpha <- Z_true * alpha0_z
  Z_disp <- as.numeric(rdirichlet(1, alpha))
  as.vector(rmultinom(1, M_z, Z_disp))
}

# gneotype frequency through time
# simulate_genotype_timecourse <- function(L, R, p, w, h, Tmax, M_z, rho_z) {
#   Z_list <- list()
#   Z_list[[1]] <- probability_genotype_fast(p = p, L = L, R = R)
#   
#   if (Tmax > 1) {
#     for (t in 2:Tmax) {
#       Z_true <- polygenic_multilocus_next_step(
#         z = Z_list[[t-1]],
#         w = w,
#         h = h,
#         L = L,
#         R = R
#       )
#       Z_list[[t]] <- Z_true
#     }
#   }
#   
#   N_list <- list()
#   for (t in seq_len(Tmax)) {
#     N_list[[t]] <- sample_genotype_counts(Z_true = Z_list[[t]], M_z = M_z, rho_z = rho_z)
#   }
#   
#   Z_mat <- do.call(cbind, Z_list)
#   N_mat <- do.call(cbind, N_list)
#   
#   list(Z_list = Z_list, N_list = N_list, Z_mat = Z_mat, N_mat = N_mat)
# }

# version greta
simulate_genotype_timecourse_greta <- function(p, w, h, Tmax,
                                               SS_mask, RR_mask, SR_mask,
                                               L, R) {
  Z_list <- vector("list", Tmax)
  Z_list[[1]] <- probability_genotype_fast_greta(p, L, R)
  if (Tmax > 1) {
    for (t in 2:Tmax) {
      Z_list[[t]] <- polygenic_multilocus_next_step_greta(
        z = Z_list[[t - 1]], w = w, h = h,
        SS_mask = SS_mask, RR_mask = RR_mask, SR_mask = SR_mask
      )
    }
  }
  Z_list
}

# Phenotype function
# compute_Ugc <- function(L, R, w, h) {
#   G <- nrow(L)
#   B <- ncol(L)
#   f <- matrix(NA, nrow = G, ncol = B)
#   
#   for (l in 1:B) {
#     SS <- (L[, l] == 1 & R[, l] == 1)
#     RR <- (L[, l] == 0 & R[, l] == 0)
#     SR <- (L[, l] == 1 & R[, l] == 0) | (L[, l] == 0 & R[, l] == 1)
#     
#     f[, l] <- 1 * SS +
#       w[l] * RR +
#       (h[l] * w[l] + (1 - h[l])) * SR
#   }
#   
#   return(f)
# }
### PHENOTYPE FUNCTION
# compute_Ugc <- function(theta, h, SS_mask, RR_mask, SR_mask) {
#   SS_mask +
#     sweep(RR_mask, 2, theta, "*") +
#     sweep(SR_mask, 2, h * theta + (1 - h), "*")
# }

compute_Ugc <- function(theta, h, SS_mask, RR_mask, SR_mask) {
  sweep(RR_mask, 2, theta, "*") +
    sweep(SR_mask, 2, h * theta, "*")
  # SS_mask contributes nothing — omitted entirely, since susceptible = 0 hazard
}

# compute_Ustar <- function(Ugc, theta, type = "", epsilon = NULL) {
#   
#   U_add <- rowSums(Ugc)
#   U_mult <- apply(Ugc, 1, prod)
#   
#   if (type == "additive") {
#     U_star <- U_add
#     
#   } else if (type == "multiplicative") {
#     U_star <- U_mult
#     
#   } else if (type == "epistatic") {
#     if (is.null(epsilon)) stop("epsilon required")
#     pairwise_epi <- rowSums(Ugc %*% epsilon * Ugc)
#     U_star <- U_add + U_mult + pairwise_epi
#   }
#   
#   U_star <- U_star / sum(theta)
#   return(U_star)
# }

# the function following is allowing us to estimate the pehenotype frequency
# from the genotype

# multiplicative (indepedent barriers) is goverm by the assumptions:

## - that genotype of differents locus are acting as "filter survival"
## that means that if we have two locus: one metabolic and one target site
## the insecticide will be metabolised by the enzyme metabolic then pass
## and will arrive in the target site gene but this gene already modified
##  the proteins
# compute_Ustar_multiplicative<- function(Ugc, theta, epsilon = NULL) {
#   
#   # U_mult <- apply(Ugc, 1, "prod")
#   # U_mult <- exp(rowSums(log(Ugc + 1e-12)))
#   U_mult <- exp(log(Ugc + 1e-12) %*% rep(1, ncol(Ugc)))
#   U_star <- U_mult # normailizing this with sum(theta was a bad idea)
#   return(U_star)
# }

# compute_Ustar_multiplicative <- function(Ugc, theta = NULL, epsilon = NULL) {
#   U_mult <- exp(log(Ugc + 1e-12) %*% rep(1, ncol(Ugc)))
#   return(U_mult)
# }
# 
# compute_p_died_multiplicative <- function(U_star) {
#   eps <- 1e-6
#   p <- eps + (1 - 2 * eps) * U_star
#   return(p)
# }
# additive (contribution for the protection)
# assumptions is govern by the assumptions:

## - that the genotype of differents locus will contribute a small on the resistant
## so here it's not like checking anymore but each resistance will
## add some protection to the mosquitoes so acting on the survival
compute_Ustar_additive <- function(Ugc, theta, epsilon = NULL) {
  # U_add <- apply(Ugc, 1, "sum")
  U_add <- Ugc %*% rep(1, ncol(Ugc))
  # U_add <- rowSums(Ugc)
  U_star <- U_add # here I some it by the total of Ugc instead of theta
  # this should be discussed with Nick
  return(U_star)
}

## computing U_star
compute_p_died_additive <- function(U_star) {
  p <- 1 - exp(-U_star)
  return(p)
}


# compute_p_died <- function(U_star) {
#   p <- 1 - U_star
#   p[p < 0] <- 0
#   p[p > 1] <- 1
#   return(p)
# } # the version is hated by greta 

# compute_p_died <- function(U_star) {
#   p <- U_star # 1- ustar give us the survival of resistant
#   p <- p * (p > 0)              # zero out anything below 0
#   p <- p * (p <= 1) + (p > 1)   # cap anything above 1 at exactly 1
#   return(p)
# }
# compute_p_died_multiplicative <- function(U_star) {
#   eps <- 1e-6
#   p <- eps + (1 - 2*eps) * U_star   # linearly rescales [0,1] into (eps, 1-eps)
#   return(p)
# }
# this is to model the phenotype based on his likelihood here

simulate_beta_binomial <- function(p, n, phi = 20) {
  alpha <- p * phi
  beta  <- (1 - p) * phi
  p_sample <- rbeta(length(p), alpha, beta)
  y <- rbinom(length(p), size = n, prob = p_sample)
  return(list(p_sample = p_sample, y = y))
}


betabinomial_p_rho <- function(N, p, rho) {
  
  # model the observation (betabinomial) sd as a multiplier on the binomial sd,
  # accounting for additional error due to nonindependent sampling of individuals
  # from the population. This is based on the INLA parameterisation
  
  # solve for a and b:
  #   p = a / (a + b)
  #   rho = 1 / (a + b + 1)
  a <- p * (1 / rho - 1)
  b <- a * (1 - p) / p
  
  # define betabinomial according to the greta interface
  beta_binomial(size = N, alpha = a, beta = b)
  
}

nearish <- function(x, y) {
  dplyr::near(x, y, tol = 1e-2)
}
# pheno <- function(Z_mat, effect_type = "additive", alpha = NULL, epsilon = NULL) {
#   p <- ncol(Z_mat)
#   
#   if (is.null(alpha)) alpha <- runif(p, 0.1, 1)
#   
#   if (effect_type == "additive") {
#     pheno <- as.vector(as.matrix(Z_mat) %*% alpha)
#     
#   } else if (effect_type == "epistatic") {
#     if (is.null(epsilon)) epsilon <- matrix(runif(p^2, -0.2, 0.2), p, p)
#     pheno <- as.vector(as.matrix(Z_mat) %*% alpha +
#                          rowSums((Z_mat %*% epsilon) * Z_mat))
#   }
#   
#   pheno <- (pheno - min(pheno)) / (max(pheno) - min(pheno))
#   return(pheno)
# }

getwd()

# Allele frequency from genotype

# allele_frequency_next_step <- function(genotype_next, L, R) {
#   G <- nrow(L); B <- ncol(L)
#   p_next <- numeric(B)
#   for (lo in seq_len(B)) {
#     a_lo <- 2L - L[, lo] - R[, lo]  # SS=0, SR=1, RR=2
#     p_next[lo] <- 0.5 * sum(genotype_next * a_lo)
#   }
#   p_next
# }

allele_frequency_from_genotype_greta <- function(Z, L, R) {
  a_mat <- 2 - L - R   # G x n_loci: SS=0, SR=1, RR=2 resistant-allele count
  p_alleles <- 0.5 * (t(a_mat) %*% Z)   # (n_loci x G) %*% (G x 1) = n_loci x 1
  return(p_alleles)
}
