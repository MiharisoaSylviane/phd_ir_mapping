

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
# Compute locus-level phenotype (U_gc), perlocus level contribution

compute_Ugc <- function(L, R, w, h) {
  G <- nrow(L)
  B <- ncol(L)
  Ugc <- matrix(NA, nrow = G, ncol = B)
  for (l in 1:B) {
    SS <- (L[, l] == 1L & R[, l] == 1L)
    RR <- (L[, l] == 0L & R[, l] == 0L)
    SR <- (L[, l] != R[, l])
    Ugc[, l] <- 1 * SS + w[l] * RR + (h[l] * w[l] + (1 - h[l])) * SR
  }
  return(Ugc)
}

# Compute genotype-level phenotype (combined effects), U_star is the genotype level score

compute_Ustar_combined <- function(Ugc, theta, epsilon) {
  U_add <- rowSums(Ugc)
  U_mult <- apply(Ugc, 1, prod)
  pairwise_epi <- apply(Ugc %*% epsilon * Ugc, 1, sum)
  
  U_star <- U_add + U_mult + pairwise_epi
  U_star <- U_star / sum(theta)   # scale by theta
  U_star <- pmin(U_star, 1)       # cap at 1
  return(U_star)
}


# Compute survival over time (population and genotype level)

compute_survival_over_time <- function(Z_list, L, R, w, h, theta, epsilon) {
  Tmax <- length(Z_list)
  G <- nrow(L)
  

  Ugc <- compute_Ugc(L, R, w, h)
  pheno_genotype <- compute_Ustar_combined(Ugc, theta, epsilon)
  
  pop_survival <- numeric(Tmax)
  geno_survival <- matrix(NA, nrow = Tmax, ncol = G)
 
  # to make the survival dynamic 
  for (t in 1:Tmax) {
    freq <- Z_list[[t]]
    geno_survival[t, ] <- pheno_genotype * freq
    pop_survival[t] <- sum(pheno_genotype * freq)
  }
  
  return(list(population = pop_survival, genotype = geno_survival))
}

# Prepare genotype dataframe for plotting
prepare_geno_df <- function(pheno_res, Tmax) {
  geno_df <- as.data.frame(pheno_res$genotype[1:Tmax, ])
  colnames(geno_df) <- paste0("G", 1:ncol(geno_df))
  geno_df$Time <- 1:Tmax
  melt(geno_df, id.vars = "Time", variable.name = "Genotype", value.name = "Survival")
}


# Plot survival (population + genotypes)
plot_survival <- function(df_pop, geno_long, title) {
  ggplot() +
    geom_line(data = geno_long, aes(x = Time, y = Survival, color = Genotype), alpha = 0.5) +
    geom_line(data = df_pop, aes(x = Time, y = Survival), color = "red", size = 1.2) +
    labs(title = title, x = "Generation", y = "Survival probability") +
    theme_minimal() +
    scale_color_viridis_d(option = "plasma") +
    theme(legend.position = "right")
}
