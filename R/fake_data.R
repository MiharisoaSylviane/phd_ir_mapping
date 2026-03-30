library(tidyverse)
library(MCMCpack)

set.seed(123)

# Dummy matrices 
source("C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/R/create_dummy_matrices.R")

n_loci <- 2
mats   <- create_dummy_matrices(n_loci)
L <- mats$left
R <- mats$right
G <- nrow(L)
B <- ncol(L)

# 3 loci : 2 kdr mutations L1014F and L1014S (WT = L1014L) + Ace1
locus_names <- c("L1014F", "L1014S")
colnames(L) <- locus_names # name of the locus will be the SNP
colnames(R) <- locus_names 

# loading Model functions 
probability_genotype_fast <- function(p, L, R) {
  prob_left  <- sweep(L, 2, 1 - p, "*") + sweep(1 - L, 2, p, "*")
  prob_right <- sweep(R, 2, 1 - p, "*") + sweep(1 - R, 2, p, "*")
  dup <- 1 + L - R
  z   <- apply(prob_left * prob_right * dup, 1, prod)
  z / sum(z)
}

polygenic_multilocus_next_step <- function(z, w, h, L, R) {
  G  <- nrow(L); B <- ncol(L)
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

sample_genotype_counts <- function(Z_true, M_z, rho_z) {
  alpha0_z <- (1 - rho_z^2) / rho_z^2
  alpha    <- pmax(Z_true * alpha0_z, 1e-8)
  Z_disp   <- as.numeric(rdirichlet(1, alpha))
  as.vector(rmultinom(1, M_z, Z_disp / sum(Z_disp)))
}

# Z_list = true genotype frequencies (biological reality) stored as list
# Z_mat  = same thing as a G x Tmax matrix
# N_list = observed counts (field reality) stored as list
# N_mat  = same thing as a G x Tmax matrix
simulate_genotype_timecourse <- function(L, R, p, w, h, Tmax, M_z, rho_z) {
  Z_list      <- vector("list", Tmax)
  Z_list[[1]] <- probability_genotype_fast(p = p, L = L, R = R)
  if (Tmax > 1) {
    for (t in 2:Tmax) {
      Z_list[[t]] <- polygenic_multilocus_next_step(
        z = Z_list[[t - 1]], w = w, h = h, L = L, R = R
      )
    }
  }
  N_list <- lapply(seq_len(Tmax), function(t)
    sample_genotype_counts(Z_true = Z_list[[t]], M_z = M_z, rho_z = rho_z)
  )
  list(
    Z_list = Z_list, N_list = N_list,
    Z_mat  = do.call(cbind, Z_list), #do.call : take each elements in z_list
    # we have a list Zlist
    N_mat  = do.call(cbind, N_list) # cbind : put column side by side
    # so saying take each element of the list with do.call and put them in colqumn side by side
  )
}

# Genotype label lookup 
genotype_label_df <- map_dfc(locus_names, function(loc) {
  tibble(!!loc := case_when(
    L[, loc] == 1 & R[, loc] == 1 ~ "SS",
    L[, loc] == 0 & R[, loc] == 0 ~ "RR",
    TRUE                           ~ "SR"
  ))
})

#  True parameters 
w_true     <- c(1.3, 1.2, 1.1)   # fitness of RR relative to SS
h_true     <- c(0.5, 0.3, 0.7)   # dominance
rho_z_true <- 0.05                # overdispersion
Tmax       <- 8                   # number of generations
M_z        <- 100                 # mosquitoes sampled per timepoint

# Village metadata 
villages <- tibble(
  village   = c("Nkolondom", "Nkolbikon", "Campo",
                "Yaounde",   "Douala",    "Bafoussam"),
  latitude  = c( 3.948,  5.602,  2.375,  3.848,  4.050,  5.478),
  longitude = c(11.505, 13.675,  9.826, 11.502,  9.700, 10.421)
)

#  Genotype lookup table 
# TIBBLE 1 : all 27 possible multilocus genotype combinations
genotype_lookup <- genotype_label_df %>%
  mutate(genotype_id = row_number(), .before = 1)


# Simulate one dataset per village 
# Helper to build mcmc_data for one village
build_village_data <- function(village_name, latitude, longitude,
                               L, R, w, h, Tmax, M_z, rho_z,
                               genotype_lookup, locus_names) {
  
  # each village gets its own random starting allele frequencies (and shouldnt be static)
  p_village <- runif(ncol(L), 0.1, 0.9)
  
  sim <- simulate_genotype_timecourse(
    L     = L,
    R     = R,
    p     = p_village,
    w     = w,
    h     = h,
    Tmax  = Tmax,
    M_z   = M_z,
    rho_z = rho_z
  )
  
  # observed counts from N_mat
  obs <- genotype_lookup %>%
    bind_cols(# add dataframe in top of each other
      as.data.frame(sim$N_mat) %>% 
        setNames(paste0("t", seq_len(Tmax))) #setnames: add the t on each sequence of the timepoint
    ) %>%
    pivot_longer(# tidy the table to be more long
      cols      = starts_with("t"),
      names_to  = "timepoint",
      values_to = "n_observed"
    ) %>%
    mutate(timepoint = as.integer(str_remove(timepoint, "t"))) # 
  
  # true frequencies from Z_mat
  true_freq <- genotype_lookup %>%
    bind_cols(
      as.data.frame(sim$Z_mat) %>%
        setNames(paste0("t", seq_len(Tmax)))
    ) %>%
    pivot_longer(
      cols      = starts_with("t"),
      names_to  = "timepoint",
      values_to = "z_true"
    ) %>%
    mutate(timepoint = as.integer(str_remove(timepoint, "t")))
  
  # n_tested per timepoint
  n_tested_df <- tibble(
    timepoint = seq_len(Tmax),
    n_tested  = as.integer(colSums(sim$N_mat))
  )
  
  # join everything together
  obs %>%
    left_join(n_tested_df, by = "timepoint") %>% #merge two dataframe
    left_join(true_freq,
              by = c("genotype_id", all_of(locus_names), "timepoint")) %>%
    mutate(
      freq_observed = round(n_observed / n_tested, 3),
      z_true        = round(z_true,                3),
      village       = village_name,
      latitude      = latitude,
      longitude     = longitude,
      p_init_L1014F = p_village[1],
      p_init_L1014S = p_village[2],
      p_init_Ace1   = p_village[3]
    ) %>%
    dplyr::select(
      village, latitude, longitude,
      p_init_L1014F, p_init_L1014S, p_init_Ace1,
      genotype_id,
      all_of(locus_names),
      timepoint,
      n_observed,
      n_tested,
      freq_observed,
      z_true
    ) %>%
    arrange(timepoint, genotype_id)
}

#  Build mcmc_data for all villages
# TIBBLE 2 : master data tibble — one row per village x genotype x timepoint
# 6 villages x 27 genotypes x 8 timepoints = 1296 rows

mcmc_data_all <- map_dfr(
  seq_len(nrow(villages)),
  function(i) {
    build_village_data(
      village_name   = villages$village[i],
      latitude       = villages$latitude[i],
      longitude      = villages$longitude[i],
      L              = L,
      R              = R,
      w              = w_true,
      h              = h_true,
      Tmax           = Tmax,
      M_z            = M_z,
      rho_z          = rho_z_true,
      genotype_lookup = genotype_lookup,
      locus_names    = locus_names
    )
  }
)

# to check if the data have latitude longitude
mcmc_data_all %>%
  count(village, latitude, longitude) %>%
  print()

# Build freq_marginal for all villages
# Collapse across other loci — marginal frequency per locus x genotype x time

freq_marginal <- bind_rows(# bind_rows: puts data frames on top of each other
  # needs same columns
  # bind_cols : puts data frames side by side
  # needs same number of rows
  mcmc_data_all %>%
    group_by(village, latitude, longitude,
             timepoint, genotype_class = L1014F) %>%
    summarise(n        = sum(n_observed),
              n_tested = first(n_tested),
              z_true   = sum(z_true),
              .groups  = "drop") %>%
    mutate(locus = "L1014F"),
  
  mcmc_data_all %>%
    group_by(village, latitude, longitude,
             timepoint, genotype_class = L1014S) %>%
    summarise(n        = sum(n_observed),
              n_tested = first(n_tested),
              z_true   = sum(z_true),
              .groups  = "drop") %>%
    mutate(locus = "L1014S"),
  
  mcmc_data_all %>%
    group_by(village, latitude, longitude,
             timepoint, genotype_class = Ace1) %>%
    summarise(n        = sum(n_observed),
              n_tested = first(n_tested),
              z_true   = sum(z_true),
              .groups  = "drop") %>%
    mutate(locus = "Ace1")
) %>%
  mutate(
    freq_observed = round(n      / n_tested, 3), # round get decimal
    z_true        = round(z_true,            3),
    col_obs       = paste0(genotype_class, "_obs_",  locus),# paste0 : combine text
    col_true      = paste0(genotype_class, "_true_", locus), # paste0 : combine text
    col_n         = paste0(genotype_class, "_n_",    locus) # paste0 : combine text
  )

# Wide format tibbles
# TIBBLE 3a : observed frequencies
# one row per village x timepoint — loci as columns

freq_wide_obs <- freq_marginal %>%
  dplyr::select(village, latitude, longitude,
         timepoint, col_name = col_obs, value = freq_observed) %>%
  pivot_wider(names_from = col_name, values_from = value) %>%
  dplyr::select(
    village, latitude, longitude, timepoint,
    SS_obs_L1014F, SR_obs_L1014F, RR_obs_L1014F,
    SS_obs_L1014S, SR_obs_L1014S, RR_obs_L1014S,
    SS_obs_Ace1,   SR_obs_Ace1,   RR_obs_Ace1
  ) %>%
  arrange(village, timepoint)%>%
  as.data.frame()

# TIBBLE 3b : true frequencies
freq_wide_true <- freq_marginal %>%
  dplyr::select(village, latitude, longitude,
         timepoint, col_name = col_true, value = z_true) %>%
  pivot_wider(names_from = col_name, values_from = value) %>%
  dplyr::select(
    village, latitude, longitude, timepoint,
    SS_true_L1014F, SR_true_L1014F, RR_true_L1014F,
    SS_true_L1014S, SR_true_L1014S, RR_true_L1014S,
    SS_true_Ace1,   SR_true_Ace1,   RR_true_Ace1
  ) %>%
  arrange(village, timepoint)

# TIBBLE 3c : raw counts
freq_wide_counts <- freq_marginal %>%
  dplyr::select(village, latitude, longitude,
         timepoint, col_name = col_n, value = n) %>%
  pivot_wider(names_from = col_name, values_from = value) %>% # transform variable to column
  dplyr::select(
    village, latitude, longitude, timepoint,
    SS_n_L1014F, SR_n_L1014F, RR_n_L1014F,
    SS_n_L1014S, SR_n_L1014S, RR_n_L1014S,
    SS_n_Ace1,   SR_n_Ace1,   RR_n_Ace1
  ) %>%
  arrange(village, timepoint)

# # Checks 
# 
# mcmc_data_all %>%
#   group_by(village, timepoint) %>%
#   summarise(
#     sum_observed = sum(n_observed),
#     n_tested     = first(n_tested),
#     ok           = sum_observed == n_tested,
#     .groups      = "drop"
#   ) %>%
#   print()  
# 
# 
# freq_wide_obs %>%
#   mutate(
#     sum_L1014F = SS_obs_L1014F + SR_obs_L1014F + RR_obs_L1014F,
#     sum_L1014S = SS_obs_L1014S + SR_obs_L1014S + RR_obs_L1014S,
#     sum_Ace1   = SS_obs_Ace1   + SR_obs_Ace1   + RR_obs_Ace1
#   ) %>%
#   dplyr::select(village, timepoint, sum_L1014F, sum_L1014S, sum_Ace1) %>%
#   print()

write.csv(freq_wide_obs, 
          "C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/dataoutput/fake.csv", 
          row.names = FALSE)

write.csv(mcmc_data_all, 
          "C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/dataoutput/mcmc_data.csv", 
          row.names = FALSE)
