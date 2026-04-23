library(tidyverse)
library(MCMCpack)
library(coda)

set.seed(123)

# Dummy matrices 
source("C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/R/create_dummy_matrices.R")

n_loci <- 3
mats   <- create_dummy_matrices(n_loci)
L <- mats$left
R <- mats$right
G <- nrow(L)
B <- ncol(L) # B = number of loci = 3

locus_names <- c("L1014F", "L1014S", "Ace1")
colnames(L) <- locus_names
colnames(R) <- locus_names

# Model functions 
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

# Log-likelihood
# Log-likelihood is allowing to assess how our parameters explain the data
# the lgamma function allow us to get the log-likelihood
# here as we have dirichlet-multinomial, so the probability of observing our data
# with parameter x for example, if the prediction is likely equal to the observation
# the contribution is positive and the model is good,
# if not we will have a negative contribution > model bad

log_lik_dm <- function(N_mat, Z_mat_pred, rho_z) {
  Tmax     <- ncol(N_mat)              # assign the generation number
  alpha0_z <- (1 - rho_z^2) / rho_z^2 # define the dirichlet parameter
  ll       <- 0                        # initialize the log of vraisemblance on 0
  
  # this function is saying for each t, calculate the contribution of each timepoint
  # and add it to ll
  for (t in seq_len(Tmax)) {
    
    z_t   <- pmax(Z_mat_pred[, t], 1e-10) # extract the column t, and replace the z_mat_t=0
    # to small number so we could keep the calculation
    
    alpha <- z_t * alpha0_z  # for each genotype recombined, the parameter
    # is proportional to the frequency z_t
    # so alpha is the weight of the model for each genotype
    
    n_t   <- N_mat[, t]  # count of mosquitoes observed with time
    # that have the genotypes recombined
    N_t   <- sum(n_t)    # number of mosquitoes captured
    
    ll    <- ll +
      lgamma(N_t + 1) +       # this is allowing us to say the way we arrange the locus
      # in the data doesn't matter, mathematically, this is
      # represented by the factorial function
      # like saying with 100 mosquitoes = we have
      # 100! = 100 manners to order the observation
      lgamma(alpha0_z) -
      lgamma(N_t + alpha0_z) + # because alpha0_z is from rho_z which define
      # the variability of the predicted data,
      # this is to prevent a bit randomness
      sum(lgamma(n_t + alpha) -
            lgamma(alpha) -
            lgamma(n_t + 1))  # this will allow us to compare model predicted and the data
    # if I observe more genotype than the model this will get greater
    # if I observe few genotype than the model predict,
    # this term will get smaller
    # if we have a bad model the sum here will get a big negative value
  }
  ll
}

# Parameter transformations 
# qlogis: remove the definition domain (0,1) so they will belong to R
# log is used here as W >1, and log gives us a strictly positive number

pack_params <- function(p, w, h, rho_z) {
  c(qlogis(p),       # p : (0,1)  → (-∞,+∞)
    log(w - 1),      # w : (>1)   → (-∞,+∞)
    qlogis(h),       # h : (0,1)  → (-∞,+∞)
    qlogis(rho_z))   # rho_z : (0,1) → (-∞,+∞)
}

# this is to get back to the real world values, the biological scale
# plogis will help us to get back the number as between 0 and 1
# we could use for the mcmc, so MCMC won't give us impossible numbers
# So unpack_params are the function to get the table of the parameters
unpack_params <- function(theta, B) {
  list(
    # this is to make the table of output in depend of the loci
    p     = plogis(theta[1:B]),                    # example 3 loci, we will have 3 values of p
    w     = 1 + exp(theta[(B + 1):(2 * B)]),       # B+1 because at 3 loci, the 3rd row will
    # be for p, so next row for the w [4:6],...
    h     = plogis(theta[(2 * B + 1):(3 * B)]),    # rows [7:9]
    rho_z = plogis(theta[3 * B + 1])               # row [10] — one shared value
  )
}

# Log posterior 
# We define here the log posterior
# log posterior = log likelihood + log priors
# from Bayes theorem : posterior = likelihood × priors
# and it is a lot of calculation so we are using the log function to transform it
# log(posterior) = log(likelihood) + log(priors)

make_log_posterior <- function(N_mat, L, R, Tmax, M_z) {
  function(theta) {
    
    # parameters
    B      <- ncol(L)
    params <- unpack_params(theta, B) # getting back to the parameters real world values
    
    # list of parameters we will extract
    p      <- params$p
    w      <- params$w
    h      <- params$h
    rho_z  <- params$rho_z
    
    if (any(!is.finite(c(p, w, h, rho_z)))) return(-Inf) # checking if there is NA, NaN or Inf
    if (rho_z <= 0 | rho_z >= 1)            return(-Inf) # checking if rho_z is between 0 and 1
    if (any(w <= 1))                         return(-Inf) # checking if w>1
    # if not the resistant won't spread
    
    # Z_list has Tmax vectors — one per generation — each with 27 genotype frequencies
    # number of vectors depends on Tmax — if Tmax=8 we get 8 vectors → matrix 27 × 8
    Z_list      <- vector("list", Tmax)
    
    # with parameters we simulate what the model predicts
    Z_list[[1]] <- probability_genotype_fast(p, L, R)  # Hardy-Weinberg at t=1
    for (t in seq_len(Tmax - 1)) {
      Z_list[[t + 1]] <- polygenic_multilocus_next_step(
        z = Z_list[[t]], w = w, h = h, L = L, R = R   # selection step
      )
    }
    
    # Z_mat_pred is the frequency given by the model
    # do.call(cbind) assembles the Tmax vectors side by side — combine genotype × timepoint
    Z_mat_pred <- do.call(cbind, Z_list)
    
    # log of likelihood of the data
    # measures how well our prediction corresponds to the data
    ll <- log_lik_dm(N_mat, Z_mat_pred, rho_z)
    if (!is.finite(ll)) return(-Inf)
    
    # priors — assess how plausible the parameters are biologically
    lp_p   <- sum(dbeta(p,     2,  2,  log = TRUE)) # density of plausible value p ~ Beta(2,2)
    lp_w   <- sum(dexp(w - 1,  1,      log = TRUE)) # w-1 ~ Exp(1)
    lp_h   <- sum(dbeta(h,     1,  1,  log = TRUE)) # h ~ Beta(1,1), no preference for dominant or recessive
    lp_rho <- dbeta(rho_z,     2,  20, log = TRUE)  # rho_z ~ Beta(2,20)
    
    # log posterior = log likelihood + log priors
    ll + lp_p + lp_w + lp_h + lp_rho
  }
}

#  Metropolis-Hastings sampler
# Manual MCMC sampler — no dependency on Hessian — always works
# Replaces MCMCmetrop1R to avoid the Hessian error

metropolis_sampler <- function(log_post, theta_init,
                               n_iter   = 5000,
                               n_burnin = 1000,
                               tune     = 0.15) {
  
  n_params <- length(theta_init)  # number of parameters = 10 with 3 loci
  chain    <- matrix(NA_real_, nrow = n_iter, ncol = n_params)  # empty matrix to store chain
  theta    <- theta_init          # start at the initial point
  lp_curr  <- log_post(theta)    # log posterior at start — used as reference
  n_accept <- 0                   # acceptance counter
  
  for (i in seq_len(n_iter)) {
    
    # propose a new point — Gaussian random walk
    # add random noise to current position to get candidate
    theta_prop <- theta + rnorm(n_params, 0, tune)
    
    # evaluate log posterior at the proposal
    # this runs the full biological simulation + likelihood + priors
    lp_prop <- log_post(theta_prop)
    
    # Metropolis-Hastings acceptance criterion
    # difference between proposal and current log posterior
    log_alpha <- lp_prop - lp_curr
    
    # accept if proposal is better (log_alpha > 0 → always accept)
    # or accept with probability exp(log_alpha) if proposal is worse
    # this prevents the chain from getting stuck in a local maximum
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      theta    <- theta_prop    # move to proposal
      lp_curr  <- lp_prop       # update current log posterior
      n_accept <- n_accept + 1  # increment acceptance counter
    }
    
    chain[i, ] <- theta  # record current position whether accepted or not
  }
  
  # remove the burn-in — first n_burnin steps are not representative of posterior
  # the chain is still searching for the high probability region during burn-in
  chain_post  <- chain[(n_burnin + 1):n_iter, ]
  accept_rate <- n_accept / n_iter
  
  cat(sprintf("  Taux d'acceptation : %.3f", accept_rate))
  if (accept_rate < 0.10) {
    cat(" ← trop bas — réduire tune\n")   # tune too large — proposals jump too far
  } else if (accept_rate > 0.60) {
    cat(" ← trop élevé — augmenter tune\n") # tune too small — proposals stay too close
  } else {
    cat(" ← OK\n")
  }
  
  list(
    chain       = chain_post,   # 4000 posterior samples — matrix 4000 × 10
    accept_rate = accept_rate
  )
}

# Read field data
# read the data
# fake.csv = field data in IR Mapper format
# SS_obs, SR_obs, RR_obs columns are PROPORTIONS (0 to 1) — not raw counts

field_data <- read.csv(
  "C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/dataoutput/fake.csv"
)

cat("=== Données chargées ===\n")
cat("Villages   :", paste(unique(field_data$village), collapse = ", "), "\n")
cat("Timepoints :", sort(unique(field_data$timepoint)), "\n\n")

# Read true parameters from mcmc_data_all.csv
# These are the real parameters used to generate fake.csv
# In a real study these would be UNKNOWN
# Here we use them only to evaluate how well the MCMC recovers them

mcmc_data_all <- read.csv(
  "C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/dataoutput/mcmc_data.csv"
)

# true biological parameters — same for all villages
w_true     <- c(1.3, 1.2, 1.1)  # true fitness per locus
h_true     <- c(0.5, 0.3, 0.7)  # true dominance per locus
rho_z_true <- 0.05               # true overdispersion

# true p_init per village — stored in mcmc_data_all when we generated the data
p_init_true <- mcmc_data_all %>%
  group_by(village) %>%
  slice(1) %>%
  dplyr::select(village,
                p_init_L1014F,
                p_init_L1014S,
                p_init_Ace1) %>%
  ungroup()

p_init_true

#  Genotype lookup table 
# To save the possible genotype in a table
# Maps each genotype_id (1 to 27) to its SS/SR/RR labels per locus
# Used in reconstruct_N_mat to know which genotype corresponds to which combination

genotype_lookup <- map_dfc(locus_names, function(loc) {
  tibble(!!loc := case_when(
    L[, loc] == 1 & R[, loc] == 1 ~ "SS",
    L[, loc] == 0 & R[, loc] == 0 ~ "RR",
    TRUE                           ~ "SR"
  ))
}) %>%
  mutate(genotype_id = row_number(), .before = 1)

# Estimate p_init
# Estimer p_init
# In a real study the true starting allele frequencies are unknown
# We estimate them from observed genotype frequencies at the first timepoint
# resistant allele frequency = (2 × freq_RR + freq_SR) / 2
# p_susceptible = 1 - resistant allele frequency
# FIX: fake.csv stores PROPORTIONS not counts — do NOT divide by n_total

estimate_p_init <- function(field_data, village_name) {
  first_t <- field_data %>%
    filter(village == village_name,
           timepoint == min(timepoint))
  
  freq_R_L1014F <- (2 * first_t$RR_obs_L1014F + first_t$SR_obs_L1014F) / 2
  freq_R_L1014S <- (2 * first_t$RR_obs_L1014S + first_t$SR_obs_L1014S) / 2
  freq_R_Ace1   <- (2 * first_t$RR_obs_Ace1   + first_t$SR_obs_Ace1)   / 2
  
  cat("  Fréq. allèle R :",
      round(c(freq_R_L1014F, freq_R_L1014S, freq_R_Ace1), 3), "\n")
  
  # p = susceptible allele frequency = 1 - resistant allele frequency
  p <- 1 - c(freq_R_L1014F, freq_R_L1014S, freq_R_Ace1)
  
  # keep away from boundaries 0 and 1 to avoid numerical issues
  pmin(pmax(p, 0.05), 0.95)
}

# Reconstruct N_mat 
# redefine N_mat
# fake.csv only has marginal frequencies per locus separately
# MCMC needs a 27 × Tmax count matrix
# We reconstruct joint counts assuming independence between loci:
# P(SR at L1014F AND SS at L1014S AND RR at Ace1)
#   ≈ P(SR at L1014F) × P(SS at L1014S) × P(RR at Ace1)

reconstruct_N_mat <- function(field_data, village_name,
                              genotype_lookup, M_z = 100) {
  obs_v <- field_data %>%
    filter(village == village_name) %>%
    arrange(timepoint)
  
  Tmax <- nrow(obs_v)  # number of timepoints
  G    <- nrow(genotype_lookup)  # number of genotypes = 27
  N_mat_approx <- matrix(0L, nrow = G, ncol = Tmax)
  
  for (t in seq_len(Tmax)) {
    
    # marginal observed frequencies — already proportions in fake.csv
    freq <- list(
      L1014F = c(SS = obs_v$SS_obs_L1014F[t],
                 SR = obs_v$SR_obs_L1014F[t],
                 RR = obs_v$RR_obs_L1014F[t]),
      L1014S = c(SS = obs_v$SS_obs_L1014S[t],
                 SR = obs_v$SR_obs_L1014S[t],
                 RR = obs_v$RR_obs_L1014S[t]),
      Ace1   = c(SS = obs_v$SS_obs_Ace1[t],
                 SR = obs_v$SR_obs_Ace1[t],
                 RR = obs_v$RR_obs_Ace1[t])
    )
    
    for (loc in names(freq)) {
      s <- sum(freq[[loc]])  # s is the sum of the genotype frequencies at the locus
      if (abs(s - 1) > 0.01) freq[[loc]] <- freq[[loc]] / s  # renormalise if not exactly 1
    }
    
    # for each of the 27 multilocus genotypes
    # multiply marginal frequencies to get approximate joint frequency
    for (g in seq_len(G)) {
      joint_freq <- freq$L1014F[ genotype_lookup$L1014F[g] ] *
        freq$L1014S[ genotype_lookup$L1014S[g] ] *
        freq$Ace1[   genotype_lookup$Ace1[g]   ]
      N_mat_approx[g, t] <- round(joint_freq * M_z)
    }  # to check for each genotype
    
    # fix rounding so column sums exactly to M_z
    diff <- M_z - sum(N_mat_approx[, t])
    if (diff != 0) {
      idx <- which.max(N_mat_approx[, t])
      N_mat_approx[idx, t] <- N_mat_approx[idx, t] + diff
    }
  }
  N_mat_approx
}

# Find best starting point 
# Find the good initial condition
# Test n_attempts random starting points and keep the one with highest log posterior
# A good starting point prevents long burn-in and convergence failures

find_best_start <- function(log_post, p_init_est, B, n_attempts = 10) {
  
  best_theta <- NULL
  best_lp    <- -Inf
  
  for (i in seq_len(n_attempts)) {
    theta_cand <- pack_params(
      # combine all parameters into one vector with pack_params
      p     = pmin(pmax(p_init_est + rnorm(B, 0, 0.05), 0.05), 0.95),
      w     = runif(B, 1.05, 1.6),   # random w between 1.05 and 1.6
      h     = runif(B, 0.1,  0.9),   # random h between 0.1 and 0.9
      rho_z = runif(1, 0.02, 0.15)   # random rho_z between 0.02 and 0.15
    )
    
    lp <- log_post(theta_cand)
    
    # keep the candidate with the highest log posterior
    if (is.finite(lp) && lp > best_lp) {
      best_lp    <- lp
      best_theta <- theta_cand
    }
  }
  
  list(theta = best_theta, lp = best_lp)
}

# MCMC settings 
# Parameters
n_iter       <- 5000   # total MCMC iterations per village
n_burnin     <- 1000   # number of iterations to discard — chain still searching
tune         <- 0.15   # proposal step size — target acceptance rate 10-60%
Tmax         <- max(field_data$timepoint)
M_z          <- 100    # mosquitoes per timepoint
village_list <- unique(field_data$village)

# parameter names — 3 p + 3 w + 3 h + 1 rho_z = 10 parameters total
param_names <- c(
  paste0("p_", locus_names),   # p_L1014F, p_L1014S, p_Ace1
  paste0("w_", locus_names),   # w_L1014F, w_L1014S, w_Ace1
  paste0("h_", locus_names),   # h_L1014F, h_L1014S, h_Ace1
  "rho_z"
)

results <- list()

# Main loop — one MCMC per village 
# Boucle principale
for (vill in village_list) {
  
  cat(sprintf("\n=== Village : %s ===\n", vill))
  
  # N_mat
  # reconstruct N_mat from marginal proportions in fake.csv
  N_mat_obs <- reconstruct_N_mat(field_data, vill, genotype_lookup, M_z)
  cat("  Dimensions N_mat :", nrow(N_mat_obs), "x", ncol(N_mat_obs), "\n")
  cat("  Sommes colonnes  :", colSums(N_mat_obs), "\n")
  
  # p_init
  # estimate p_init from first timepoint
  p_init_est <- estimate_p_init(field_data, vill)
  cat("  p_init estimé    :", round(p_init_est, 3), "\n")
  
  # log posterior
  # build log posterior locked to this village's data
  log_post <- make_log_posterior(N_mat_obs, L, R, Tmax, M_z)
  
  # verification
  # verify log posterior is finite at the starting point
  theta_test <- pack_params(p_init_est, rep(1.2, B), rep(0.5, B), 0.05)
  lp_test    <- log_post(theta_test)
  cat("  Log posterior    :", round(lp_test, 2),
      if (!is.finite(lp_test)) " ← PROBLÈME" else " ← OK", "\n")
  
  if (!is.finite(lp_test)) {
    cat("  Village ignoré — log posterior non fini\n")
    next
  }
  
  # meilleur point de départ
  # find best starting point among 10 random candidates
  start <- find_best_start(log_post, p_init_est, B, n_attempts = 10)
  
  if (is.null(start$theta)) {
    cat("  Aucun point de départ valide\n")
    next
  }
  
  cat(sprintf("  Meilleur log posterior départ : %.2f\n", start$lp))
  
  # Étape 6 — MCMC manuel — pas de Hessian — fonctionne toujours
  # manual sampler — no Hessian dependency
  cat("  Lancement du MCMC...\n")
  
  mcmc_res <- metropolis_sampler(
    log_post   = log_post,
    theta_init = start$theta,
    n_iter     = n_iter,
    n_burnin   = n_burnin,
    tune       = tune
  )
  
  # transformer vers l'échelle biologique
  # transform from MCMC scale back to biological scale
  # plogis() for p and h, 1+exp() for w, plogis() for rho_z
  post_samples <- t(apply(mcmc_res$chain, 1, function(theta) {
    par <- unpack_params(theta, B)
    c(par$p, par$w, par$h, par$rho_z)
  }))
  colnames(post_samples) <- param_names
  
  # résumer le posterior
  # summarise posterior for each parameter
  post_df <- as.data.frame(post_samples)
  
  summary_village <- map_dfr(param_names, function(par) {
    samp <- post_df[[par]]
    tibble(
      village   = vill,
      parameter = par,
      post_mean = round(mean(samp),            3),
      post_med  = round(median(samp),          3),
      post_sd   = round(sd(samp),              3),
      ci_lo_95  = round(quantile(samp, 0.025), 3),
      ci_hi_95  = round(quantile(samp, 0.975), 3),
      ci_lo_50  = round(quantile(samp, 0.25),  3),
      ci_hi_50  = round(quantile(samp, 0.75),  3)
    )
  })
  
  results[[vill]] <- summary_village
  cat("  Terminé.\n")
}


# Results and comparison table 
# Output
if (length(results) == 0) {
  cat("\nAucun village n'a convergé.\n")
} else {
  
  all_results <- bind_rows(results)
  
  cat("\n=== Estimations postérieures ===\n")
  print(all_results, n = 200)
  
  # ── Comparison table : true parameters vs posterior estimates
  # This is the simulation estimation study evaluation
  # true_value = what we used to generate fake.csv — the ground truth
  # post_mean  = what the MCMC estimated from the data
  # bias       = difference between estimate and truth (0 = perfect recovery)
  # covered    = is the true value inside the 95% credible interval?
  
  comparison_long <- map_dfr(village_list, function(vill) {
    
    # true p_init for this village from mcmc_data_all
    p_true_vill <- p_init_true %>%
      filter(village == vill) %>%
      dplyr::select(p_init_L1014F, p_init_L1014S, p_init_Ace1) %>%
      unlist() %>%
      as.numeric()
    
    # full true parameter vector for this village
    # same order as param_names : p, w, h, rho_z
    true_vals <- c(
      p_true_vill,    # p_L1014F, p_L1014S, p_Ace1
      w_true,         # w_L1014F, w_L1014S, w_Ace1
      h_true,         # h_L1014F, h_L1014S, h_Ace1
      rho_z_true      # rho_z
    )
    names(true_vals) <- param_names
    
    # posterior estimates for this village
    post_vill <- all_results %>%
      filter(village == vill) %>%
      dplyr::select(parameter, post_mean, post_sd, ci_lo_95, ci_hi_95)
    
    # combine true values with posterior estimates
    post_vill %>%
      mutate(
        village    = vill,
        true_value = true_vals[parameter]  # true parameter used to generate data
      ) %>%
      dplyr::select(
        village,
        parameter,
        true_value,   # ground truth — would be unknown in a real study
        post_mean,    # MCMC estimate
        post_sd,      # uncertainty around the estimate
        ci_lo_95,
        ci_hi_95
      ) %>%
      mutate(
        # bias = difference between MCMC estimate and truth
        # 0 = perfect recovery
        bias    = round(post_mean - true_value, 3),
        # coverage = is the true value inside the 95% credible interval?
        # TRUE for most parameters = model works well
        covered = true_value >= ci_lo_95 & true_value <= ci_hi_95,
        # qualitative summary of recovery quality
        quality = case_when(
          abs(bias) <= 0.05            ~ "well recovered",
          abs(bias) >  0.05 & covered  ~ "slight bias but covered",
          abs(bias) >  0.05 & !covered ~ "poorly recovered",
          TRUE                         ~ "check"
        )
      )
  })
  
  cat("\n=== Comparison : true parameters vs posterior estimates ===\n")
  print(comparison_long, n = 200)
  
  # Summary across all villages 
  study_summary <- comparison_long %>%
    group_by(parameter, true_value) %>%
    summarise(
      n_villages   = n(),
      mean_bias    = round(mean(bias),         4),
      mean_rmse    = round(sqrt(mean(bias^2)), 4),
      coverage_95  = round(mean(covered),      3),
      mean_post_sd = round(mean(post_sd),      4),
      .groups = "drop"
    )
  
  cat("\n=== Summary across all villages ===\n")
  print(study_summary, n = 50)
  
  # true value vs posterior mean 
  # points on the diagonal = perfect recovery
  p1 <- comparison_long %>%
    ggplot(aes(x = true_value, y = post_mean)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey60") +
    geom_errorbar(aes(ymin = ci_lo_95, ymax = ci_hi_95),
                  width = 0, alpha = 0.4, colour = "#378ADD") +
    geom_point(aes(colour = covered), size = 2.5) +
    scale_colour_manual(
      values = c("TRUE"  = "#1D9E75",
                 "FALSE" = "#D85A30"),
      labels = c("TRUE"  = "CI covered truth",
                 "FALSE" = "CI missed truth"),
      name = NULL
    ) +
    facet_wrap(~ parameter, scales = "free", ncol = 4) +
    labs(
      title    = "Parameter recovery — true vs estimated",
      subtitle = "Points on diagonal = perfect recovery — bars = 95% CI",
      x        = "True value",
      y        = "Posterior mean + 95% CI"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  print(p1)
  
  # ── Graphique 2 — bias distribution 
  p2 <- comparison_long %>%
    ggplot(aes(x = bias)) +
    geom_vline(xintercept = 0,
               linetype = "dashed", colour = "grey60") +
    geom_histogram(bins = 8, fill = "#378ADD", alpha = 0.7) +
    facet_wrap(~ parameter, scales = "free_x", ncol = 4) +
    labs(
      title    = "Bias distribution across villages",
      subtitle = "Centred at 0 = unbiased estimation",
      x        = "Bias (posterior mean - true value)",
      y        = "Count"
    ) +
    theme_minimal(base_size = 11)
  print(p2)
  
  # ── Graphique 3 — 95% credible interval coverage 
  p3 <- study_summary %>%
    ggplot(aes(x = reorder(parameter, coverage_95),
               y = coverage_95)) +
    geom_hline(yintercept = 0.95,
               linetype = "dashed", colour = "grey60") +
    geom_col(fill = "#1D9E75", alpha = 0.8, width = 0.6) +
    geom_text(aes(label = scales::percent(coverage_95, accuracy = 1)),
              hjust = -0.1, size = 3) +
    coord_flip(ylim = c(0, 1.1)) +
    labs(
      title    = "95% credible interval coverage per parameter",
      subtitle = "Dashed line = target 0.95",
      x        = NULL,
      y        = "Coverage probability"
    ) +
    theme_minimal(base_size = 11)
  print(p3)
  
  # ── Graphique 4 — fitness w per locus per village
  # open circle = true value   filled point = posterior mean + 95% CI
  p4 <- comparison_long %>%
    filter(str_starts(parameter, "w_")) %>%
    ggplot(aes(x = village, y = post_mean, colour = parameter)) +
    geom_point(aes(y = true_value), shape = 1,
               size = 3, colour = "grey40") +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_lo_95, ymax = ci_hi_95),
                  width = 0.2, linewidth = 0.7) +
    geom_hline(yintercept = 1,
               linetype = "dashed", colour = "grey60") +
    facet_wrap(~ parameter, ncol = 3) +
    scale_colour_manual(values = c(w_L1014F = "#D85A30",
                                   w_L1014S = "#378ADD",
                                   w_Ace1   = "#1D9E75")) +
    labs(
      title    = "Fitness estimée par locus et par village",
      subtitle = "Moyenne + IC 95% — ligne = pas de résistance (w=1)",
      x = "Village", y = "Fitness w", colour = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  print(p4)
  
  # allele frequency p 
  p5 <- comparison_long %>%
    filter(str_starts(parameter, "p_")) %>%
    ggplot(aes(x = village, y = post_mean, colour = parameter)) +
    geom_point(aes(y = true_value), shape = 1,
               size = 3, colour = "grey40") +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_lo_95, ymax = ci_hi_95),
                  width = 0.2, linewidth = 0.7) +
    geom_hline(yintercept = 0.5,
               linetype = "dashed", colour = "grey60") +
    facet_wrap(~ parameter, ncol = 3) +
    scale_colour_manual(values = c(p_L1014F = "#D85A30",
                                   p_L1014S = "#378ADD",
                                   p_Ace1   = "#1D9E75")) +
    labs(
      title    = "Fréquence allélique susceptible initiale",
      subtitle = "Moyenne + IC 95%",
      x = "Village", y = "p", colour = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  print(p5)
  
  # Graph dominance h 
  p6 <- comparison_long %>%
    filter(str_starts(parameter, "h_")) %>%
    ggplot(aes(x = village, y = post_mean, colour = parameter)) +
    geom_point(aes(y = true_value), shape = 1,
               size = 3, colour = "grey40") +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_lo_95, ymax = ci_hi_95),
                  width = 0.2, linewidth = 0.7) +
    geom_hline(yintercept = 0.5,
               linetype = "dashed", colour = "grey60") +
    facet_wrap(~ parameter, ncol = 3) +
    scale_colour_manual(values = c(h_L1014F = "#D85A30",
                                   h_L1014S = "#378ADD",
                                   h_Ace1   = "#1D9E75")) +
    labs(
      title    = "Dominance estimée par locus et par village",
      subtitle = "0 = récessif   0.5 = additif   1 = dominant",
      x = "Village", y = "h", colour = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  print(p6)
  
  # surdispersion rho_z 
  p7 <- comparison_long %>%
    filter(parameter == "rho_z") %>%
    ggplot(aes(x = village, y = post_mean)) +
    geom_point(aes(y = true_value), shape = 1,
               size = 3, colour = "grey40") +
    geom_point(size = 3, colour = "#534AB7") +
    geom_errorbar(aes(ymin = ci_lo_95, ymax = ci_hi_95),
                  width = 0.2, linewidth = 0.7, colour = "#534AB7") +
    labs(
      title    = "Surdispersion estimée par village",
      subtitle = "Proche de 0 = peu de variabilité entre captures",
      x = "Village", y = "rho_z"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p7)
  
  # Save
  output_path <- "C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/dataoutput/"
  
  write.csv(all_results,
            paste0(output_path, "field_estimation_results.csv"),
            row.names = FALSE)
  
  write.csv(comparison_long,
            paste0(output_path, "comparison_true_vs_posterior.csv"),
            row.names = FALSE)
  
  write.csv(study_summary,
            paste0(output_path, "sim_study_summary.csv"),
            row.names = FALSE)
  
  cat("\nSauvegardé :\n")
  cat("  field_estimation_results.csv      — full posterior estimates\n")
  cat("  comparison_true_vs_posterior.csv  — true vs estimated per village\n")
  cat("  sim_study_summary.csv             — summary across villages\n")
  
  cat("\n=== Column meanings ===\n")
  cat("  true_value : parameter used to generate fake.csv — the ground truth\n")
  cat("  post_mean  : MCMC estimate — best guess from the data\n")
  cat("  post_sd    : uncertainty around the estimate\n")
  cat("  ci_lo_95   : lower bound 95% credible interval\n")
  cat("  ci_hi_95   : upper bound 95% credible interval\n")
  cat("  bias       : post_mean - true_value  (0 = perfect recovery)\n")
  cat("  covered    : TRUE if true value is inside the 95% CI\n")
  cat("  quality    : qualitative summary of parameter recovery\n")
}
 all_results

 
 # In comparison_long replace post_mean with post_med
 comparison_long <- map_dfr(village_list, function(vill) {
   
   p_true_vill <- p_init_true %>%
     filter(village == vill) %>%
     dplyr::select(p_init_L1014F, p_init_L1014S, p_init_Ace1) %>%
     unlist() %>%
     as.numeric()
   
   true_vals        <- c(p_true_vill, w_true, h_true, rho_z_true)
   names(true_vals) <- param_names
   
   post_vill <- all_results %>%
     filter(village == vill) %>%
     dplyr::select(parameter, post_med, post_sd, ci_lo_95, ci_hi_95)
   
   post_vill %>%
     mutate(
       village    = vill,
       true_value = true_vals[parameter]
     ) %>%
     dplyr::select(
       village, parameter,
       true_value,   # what we used to generate the data — the truth
       post_med,     # our best estimate from the MCMC
       post_sd,      # uncertainty
       ci_lo_95,
       ci_hi_95
     ) %>%
     mutate(
       # bias = how far our estimate is from the truth
       # 0 = perfect recovery
       bias = round(post_med - true_value, 3),
       
       # covered = is the true value inside the 95% credible interval?
       # TRUE = the model is honest about its uncertainty
       # FALSE = the model missed the truth — overconfident
       covered = true_value >= ci_lo_95 & true_value <= ci_hi_95,
       
       # qualitative summary
       quality = case_when(
         abs(bias) <= 0.05            ~ "well recovered",
         abs(bias) >  0.05 & covered  ~ "slight bias but covered",
         abs(bias) >  0.05 & !covered ~ "poorly recovered",
         TRUE                         ~ "check"
       )
     )
 })

 
 study_summary <- comparison_long %>%
   group_by(parameter, true_value) %>%
   summarise(
     n_villages    = n(),
     
     # how far is the posterior median from the truth on average
     mean_bias     = round(mean(bias),            4),
     
     # penalises large errors more — the key metric
     rmse          = round(sqrt(mean(bias^2)),    4),
     
     # proportion of villages where true value is inside 95% CI
     # target = 0.95 — if lower the model is overconfident
     coverage_95   = round(mean(covered),         3),
     
     # average uncertainty
     mean_post_sd  = round(mean(post_sd),         4),
     
     # average width of the 95% CI
     mean_ci_width = round(mean(ci_hi_95 - ci_lo_95), 4),
     
     .groups = "drop"
   )
 
 print(study_summary, n = 50)















