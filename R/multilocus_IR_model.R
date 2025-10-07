library(greta)
library(MCMCpack) 
##### 
set.seed(123)
## Prior
# p is the initial allele frequency

n_loci <- 4

# call the dummy function
source("C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/R/create_dummy_matrices.R")
genotype_combination <- create_dummy_matrices(n_loci)
# defining the matrix name for the hand of genotype
L <- genotype_combination$left
L
R <- genotype_combination$right
L[3, 2]
G <- nrow(L)
B <- ncol(L)
## Multilocus genotype frequencies
# we will use the dummy matrice of the left and right hand of the allele, 
# we define a function F which is the genotype probability across the n_loci
# sweep(x, MARGIN, STATS, FUN) : it allows Sylviane to operate value in matrix,
# x= matrix name, margin = to perform on Row = 1, on column =2, stats = the value
# to use (formula), function is the * or / or +, 
# I- NO TIMESERIES
# PERLOCUS PROBABILITY THAT MATCHES THE GENOTYPES 
### 1- First way
# probability_genotype <- function(p){
  # G <- nrow(L)
  # B <- ncol(L)
  # if(length(p) != B){stop("locus probability dimension mismatch")}
  # F <- matrix(data=NA, nrow = G, ncol= B)
  # for (g in seq_len(G)) {
  #   for (lo in seq_len(B)) {
  #     prob_left <- (1-p[lo]) * L[g,lo] + p[lo]* (1-L[g,lo])
  #     prob_right <- (1-p[lo]) * R[g,lo] + p[lo]* (1-R[g,lo])
  #     heterozygote_duplicates <- 1 + L[g,lo] - R[g,lo]
  #     F[g,lo] <- prob_left * prob_right * heterozygote_duplicates
  #     
  #   }
  #   z[g] <- prod(F[g, ])
    # # or this:
    # z[g] <- 1
    # for(lo in seq_len(B)) {
    #   z[g] <- z[g] * F[g, lo]
    # }
#   }
#   z[g] <- prod(F[g, ])
# }
# p <- runif(B)

# ## 2- Second way : fast way to get genotype perlocus
# probability_genotype_fast <- function(p){
#   G <- nrow(L)
#   B <- ncol(L)
#   if(length(p) != B){stop("locus probability dimension mismatch")}
#   F <- matrix(data=NA, nrow = G, ncol= B)
# 
#   prob_left <- sweep(L, 2, 1 - p, FUN = "*") +
#     sweep(1 - L, 2, p, FUN = "*")
# 
#   prob_right <- sweep(R, 2, 1 - p, FUN = "*") + 
#     sweep(1 - R, 2, p, FUN = "*")
# 
#   heterozygote_duplicates <- 1 + L - R
# 
#   F <- prob_left * prob_right * heterozygote_duplicates
#   z <- apply(F, 1, prod)
# 
#   return(z)
# }

### 3- Third way
# p : allele frequency at each locus for each genotype (t)

p <- runif(B, 0, 1) 

probability_genotype_fast <- function(p){
  stopifnot(is.matrix(L), is.matrix(R), all(dim(L) == dim(R)))
  G <- nrow(L); B <- ncol(L)
  stopifnot(length(p) == B)
  
  prob_left  <- sweep(L, 2, 1 - p, "*") + sweep(1 - L, 2, p, "*")
  prob_right <- sweep(R, 2, 1 - p, "*") + sweep(1 - R, 2, p, "*")
  dup <- 1 + L - R   
  
  F <- prob_left * prob_right * dup
  z_unnorm <- apply(F, 1, prod)
  z <- z_unnorm / sum(z_unnorm)  
  z
}
## to test our function
z <- probability_genotype_fast(p)
z

# SELECTION PRESSURE AND RELATIVE FITNESS AT LOCUS LEVEL 
## define the covariates X with number K
K <- 4
X <- runif(K)
betamat <- matrix(rnorm(B * K), nrow= B, ncol= K)
gammax <- exp(betamat)
## selection pressure and relative fitness
s <- rep(0, B)
w <- rep(NA, B)
for(lo in seq_len(B)) {
  for (k in seq_len(K)) {
    s[lo] <- s[lo] + X[k] * gammax[lo, k]
  }
  w[lo] <- 1 + s[lo]
}

# # or
# for(lo in seq_len(B)) {
#   s[lo] <- sum(X * gammax[lo, ])
# }
# 
# # or
# s <- rowSums(sweep(gammax, 2, X, FUN = "*"))
# # or
# s <- t(gammax %*% X)
# # or
# s <- X %*% t(gammax)

# selection pressure
# st <- gammax*itns_covariates
# # relative fitness of mosquitoes with allele resistant R/S
# w <- 1 + st

#=======================================================
# GENOTYPE FREQUENCY AT NEXT STEP
## step to get the genotype frequency at next step without using time series
# polygenic_multilocus_next_step <- function(w, h){
#   stopifnot(is.matrix(L), is.matrix(R), all(dim(L) == dim(R)))
#  # define G and B 
#   G <- nrow(L); B <- ncol(L)
#   stopifnot(length(w) == B, length(h) == B)
#   stopifnot(exists("z"), length(z) == G)
#  # define an empty matrix to put the value of the fitness polygenic in depend of 
#   # the zygosity (homo - hetero) of our hand of alleles
#   Gw <- matrix(NA_real_, nrow = G, ncol = B)
#   # define the combination SS and RR and RS to give more sense to the calculation
#   for (lo in seq_len(B)) {
#     SS  <- (L[, lo] == 1L) & (R[, lo] == 1L)      # (1,1)
#     RR  <- (L[, lo] == 0L) & (R[, lo] == 0L)      # (0,0)
#     SR  <- (L[, lo] == 1L) & (R[, lo] == 0L)      # (1,0)
#     Gw[, lo] <- 1 * SS +
#       w[lo] * RR +
#       (h[lo] * w[lo] + (1 - h[lo])) * SR
#   }
#   
#   #r_vec  <- apply(Gw, 1, prod) # the 1 here is because we are operating on row
#   r_vec <- exp(rowSums(log(Gw)))
#   unpost <- z * r_vec
#   genotype_post   <- unpost / sum(unpost)
#   
#   list(Gw = Gw,
#        r = matrix(r_vec, nrow = G, ncol = 1),
#        genotype_next = matrix(genotype_post, nrow = G, ncol = 1),
#        z = matrix(z, nrow = G, ncol = 1))
# }
# # To test our function returning values
# test1 <- polygenic_multilocus_next_step(w, h)
# 
# 
# # Allele frequency at next time series p(t+1)
# allele_frequency_next_step <- function(polygenic_multilocus_next_step) {
#   stopifnot(is.matrix(L), is.matrix(R), all(dim(L) == dim(R)))
#   G <- nrow(L); B <- ncol(L)
#   stopifnot(length(polygenic_multilocus_next_step) == G)
#   
#   allele_post <- numeric(B)
#   for (lo in seq_len(B)) {
#     a_lo <- 2L - L[, lo] - R[, lo]          # length G
#     allele_post[lo] <- 0.5 * sum(genotype_next * a_lo)
#   }
#   return(allele_post)                         # length B
# }
# Draw one set of genotype proportions, here 1 represent the number of drawn we want
# we could change it to nsim (depend on how many draws we want)
# defining dominance
h <- runif(B, 0, 1)
## final function for genotype frequency at next step
polygenic_multilocus_next_step <- function(w, h){
  stopifnot(is.matrix(L), is.matrix(R), all(dim(L) == dim(R)))
  G <- nrow(L); B <- ncol(L)
  stopifnot(length(w) == B, length(h) == B)
  stopifnot(exists("z"), length(z) == G)
  Gw <- matrix(NA_real_, nrow = G, ncol = B)
  for (lo in seq_len(B)) {
    SS <- (L[, lo] == 1L) & (R[, lo] == 1L)
    RR <- (L[, lo] == 0L) & (R[, lo] == 0L)
    SR <- (L[, lo] == 1L) & (R[, lo] == 0L)
    Gw[, lo] <- 1 * SS +
      w[lo] * RR +
      (h[lo] * w[lo] + (1 - h[lo])) * SR
  }
  
  # r_vec <- exp(rowSums(log(pmax(Gw, .Machine$double.eps))))
  r_vec <- apply(Gw, 1, prod) 
  unpost <- z * r_vec
  genotype_post <- unpost / sum(unpost)
  
  list(
    Gw = Gw,
    r = matrix(r_vec, nrow = G, ncol = 1),
    genotype_next = matrix(genotype_post, nrow = G, ncol = 1), # <- this is what you use
    z = matrix(z, nrow = G, ncol = 1)
  )
}

# allele frequency at next time step
allele_frequency_next_step <- function(genotype_next) {
  stopifnot(is.matrix(L), is.matrix(R), all(dim(L) == dim(R)))
  G <- nrow(L); B <- ncol(L)
  stopifnot(length(genotype_next) == G)
  
  p_next <- numeric(B)
  for (lo in seq_len(B)) {
    a_lo <- 2L - L[, lo] - R[, lo]   # SS=0, SR=1, RR=2
    p_next[lo] <- 0.5 * sum(genotype_next * a_lo)
  }
  p_next
}
# assume z (length G) is already defined
step1 <- polygenic_multilocus_next_step(w, h)

# extract the Gx1 matrix as a vector
genotype_next <- as.vector(step1$genotype_next)   # length G

# compute allele frequencies at next step
p_next <- allele_frequency_next_step(genotype_next)
print(p_next)


# II - Loop WITH TIMESERIES
# known variables
times <- 10

# matrix storage
F_arr <- array(NA_real_, dim = c(G, B, times))  
z_mat <- matrix(NA_real_, nrow = G, ncol = times)
r_mat <- matrix(NA_real_, nrow = G, ncol = times)
Gw_arr <-array(NA_real_, dim = c(G, B, times))
p_series <- matrix(NA_real_, nrow = B, ncol = times)  
Z_post <- matrix(NA_real_, nrow = G, ncol = times)
probability_genotype_fast_t <- function(p_t){
  stopifnot(is.matrix(L), is.matrix(R), all(dim(L) == dim(R)))
  G <- nrow(L); B <- ncol(L); stopifnot(length(p_t) == B)
  
  prob_left  <- sweep(L, 2, 1 - p_t, "*") + sweep(1 - L, 2, p_t, "*")
  prob_right <- sweep(R, 2, 1 - p_t, "*") + sweep(1 - R, 2, p_t, "*")
  dup <- 1 + L - R                     
  
  F_t <- prob_left * prob_right * dup  
  z_u <- apply(F_t, 1, prod)          
  s <- sum(z_u)
  z_t <- z_u / s
  list(F = F_t, z = z_t)
}

# run over time
for (t in 1:times) {
  p_t <- runif(B, 0, 1)
  # define the allele frequency over time
  p_series[, t] <- p_t           
  # and define new matrix by adding t as indice
  out <- probability_genotype_fast_t(p_t)
  F_arr[, , t] <- out$F
  z_mat[, t]   <- out$z
  
  # selection step
  step <- polygenic_multilocus_next_step(w,h)
  Gw_arr[, , t] <- step$Gw
  r_mat[, t] <- step$r[ ,1]
  Z_post[, t] <- step$genotype_next[, 1]
}
 # check our z across tine
# p0 <- runif(B, 0, 1)
# p_series <- matrix(NA_real_, nrow = B, ncol = T)
# 
# p_t <- p0  
# for (t in 1:T) {# same prior each time
#      p_series[, t] <- p_t
#      out <- probability_genotype_fast_t(p_t)
#      F_arr[, , t] <- out$F
#      z_mat[, t]   <- out$z
# }

# to summary our results in depend of times, we will use the tranpose function t()
# to change the matrix 
# # allele frequency P (T * B)
# P_all <- t(p_series)
# colnames(P_all) <- paste0("locus_", seq_len(ncol(P_all)))
# rownames(P_all) <- paste0("t_", seq_len(nrow(P_all)))
# P_all
# # genotype frequency Z (T * G)
# Z_all <- t(z_mat)   
# colnames(Z_all) <- paste0("geno_", seq_len(ncol(Z_all)))
# rownames(Z_all) <- paste0("t_", seq_len(nrow(Z_all)))
# Z_all

# SIMULATION
# defining the matrix for our time_series in depend of locus
p_series <- matrix(NA_real_, nrow = B, ncol = times)
# giving some values for the initial allele frequencies and r
# p_series and r  should be equal to n_loci
p_series[, 1] <- p0 <- c(0.10, 0.30, 0.4, 0.02)        
r <- c(0.20, 0.35, 0.2, 0.1)                          
# definining the loop for those allele frequency at, lzed fraction of genotype frequency
# z_u, and the posterior genotype frequency z_t
for (t in 1:times) {
  # 1) allele freq at time t
  if (t == 1) {
    p_t <- p_series[, 1]
  } else {
    p_prev <- p_series[, t-1]
    p_t <- p_prev + (1 - p_prev) * r       
    p_t <- pmin(p_t, 0.999999)
  }
  p_series[, t] <- p_t
  
  # 2) build prior per-locus probabilities
  prob_left  <- sweep(L, 2, 1 - p_t, "*") + sweep(1 - L, 2, p_t, "*")
  prob_right <- sweep(R, 2, 1 - p_t, "*") + sweep(1 - R, 2, p_t, "*")
  dup <- 1 + L - R
  
  # 3) F_t and z_t
  F_t <- prob_left * prob_right * dup
  z_u <- apply(F_t, 1, prod)
  z_t <- z_u / sum(z_u)
  
  # 4) store
  F_arr[, , t] <- F_t
  z_mat[, t]   <- z_t
  
  # 5) selection step
  z <<- z_t
  step <- polygenic_multilocus_next_step(w, h)
  Gw_arr[, , t] <- step$Gw
  r_mat[, t]    <- step$r[, 1]
  Z_post[, t]   <- step$genotype_next[, 1]
}

###### playing and testing with simulation
# Genotype for RR over time only
# is_RR_all <- (rowSums(1 - L) == ncol(L)) & (rowSums(1 - R) == ncol(R))  # all loci RR
# g_RR <- which(is_RR_all)                    
# RR_geno_freq <- Z_post[g_RR, ]              
# plot(RR_geno_freq, type = "b", ylim = c(0,1),
#      xlab = "time t", ylab = "Pr(RR…RR genotype)",
#      main = "Fully resistant genotype over time")
# # posterior resistant-allele frequency
# R_count <- t(Z_post) %*% ((1 - L) + (1 - R))   
# R_freq  <- R_count / 2                          
# # Plot (type ="b": both lines and point connected
# # col = 1 to n give us differents colors
# plot(RR_geno_freq, type= "b", ylim= c(0,1), col = 2,
#      xlab= "time t", ylab= "Probability",
#      main= "RR…RR genotype & resistant-allele frequency")
### now Allele frequency from RR
# lines () is a function that allow us the overlay the new plot on the old plo
# lty is the line type, here lty = 2 : give us a curve in dash
# lty could go from 1 to 6 (noted in Exercise Bayes)
# lines(R_freq[, 1], lty=2, lwd = 2, col =3)
# # here we have two or more than two locus
# if (ncol(R_freq) >= 2) lines(R_freq[, 2], lty=3)   
# # in legend, bty could be = "n" when we don't want the legend to be in box, 
# # and "o" if we wish to have it on a box
# legend( "topleft",
#        c("Pr(RR…RR genotype)", "Resistant allele (locus 1)", "Resistant allele (locus 2)"),
#        lty=c(1,2,3), bty="n")

# genotype next step per locus for homo resistant genotype
## Plot: 3 side-by-side panels
op <- par(mfrow = c(1, 3), mar = c(4,4,2,1))
cols <- seq_len(ncol(L))

 RR_locus <- (L == 0) & (R == 0)                 
RR_per_locus <- t(Z_post) %*% RR_locus          
# to plot a column of matrice, 
matplot(RR_per_locus, type="l", lty=1, lwd=2, lend = par("lend"), ylim=c(0,1),
        xlab="time t", ylab="Pr(RR at locus)",
        main="Posterior RR genotype per locus")
legend("topleft", legend=paste0("locus_", 1:ncol(L)), col=1:ncol(L), lty=1, bty="n")
# homo_suceptible genotype
SS_locus <- (L==1) & (R == 1) 
SS_per_locus <- t(Z_post) %*% SS_locus 
matplot(SS_per_locus, type="l", lty=1, lwd = 2, ylim=c(0,1),
        xlab="time t", ylab="Pr(SS at locus)",
        main="Posterior SS genotype per locus")
legend("topleft", legend=paste0("locus_", 1:ncol(L)), col=1:ncol(L), lty=1, bty="n")
# hetero_resistant genotype
SR_locus <- (L == 1) & (R == 0)  
SR_per_locus <- t(Z_post) %*% SR_locus 
matplot(SR_per_locus, type="l", lty=1, lwd = 2, ylim=c(0,1),
        xlab="time t", ylab="Pr(SR at locus)",
        main="Posterior SR genotype per locus")
legend("topleft", legend=paste0("locus_", 1:ncol(L)), col=1:ncol(L), lty=1, bty="n")
## Quick check: per time & locus, RR+SR+SS ≈ 1
chk <- SR_locus  + SS_locus + RR_locus
chk

par(op)
# allele frequency from our posterior and
new_allele_freq_rr <- (t(Z_post) %*% ((1 - L) + (1 - R))) / 2 
new_allele_freq_sr <- (t(Z_post)  %*% (L* (1 - R)) / 2)
new_allele_freq_ss <- (t(Z_post)  %*% (L * R) / 2)
# let's compare by HWE and our new allele frequency for heterozygote
hetero_SR_hwe   <- 2 * new_allele_freq_rr * (1 - new_allele_freq_rr)
# plot by using matplot 
matplot(hetero_SR_hwe, type="l", lty=1, ylim=c(0,1),
        xlab="time t", ylab="Pr(SR at locus hwe)",
        main="Posterior heterozygote allele frequency and heterozygote from HWE")
matlines(new_allele_freq_sr, lty=2, ylim= c(0,1))


## Phenotype Q* — Beta–Binomial link  
## We model between-batch heterogeneity in the susceptible phenotype.
## - p   : average probability of being susceptible (mean mortality) across batches.
## - rho : overdispersion / intraclass correlation (how much batches differ beyond Binomial).
## The latent batch-level probability θ ~ Beta(alpha, beta) with:
##   m = alpha + beta (precision; larger m = less heterogeneity),heterogeneity mean
## that there is difference between the conditions where our sample are from (temperature,
## locality, preicipitation, ....)
##   rho = 1 / (m + 1)  =>  m = 1/rho - 1,
##   alpha = p * m,  beta = (1 - p) * m.

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
## SIMULATION
## - Builds G genotypes over B loci (L,R) and time-varying genotype frequencies Z[,t] (sum_g Z_{g,t}=1).
## - Locus effects (θ_{b,c}) and dominance (h_b) define, for each genotype g and insecticide c,
##     U_{g,c} = v_s[c] * ∏_b u_{g b c}
##   = probability that genotype g is SUSCEPTIBLE under c (per-genotype “penetrance”).
## - The baseline susceptible mortality v_s[c]  over time: v_s^{(t)}[c]
##   enforcing a trend (e.g., decreasing susceptibility across time).
## - Overall susceptible phenotype at time t:
##     Q*_{t,c} = Σ_g Z_{g,t} * U_{g,c}^{(t)}
## - Observed deaths are drawn with Beta–Binomial:
##     X_{t,c} ~ BetaBinomial(N, p=Q*_{t,c}, ρ)
##   where ρ captures BETWEEN-batch heterogeneity (θ varies across batches); ρ=0 reduces to Binomial.
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
# # plotting q_star acroos time
# plot_qstar <- function(df, one_minus = FALSE) {
#   ylab <- if (one_minus) "1 − Q*" else "Q*"
#   ymap <- if (one_minus) (1 - df$Qstar) else df$Qstar
#   p <- ggplot2::ggplot(transform(df, y = ymap),
#                        ggplot2::aes(x = as.numeric(Time), y = y,
#                                     color = InsecticideName, group = InsecticideName)) +
#     ggplot2::geom_line(linewidth = 1.2) +
#     ggplot2::geom_point(size = 2) +
#     ggplot2::coord_cartesian(ylim = c(0, 1)) +
#     ggplot2::labs(x = "Time", y = ylab, color = "Insecticide") +
#     ggplot2::theme_minimal(base_size = 13) +
#     ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())
#   print(p); invisible(p)
# }
# # simulation with values
sim <- simulate_fake_data(
 G = 9, B = 2, Ttime = 8, C = 2,
  N_per_cell = 120, rho = 0.03, 
 v_s_decay = 0.97   
)
df <- sim$data
head(df)

# Plot Q* (or set one_minus = TRUE for 1 − Q*)
# plot_qstar(df, one_minus = FALSE)

# #  our prior allele frequency vs the posterior
# P_prior <- t(p_series)
# matplot(P_prior, type="l", lty=2, ylim=c(0,1),
#         xlab="time t", ylab="Allele frequency",
#         main="Resistant allele: prior (p_t) vs posterior (from Z_post)")
# matlines(R_freq, lty=1)
# legend("bottomright",
#        c(paste0("prior locus ", 1:ncol(L)), paste0("posterior locus ", 1:ncol(L))),
#        lty=c(rep(2, ncol(L)), rep(1, ncol(L))), bty="n", cex=0.8)

# ## Build indicators once
# RR_locus <- (L == 0L) & (R == 0L)   # G x B
# SR_locus <- (L == 1L) & (R == 0L)   # G x B   
# SS_locus <- (L == 1L) & (R == 1L)   # G x B
# 
# ## we are taking the genotype at each locus by transposing 
# RR_per_locus <- t(Z_post) %*% RR_locus
# SR_per_locus <- t(Z_post) %*% SR_locus
# SS_per_locus <- t(Z_post) %*% SS_locus
# 
# ## allele frequency
# new_allele_freq<- (t(Z_post) %*% ((1 - L) + (1 - R))) / 2 


# or 
# matplot(RR_per_locus, type = "l", lty = 1, col = cols,
#         ylim = c(0,1), xlab = "time t", ylab = "Pr(RR at locus)",
#         main = "RR per locus")
# matlines()
# legend("topleft", legend = paste0("locus_", 1:ncol(L)),
#        col = cols, lty = 1, bty = "n", cex = 0.9)
# 
# matplot(SR_per_locus, type = "l", lty = 1, col = cols,
#         ylim = c(0,1), xlab = "time t", ylab = "Pr(SR at locus)",
#         main = "SR per locus")
# legend("topleft", legend = paste0("locus_", 1:ncol(L)),
#        col = cols, lty = 1, bty = "n", cex = 0.9)
# 
# matplot(SS_per_locus, type = "l", lty = 1, col = cols,
#         ylim = c(0,1), xlab = "time t", ylab = "Pr(SS at locus)",
#         main = "SS per locus")
# legend("topleft", legend = paste0("locus_", 1:ncol(L)),
#        col = cols, lty = 1, bty = "n", cex = 0.9)

# phenotype frequency
## Compute P*_{l,t,c} from (L, R, Z, theta, h, v_s)
## L, R: G x B (0/1)
## Z:    G x L x T  (genotype frequencies per location l and time t)
## theta: B x C
## h:     length B
## v_s:   length C

# P*_l,t,c from (L, R, Z, theta, h, v_s)
# L, R: G x B (0/1); B = n_loci
# Z:    G x Lloc x T  (genotype frequency cube)
# theta: B x C        (locus x insecticide)
# h:     length B     (dominance per locus)
# v_s:   length C     (baseline susceptible mortality per insecticide)
# before calling compute_Pstar(...)
# compute_Qstar <- function(L, R, z, theta, h, v_s) {
#   stopifnot(is.matrix(L), is.matrix(R), all(dim(L) == dim(R)))
#   G <- nrow(L); B <- ncol(L)
#   stopifnot(is.matrix(theta) && nrow(theta) == B)
#   C <- ncol(theta); stopifnot(length(h) == B, length(v_s) == C)
#   
#   # ----- per-genotype × insecticide utility U_gc (G × C) -----
#   SS <- L * R
#   RR <- (1L - L) * (1L - R)
#   SR <- L * (1L - R)
#   theta_3 <- array(theta, c(1,B,C)); theta_3 <- theta_3[rep(1,G), , , drop=FALSE]
#   h_3     <- array(h,     c(1,B,1)); h_3     <- h_3[rep(1,G), , rep(1,C), drop=FALSE]
#   SS_3 <- array(SS, c(G,B,1))[, , rep(1,C), drop=FALSE]
#   RR_3 <- array(RR, c(G,B,1))[, , rep(1,C), drop=FALSE]
#   SR_3 <- array(SR, c(G,B,1))[, , rep(1,C), drop=FALSE]
#   eps  <- .Machine$double.eps
#   u_3  <- SS_3 + RR_3 * theta_3 + SR_3 * (h_3 * theta_3 + (1 - h_3))
#   U_gc <- exp(apply(log(pmax(u_3, eps)), c(1,3), sum))      # product over loci
#   U_gc <- sweep(U_gc, 2, v_s, "*")                          # scale by v_s[c]
#   
#   # ----- combine with Z over genotypes -----
#   if (length(dim(z)) == 2L) {
#     # Z is G × T  -> return T × C
#     stopifnot(nrow(z) == G)
#     Ttime <- ncol(z)
#     Qstar <- matrix(NA_real_, nrow=Ttime, ncol=C)
#     for (t in seq_len(Ttime)) Qstar[t, ] <- as.numeric(crossprod(z[, t], U_gc))
#     return(pmin(pmax(Qstar, 0), 1))
#   } else {
#     # Z is G × L × T  -> return L × T × C
#     stopifnot(dim(z)[1] == G)
#     Lloc <- dim(z)[2]; Ttime <- dim(z)[3]
#     Qstar <- array(NA_real_, c(Lloc, Ttime, C))
#     for (t in seq_len(Ttime)) Qstar[, t, ] <- t(z[, , t]) %*% U_gc
#     return(pmin(pmax(Qstar, 0), 1))
#   }
# }
# 
# # PRIOR
# # Resistance per locus & insecticide: theta (B x C), dominance h (length B),
# # baseline susceptible mortality v_s (length C)
# C <- 2
# theta <- matrix(runif(B * C), nrow = B, ncol = C)  # 2 x 2
# h     <- runif(B, 0, 1)                                      # 
# v_s   <- runif(C)                                  # length 2
# 
# # Compute P* (Lloc x Ttime x C) -> here (1 x 3 x 2)
# Qstar <- compute_Qstar(L, R, z, theta, h, v_s)
# Qstar
# 
# simulate_fake_data <- function(
    #     G = 9, B = 2, Lloc = 3, Ttime = 10, C = 2,
#     N_per_cell = 100, rho = 0.05, seed = 1,
#     .assign_df = NULL, .envir = parent.frame()
# ) {
#   set.seed(seed)
#   L <- matrix(sample(0:1, G*B, TRUE), nrow=G)
#   R <- matrix(sample(0:1, G*B, TRUE), nrow=G)
#   
#   Z <- array(NA_real_, c(G, Lloc, Ttime))
#   for (t in 1:Ttime) for (l in 1:Lloc) {
#     v <- runif(G); Z[, l, t] <- v / sum(v)
#   }
#   
#   theta <- matrix(runif(B*C, 0.3, 0.9), nrow=B, ncol=C)
#   h     <- runif(B, 0, 1)
#   v_s   <- runif(C, 0.6, 1.0)
#   
#   Qstar <- compute_Qstar(L, R, Z, theta, h, v_s)
#   
#   df <- expand.grid(Location=seq_len(Lloc), Time=seq_len(Ttime), Insecticide=seq_len(C),
#                     KEEP.OUT.ATTRS=FALSE)
#   # Correctly index Qstar to avoid flat lines
#   df$Qstar <- Qstar[cbind(df$Location, df$Time, df$Insecticide)]
#   df$InsecticideName <- paste0("Ins_", df$Insecticide)
#   df$TimeLabel <- df$Time
#   df <- tibble::as_tibble(df)
#   
#   df$N <- N_per_cell
#   df$Deaths <- mapply(function(n,q) rbetabinom1(n,q,rho), df$N, df$Qstar)
#   df$Observed <- df$Deaths / df$N
#   
#   if (!is.null(.assign_df)) assign(.assign_df, df, envir=.envir)
#   
#   list(
#     data = df,
#     Qstar = Qstar,
#     L = L, R = R, Z = Z,                 # <- Z (not Z_post)
#     theta = theta, h = h, v_s = v_s,
#     params = list(G=G, B=B, Lloc=Lloc, Ttime=Ttime, C=C,
#                   N_per_cell=N_per_cell, rho=rho, seed=seed)
#   )
# }
# 
# sim <- simulate_fake_data(G = 9, B = 2, #Lloc = 2, 
#                           Ttime = 8, C = 2,
#                           N_per_cell = 120, rho = 0.0001, seed = 123)
# df <- sim$data
# 
# 
# ggplot(df, aes(x = Time, y = Qstar, color = InsecticideName)) +
#   geom_line() + geom_point() +
#   facet_wrap(~ InsecticideName) +
#   coord_cartesian(ylim = c(0, 1)) +
#   labs(x = "Time", y = "Q*", color = "Insecticide")
# 
# 
# pretty_qstar <- function(df, one_minus = FALSE) {
#   ycol <- if (one_minus) 1 - df$Qstar else df$Qstar
#   ylab <- if (one_minus) "1 − Q*" else "Q*"
#   
#   ggplot(
#     transform(df, y = ycol),
#     aes(x = as.numeric(Time), y = y, color = InsecticideName, group = InsecticideName)
#   ) +
#     geom_line(linewidth = 1.2) +
#     geom_point(size = 2.4, alpha = 0.9) +
#     scale_x_continuous(breaks = sort(unique(as.numeric(df$Time))), expand = expansion(add = 0.1)) +
#     scale_y_continuous(limits = c(0, 1)) +
#     labs(x = "Time", y = ylab, color = "Insecticide") +
#     theme_minimal(base_size = 14) +
#     theme(
#       legend.position = "top",
#       panel.grid.minor = element_blank(),
#       panel.grid.major.x = element_blank(),
#       axis.title.x = element_text(margin = margin(t = 6)),
#       axis.title.y = element_text(margin = margin(r = 6))
#     )
# }
# 
# # Q*
# pretty_qstar(df)
# 
# # 1 − Q*
# pretty_qstar(df, one_minus = TRUE)

# Plot if we will insert location
# plot_lines <- function(df) {
#   ggplot(df, aes(x = as.numeric(TimeLabel),
#                  y = Qstar,
#                  color = InsecticideName,
#                  group = InsecticideName)) +
#     geom_line() +
#     geom_point() +
#     facet_wrap(~ LocationName, ncol = 2) +
#     coord_cartesian(ylim = c(0, 1)) +
#     labs(x = "Time", y = "Q*", color = "Insecticide",
#          title = "Q* over time by location") +
#     theme_minimal(base_size = 12)
# }
# 
# # usage
# plot_lines(df) 
# 
# # Inspect the tidy fake dataset (one row per Location x Time x Insecticide)
# head(sim$data)
# 
# 
# df$OneMinusQ <- 1 - df$Qstar
# ggplot(df, aes(as.numeric(TimeLabel), OneMinusQ, color = InsecticideName, group = InsecticideName)) +
#   geom_line() + geom_point() + facet_wrap(~ LocationName, ncol = 2) +
#   coord_cartesian(ylim = c(0,1)) + labs(y = "1 − Q*")
# compute_Qstar <- function(L, R, z, theta, h, v_s) {
#   stopifnot(is.matrix(L), is.matrix(R), all(dim(L) == dim(R)))
#   G <- nrow(L); B <- ncol(L)
#   stopifnot(is.matrix(theta) && nrow(theta) == B)
#   C <- ncol(theta); stopifnot(length(h) == B, length(v_s) == C)
#   
#   # ----- per-genotype × insecticide utility U_gc (G × C) -----
#   SS <- L * R
#   RR <- (1L - L) * (1L - R)
#   SR <- L * (1L - R)
#   
#   theta_3 <- array(theta, c(1,B,C)); theta_3 <- theta_3[rep(1,G), , , drop = FALSE]
#   h_3     <- array(h,     c(1,B,1)); h_3     <- h_3[rep(1,G), , rep(1,C), drop = FALSE]
#   
#   SS_3 <- array(SS, c(G,B,1))[, , rep(1,C), drop = FALSE]
#   RR_3 <- array(RR, c(G,B,1))[, , rep(1,C), drop = FALSE]
#   SR_3 <- array(SR, c(G,B,1))[, , rep(1,C), drop = FALSE]
#   
#   eps  <- .Machine$double.eps
#   u_3  <- SS_3 + RR_3 * theta_3 + SR_3 * (h_3 * theta_3 + (1 - h_3))
#   U_gc <- exp(apply(log(pmax(u_3, eps)), c(1,3), sum))   # G × C
#   U_gc <- sweep(U_gc, 2, v_s, "*")                       # apply baseline
#   
#   # ----- combine with z over genotypes -----
#   if (is.null(dim(z))) {
#     # z: length G
#     stopifnot(length(z) == G)
#     out <- as.numeric(crossprod(z, U_gc))                # 1 × C
#     return(pmin(pmax(out, 0), 1))
#   } else if (length(dim(z)) == 2L) {
#     # z: G × T  -> T × C
#     stopifnot(nrow(z) == G)
#     Ttime <- ncol(z)
#     Q <- matrix(NA_real_, nrow = Ttime, ncol = C)
#     for (t in seq_len(Ttime)) Q[t, ] <- as.numeric(crossprod(z[, t], U_gc))
#     return(pmin(pmax(Q, 0), 1))
#   } else {
#     # z: G × L × T  -> L × T × C
#     stopifnot(dim(z)[1] == G)
#     Lloc <- dim(z)[2]; Ttime <- dim(z)[3]
#     Q <- array(NA_real_, c(Lloc, Ttime, C))
#     for (t in seq_len(Ttime)) Q[, t, ] <- t(z[, , t]) %*% U_gc
#     return(pmin(pmax(Q, 0), 1))
#   }
# }
# 
# simulate_fake_data <- function(
    #     G = 9, B = 2, Lloc = 3, Ttime = 10, C = 2,
#     N_per_cell = 100, rho = 0.05, seed = 1,
#     .assign_df = NULL, .envir = parent.frame()
# ) {
#   set.seed(seed)
#   
#   # Force scalar integers (prevents seq_len() complaints)
#   Lloc <- as.integer(Lloc)[1]
#   Ttime <- as.integer(Ttime)[1]
#   C <- as.integer(C)[1]
#   G <- as.integer(G)[1]
#   B <- as.integer(B)[1]
#   
#   L <- matrix(sample(0:1, G*B, TRUE), nrow = G)
#   R <- matrix(sample(0:1, G*B, TRUE), nrow = G)
#   
#   Z <- array(NA_real_, c(G, Lloc, Ttime))
#   for (t in seq_len(Ttime)) for (l in seq_len(Lloc)) {
#     v <- runif(G); Z[, l, t] <- v / sum(v)
#   }
#   
#   theta <- matrix(runif(B*C, 0.3, 0.9), nrow = B, ncol = C)
#   h     <- runif(B, 0, 1)
#   v_s   <- runif(C, 0.6, 1.0)
#   
#   Qstar <- compute_Qstar(L, R, Z, theta, h, v_s)  # L × T × C
#   
#   df <- expand.grid(Location = seq_len(Lloc),
#                     Time     = seq_len(Ttime),
#                     Insecticide = seq_len(C),
#                     KEEP.OUT.ATTRS = FALSE)
#   df$Qstar <- Qstar[cbind(df$Location, df$Time, df$Insecticide)]
#   df$InsecticideName <- paste0("Ins_", df$Insecticide)
#   df <- tibble::as_tibble(df)
#   
#   df$N <- N_per_cell
#   df$Deaths <- mapply(function(n, q) rbetabinom1(n, q, rho), df$N, df$Qstar)
#   df$Observed <- df$Deaths / df$N
#   
#   if (!is.null(.assign_df)) assign(.assign_df, df, envir = .envir)
#   
#   list(
#     data = df,
#     Qstar = Qstar,
#     L = L, R = R, Z = Z,
#     theta = theta, h = h, v_s = v_s,
#     params = list(G=G, B=B, Lloc=Lloc, Ttime=Ttime, C=C,
#                   N_per_cell=N_per_cell, rho=rho, seed=seed)
#   )
# }


