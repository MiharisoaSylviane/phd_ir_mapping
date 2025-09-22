library(greta)
##### 
set.seed(123)
## Prior
# p is the initial allele frequency

n_loci <- 5
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
p <- runif(B) 
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
G <- nrow(L);
B <- ncol(L)
T <- 10

# matrix storage
F_arr <- array(NA_real_, dim = c(G, B, T))  
z_mat <- matrix(NA_real_, nrow = G, ncol = T)
r_mat <- matrix(NA_real_, nrow = G, ncol = T)
Gw_arr <-array(NA_real_, dim = c(G, B, T))
p_series <- matrix(NA_real_, nrow = B, ncol = T)  
Z_post <- matrix(NA_real_, nrow = G, ncol = T)
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
for (t in 1:T) {
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
p_series <- matrix(NA_real_, nrow = B, ncol = T)
# giving some values for the initial allele frequencies and r
p_series[, 1] <- p0 <- c(0.10, 0.30, 0.4, 0.02, 0.03)        
r <- c(0.20, 0.35, 0.2, 0.1, 0.02)                          
# definining the loop for those allele frequency at t and calculate 
# fraction of each multilocus genotype F_t, unormailzed fraction of genotype frequency
# z_u, and the posterior genotype frequency z_t
for (t in 1:T) {
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

###### playing and testing the simulation
# Genotype for RR over time only
is_RR_all <- (rowSums(1 - L) == ncol(L)) & (rowSums(1 - R) == ncol(R))  # all loci RR
g_RR <- which(is_RR_all)                    
RR_geno_freq <- Z_post[g_RR, ]              
plot(RR_geno_freq, type = "b", ylim = c(0,1),
     xlab = "time t", ylab = "Pr(RR…RR genotype)",
     main = "Fully resistant genotype over time")

# posterior resistant-allele frequency
R_count <- t(Z_post) %*% ((1 - L) + (1 - R))   
R_freq  <- R_count / 2                          

# Plot (type ="b": both lines and point connected
# col = 1 to n give us differents colors
plot(RR_geno_freq, type="b", ylim=c(0,1), col = 2,
     xlab="time t", ylab="Probability",
     main="RR…RR genotype & resistant-allele frequency")
### now Allele frequency from RR
# lines () is a function that allow us the overlay the new plot on the old plo
# lty is the line type, here lty = 2 : give us a curve in dash
# lty could go from 1 to 6 (noted in Exercise Bayes)
lines(R_freq[, 1], lty=2, lwd = 2, col =3)
# here we have two or more than two locus
if (ncol(R_freq) >= 2) lines(R_freq[, 2], lty=3)   
# in legend, bty could be = "n" when we don't want the legend to be in box, 
# and "o" if we wish to have it on a box
legend("bottomright",
       c("Pr(RR…RR genotype)", "Resistant allele (locus 1)", "Resistant allele (locus 2)"),
       lty=c(1,2,3), bty="n")

# genotype next step per locus for homo resistant genotype
RR_locus <- (L == 0) & (R == 0)                 
RR_per_locus <- t(Z_post) %*% RR_locus          
# to plot a column of matrice, 
matplot(RR_per_locus, type="l", lty=1, lwd=3, lend = par("lend"), ylim=c(0,1),
        xlab="time t", ylab="Pr(RR at locus)",
        main="Posterior RR genotype per locus")
legend("bottomright", legend=paste0("locus_", 1:ncol(L)), col=1:ncol(L), lty=1, bty="n")
# homo_suceptible genotype
SS_locus <- (L==1) & (R == 1) 
SS_per_locus <- t(Z_post) %*% SS_locus 
matplot(SS_per_locus, type="l", lty=1, ylim=c(0,1),
        xlab="time t", ylab="Pr(SS at locus)",
        main="Posterior SS genotype per locus")
legend("bottomright", legend=paste0("locus_", 1:ncol(L)), col=1:ncol(L), lty=1, bty="n")
# hetero_resistant genotype
SR_locus <- (L == 1) & (R == 0)  
SR_per_locus <- t(Z_post) %*% SR_locus 
matplot(SR_per_locus, type="l", lty=1, ylim=c(0,1),
        xlab="time t", ylab="Pr(SR at locus)",
        main="Posterior SR genotype per locus")
legend("bottomright", legend=paste0("locus_", 1:ncol(L)), col=1:ncol(L), lty=1, bty="n")

# allele frequency from our posterior and
new_allele_freq_rr <- (t(Z_post) %*% ((1 - L) + (1 - R))) / 2 
new_allele_freq_sr <- (t(Z_post)  %*% (L* (1 - R)) / 2)
# let's compare by HWE and our new allele frequency for heterozygote
hetero_SR_hwe   <- 2 * new_allele_freq_rr * (1 - new_allele_freq)
# plot by using matplot 
matplot(hetero_SR_hwe, type="l", lty=1, ylim=c(0,1),
        xlab="time t", ylab="Pr(SR at locus hwe)",
        main="Posterior heterozygote allele frequency and heterozygote from HWE")
matlines(new_allele_freq_sr, lty=2, ylim= c(0,1))


# #  our prior allele frequency vs the posterior
# P_prior <- t(p_series)
# matplot(P_prior, type="l", lty=2, ylim=c(0,1),
#         xlab="time t", ylab="Allele frequency",
#         main="Resistant allele: prior (p_t) vs posterior (from Z_post)")
# matlines(R_freq, lty=1)
# legend("bottomright",
#        c(paste0("prior locus ", 1:ncol(L)), paste0("posterior locus ", 1:ncol(L))),
#        lty=c(rep(2, ncol(L)), rep(1, ncol(L))), bty="n", cex=0.8)


## Build indicators once
RR_locus <- (L == 0L) & (R == 0L)   # G x B
SR_locus <- (L == 1L) & (R == 0L)   # G x B   
SS_locus <- (L == 1L) & (R == 1L)   # G x B

## we are taking the genotype at each locus by transposing 
RR_per_locus <- t(Z_post) %*% RR_locus
SR_per_locus <- t(Z_post) %*% SR_locus
SS_per_locus <- t(Z_post) %*% SS_locus

## allele frequency
new_allele_freq<- (t(Z_post) %*% ((1 - L) + (1 - R))) / 2 

## Quick check: per time & locus, RR+SR+SS ≈ 1
chk <- RR_per_locus + SR_per_locus + SS_per_locus


## Plot: 3 side-by-side panels
op <- par(mfrow = c(1, 3), mar = c(4,4,2,1))
cols <- seq_len(ncol(L))
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

par(op)

                  
