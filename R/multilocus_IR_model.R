set.seed(123)
## Prior
# p is the initial allele frequency

n_loci <- 2
# call the dummy function
source("C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/R/create_dummy_matrices.R")
genotype_combination <- create_dummy_matrices(n_loci=2)

L <- genotype_combination$left
L
R <- genotype_combination$right
L[3, 2]
## Multilocus genotype frequencies
# we will use the dummy matrice of the left and right hand of the allele, 
# we define a function F which is the genotype probability across the n_loci
# sweep(x, MARGIN, STATS, FUN) : it allows Sylviane to operate value in matrix,
# x= matrix name, margin = to perform on Row = 1, on column =2, stats = the value
# to use (formula), function is the * or / or +, 
## per-locus probability that matches the genotype of g 
probability_genotype <- function(p){
  G <- nrow(L)
  B <- ncol(L)
  if(length(p) != B){stop("locus probability dimension mismatch")}
  F <- matrix(data=NA, nrow = G, ncol= B)
  for (g in seq_len(G)) {
    for (lo in seq_len(B)) {
      prob_left <- (1-p[lo]) * L[g,lo] + p[lo]* (1-L[g,lo])
      prob_right <- (1-p[lo]) * R[g,lo] + p[lo]* (1-R[g,lo])
      heterozygote_duplicates <- 1 + L[g,lo] - R[g,lo]
      F[g,lo] <- prob_left * prob_right * heterozygote_duplicates
      
    }
    z[g] <- prod(F[g, ])
    # # or this:
    # z[g] <- 1
    # for(lo in seq_len(B)) {
    #   z[g] <- z[g] * F[g, lo]
    # }
  }
  z[g] <- prod(F[g, ])
}
p <- runif(B)

### fast way to get genotype perlocus
probability_genotype_fast <- function(p){
  G <- nrow(L)
  B <- ncol(L)
  if(length(p) != B){stop("locus probability dimension mismatch")}
  F <- matrix(data=NA, nrow = G, ncol= B)
  
  prob_left <- sweep(L, 2, 1 - p, FUN = "*") +
    sweep(1 - L, 2, p, FUN = "*")
  
  prob_right <- sweep(R, 2, 1 - p, FUN = "*") +
    sweep(1 - R, 2, p, FUN = "*")
  
  heterozygote_duplicates <- 1 + L - R
  
  F <- prob_left * prob_right * heterozygote_duplicates
  z <- apply(F, 1, prod)

  return(z) 
}


# p : allele frequency at each locus 
p <- runif(B) 
# genotype frequency as p is a matrix and returning a vector z
probability_genotype_fast <- function(p){
  stopifnot(is.matrix(L), is.matrix(R), all(dim(L) == dim(R)))
  G <- nrow(L); B <- ncol(L)
  stopifnot(length(p) == B)
  
  prob_left  <- sweep(L, 2, 1 - p, "*") + sweep(1 - L, 2, p, "*")
  prob_right <- sweep(R, 2, 1 - p, "*") + sweep(1 - R, 2, p, "*")
  dup <- 1 + L - R   # =2 pour (1,0), =1 sinon
  
  F <- prob_left * prob_right * dup
  z_unnorm <- apply(F, 1, prod)
  z <- z_unnorm / sum(z_unnorm)  # vecteur longueur G
  z
}
# to test our function
test <- probability_genotype_fast(p)
test

## selection at locus level
## the effect of intervention 
K <- 4
probability_genotype_fast <- function(p){
  G <- nrow(L); B <- ncol(L)
  if (length(p) != B) stop("locus probability dimension mismatch")
  
  # proba d'observer l'allèle encodé en L et R (colonne par colonne)
  prob_left  <- sweep(L, 2, 1 - p, "*") + sweep(1 - L, 2, p, "*")  # L==1 -> 1-p ; L==0 -> p
  prob_right <- sweep(R, 2, 1 - p, "*") + sweep(1 - R, 2, p, "*")  # R==1 -> 1-p ; R==0 -> p
  
  # facteur de doublon (2 pour hétérozygote gardé (1,0), 1 sinon)
  dup <- 1 + L - R
  
  # probabilité par locus, puis produit sur les loci -> masse par génotype (ligne)
  F <- prob_left * prob_right * dup
  z_unnorm <- apply(F, 1, prod)     # <- longueur G (une par combinaison L/R)
  
  # normalisation -> vecteur de proba
  s <- sum(z_unnorm)
  if (!is.finite(s) || s == 0) stop("Somme nulle/non-finie pour z; vérifiez L/R/p")
  z <- z_unnorm / s
  return(z)                          # <-- VECTEUR longueur G
}

# Xn <- c(0.1, 0.2,0.3,0.4)
X <- runif(K)
betamat <- matrix(rnorm(B * K), nrow=B, ncol=K)
gammax <- exp(betamat)
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



## GENOTYPE NEXT STEP
# step to get the genotype frequency at next step without using time series
polygenic_multilocus_next_step <- function(w, h){
  stopifnot(is.matrix(L), is.matrix(R), all(dim(L) == dim(R)))
 # define G and B 
  G <- nrow(L); B <- ncol(L)
  stopifnot(length(w) == B, length(h) == B)
  stopifnot(exists("z"), length(z) == G)
 # define an empty matrix to put the value of the fitness polygenic in depend of 
  # the zygosity (homo - hetero) of our hand of alleles
  Gw <- matrix(NA_real_, nrow = G, ncol = B)
  # define the combination SS and RR and RS to give more sense to the calculation
  for (lo in seq_len(B)) {
    SS  <- (L[, lo] == 1L) & (R[, lo] == 1L)      # (1,1)
    RR  <- (L[, lo] == 0L) & (R[, lo] == 0L)      # (0,0)
    SR  <- (L[, lo] == 1L) & (R[, lo] == 0L)      # (1,0)
    
    Gw[, lo] <- 1 * SS +
      w[lo] * RR +
      (h[lo] * w[lo] + (1 - h[lo])) * SR
  }
  
  r_vec  <- apply(Gw, 1, prod)
  unpost <- z * r_vec
  post   <- unpost / sum(unpost)
  
  list(Gw = Gw,
       r = matrix(r_vec,  nrow = G, ncol = 1),
       post_next = matrix(post, nrow = G, ncol = 1),
       z = matrix(z,     nrow = G, ncol = 1))
}











