set.seed(123)
## Prior
# p is the initial allele frequency
p <- runif(0,1)
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
  return(z) 
}
 
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

## the effect of intervention 
K <- 4
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
h <- runif(B, 0,1)
time_series <- 10
polygenic_multilocus_next_step <- function(w, h){
   G <- nrow(L)
   B <- ncol(L)
   Gw <- matrix(data=NA, nrow = G, ncol= B)
   if(length(w) != B){stop("locus probability dimension mismatch")}
   for (t in 1:(time_series)){
   for (g in seq_len(G)) {
     for (lo in seq_len(B)) {
       # polygenic fitness
       Gw[g, lo, t] <- 1 * R[g, lo, t] +
         w[lo, t] * (1 - L[g, lo, t]) +
         (w[lo,t] * h[lo,t] + 1 - h[lo,t]) * (L[g, lo,t] - R[g, lo,t])
         heterozygote_duplicates <- 1 + L - R
      r <- apply(Gw, 1, prod) 
     }
     unpost [g,t+1] <- z[g,t] * r[g,t]
   }
     post[g,t+1] / sum(post[g,t+1])
   }
}



# assumes globals: L, R, z













# multilocus genotype at next step with all parameters
p <- rbeta(B, 2, 2)         # prior per locus (safer than runif)

polygenic_multilocus_next_step<- function( w, h) {
   method <- match.arg(method)
   #checking for each matrix 
   stopifnot(is.matrix(L), 
             is.matrix(R), 
             all(dim(L) == dim(R)))
             G <- nrow(L)
             B <- ncol(L)
             stopifnot(length(w) == B, length(h) == B)
             if (is.null(dim(w))) w <- matrix(rep(w, T), nrow = B, ncol = T)
             if (is.null(dim(h))) h <- matrix(rep(h, T), nrow = B, ncol = T)
             Gw     <- array(NA_real_, dim = c(G, B, T))
             r      <- matrix(NA_real_, nrow = G, ncol = T)
             unpost <- matrix(NA_real_, nrow = G, ncol = T + 1)
             post   <- matrix(NA_real_, nrow = G, ncol = T + 1)
   # assign a matrix to 
   p <- matrix(p, G, B, byrow = TRUE)
   W <- matrix(w, G, B, byrow = TRUE)
   H <- matrix(h, G, B, byrow = TRUE)
   
   # genotype
   prob_left  <- L * (1 - p) + (1 - L) * p
   prob_right <- R * (1 - p) + (1 - R) * p
   z <- apply(prob_left * prob_right * (1 + L - R), 1, prod)
   
   # polygenic fitness
   Gw <- R + W * (1 - L) + (W * H + 1 - H) * (L - R)
   r <- apply(Gw, 1, prod)
   
   # posterior over genotypes
   unpost <- z * r
   post <- unpost/ sum(unpost)
 }

post <- posterior_genotype(w, h)
 
### Allele frequency at next step a(t+1) in depend of z



 



 

 
 
 
 
 
 
 
 
 
 
 
 
 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ## this is for gneotype and locus and times scale
time
probability_genotype <- function(p){
  G <- nrow(L)
  B <- ncol(L)
  p <- matrix(p, nrow = B, ncol = time_series) 
  Tn <- ncol(p)
  if(length(p) != B){stop("locus probability dimension mismatch")}
  F <- matrix(data=NA, dim = c(G, B, Tn)) #F <- matrix(data=NA, nrow = G, ncol= B)
  for (t in seq_len(Tn))
    for (g in seq_len(G)) {
      for (lo in seq_len(B)) {
        prob_left <- (1-p[lo,t]) * L[g,lo,t] + p[lo,t]* (1-L[g,lo,t])
        prob_right <- (1-p[lo,t]) * R[g,lo,t] + p[lo,t]* (1-R[g,lo,t])
        F[g,lo,t] <- prob_left * prob_right
      }
    }
  return(F) 
}













