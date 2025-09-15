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
##multilocus genotype frequencies
# we will use the dummy matrice of the left and right hand of the allele, 
# we define a function F which is the genotype probability across the n_loci
# sweep(x, MARGIN, STATS, FUN) : it allows Sylviane to operate value in matrix,
# x= matrix name, margin = to perform on Row = 1, on column =2, stats = the value
# to use (formula), function is the * or / or +, 
# 
 probability_genotype <- function(p){
   G <- nrow(L)
   B <- ncol(L)
   if(length(p) != B){stop("locus probability dimension mismatch")}
   F <- matrix(data=NA, nrow = G, ncol= B)
   for (g in seq_len(G)) {
     for (lo in seq_len(B)) {
    prob_left <- (1-p[lo]) * L[g,lo] + p[lo]* (1-L[g,lo])
    prob_right <- (1-p[lo]) * R[g,lo] + p[lo]* (1-R[g,lo])
     F[g,lo] <- prob_left * prob_right
     }
   }
   return(F) 
 }
 


 
