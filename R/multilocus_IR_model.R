set.seed(123)
## Prior
# p is the initial allele frequency
p <- runif(0,1)
n_loci <- 2
# call the dummy function
source("C:/Users/Sylviane/Desktop/training_perth_2024/phd_ir_mapping/R/create_dummy_matrices.R")
genotype_combination <- create_dummy_matrices()
##multilocus genotype frequencies
# we will use the dummy matrice of the left and right hand of the allele, 
# we define a function F which is the genotype probability across the n_loci
# sweep(x, MARGIN, STATS, FUN) : it allows Sylviane to operate value in matrix,
# x= matrix name, margin = to perform on Row = 1, on column =2, stats = the value
# to use (formula), function is the * or / or +
# 
 genotype_probabilty <- function(p) {
  #combination of the alleles
   comb_gen <- genotype_combination(n_loci)
  # A= (1-p)*L + p*(1-L)
  left_prob  <- sweep(left,  2, 1 - p, "*") + sweep(1 - left,  2, p, "*")
  # B= (1-p)*R + p*(1-R)
  right_prob <- sweep(right, 2, 1 - p, "*") + sweep(1 - right, 2, p, "*")
  # A*B
  F <- left_prob * right_prob
  return(F)
}
genotype_probabilty(n_loci = 2, p)
