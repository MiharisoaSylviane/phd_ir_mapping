# mapping from the allele fequency of different markers into the survival proportion
# for one insecticide resistance
# load library
library(greta)
# library(tidyverse)
# n<-100
# data_template <- expand_grid(
#   allele = seq_len(3),
#   anopheles_indiv = seq_len(n),
#   # ratio_susceptible =0,
#   # ratio_resistant = 0,
#   time_series = seq_len(4),
#   allele_frequency = NA) %>%
#   mutate(
#     before_after = case_when(
#       time_series %in% 1:2 ~ "before",  # Time points 1 and 2 are "before"
#       time_series %in% 3:4 ~ "after"   # Time points 3 and 4 are "after"
#     )
#   )


# define the inital value of allele frequency
# set.seed(123)
# data_template <- data_template %>%
#   mutate(
#     allele_frequency = ifelse(
#       before_after == "before",
#       runif(n()),  
#       allele_frequency
#     )
#   )

# time series
time_series <- 10

# itns_covariates <- seq(0.1, 0.9, length.out = time_series)
intervention_times <- 4
itns_covariates <- c(rep(0.9, length.out = intervention_times),
                     rep(0.0, length.out = time_series - intervention_times))

# # "before" intervention, 10% of population using nets. "after" intervention, 90%
# # using nets
# itns_covariates <- c(0.1, 0.1, 0.9, 0.9)


# Define greta model
# Prior for proportion of susceptible
p0 <- uniform(0, 1)


# # Prior for selection pressure average
# s <- normal(0, 1, truncation =c(0,1))  

# coefficient of itns effect resistance
# beta <- log(gammax)

beta <- normal(0,1)
gammax <- exp(beta)

 
# selection pressure
st <- gammax*itns_covariates
 
# # relative fitness of mosquitoes with allele resistant R/S
w <- 1 + st

# calculate(w, nsim = 1)
 
  
# Simulate allele frequencies
# to assign it as varying with time
# p <- numeric(time_series)

# # proportion without resistance p where we will define the fraction of resistant
# for (t in 1:(time_series)) {
#   if (t >2) {
# p <- p[t] / p[t] + (1 - p[t]) * w 
# p_resist <- 1 - p # proportion of resistant
#   }
# }


# Probability of an allele to occur
p_susc <- zeros(time_series)
# p_susc <- as_data(rep(0, 4))
p_susc[1] <- p0
for (t in 1:(time_series-1)) {
  p_susc[t + 1] <- p_susc[t]/ (p_susc[t] + (1- p_susc[t]) * w[t])
  
}

p_resist <- 1 - p_susc


# simulate genotypic allele frequency data
n_tested <- rep(100, time_series)
n_positive <- binomial(size = n_tested, prob = p_resist)

# distribution(n_positive) <- binomial(size = n_tested, prob = p_resist)


observed_frequency <- n_positive / n_tested
sim<- calculate(p_resist, nsim=1)


sims <- calculate(gammax, p_resist,
                  n_positive,
                  observed_frequency,
                  nsim = 1,
                  values = list(p0 = 0.9))
sims

plot(sims$p_resist[1, , 1],
     ylim = c(0, 1),
     type = "b")

points(sims$observed_frequency,
       pch = 16)


### Sylviane has to determine the phenotype proportion by using the allele frequency

resistant_phenotype <- p_resist^2 + 2 * p_susc * p_resist

# simulation
sims <- calculate(gammax, p_resist, resistant_phenotype,
                  n_positive, observed_frequency,
                  nsim = 1,
                  values = list(p0 = 0.9))

sims

plot(sims$resistant_phenotype [1, , 1],
     ylim = c(0, 1),
     type = "b")


# Cases when we have two alleles resistant in the population, homozygous allele and heterygous allele
# we need to find the mean_fitness
# degree of dominance h= 0 if resistant allele recessive
# h= 0.5 if no dominance then additive
# h = 1 if allele completely dominant
library(tidyverse)
h <- uniform(0,1)


p <- function(p, s, h) {
  q <- 1 - p
  numerator <- p^2 * (1 + s) + p * q * (1 + h * s)
  denominator <- 1 + s * (p^2 + 2 * h * p * q)
  p_next <- numerator / denominator
  return(p_next)
}


p_resist <- 1 - p_susc


# simulate genotypic allele frequency data
n_tested <- rep(100, time_series)
n_positive <- binomial(size = n_tested, prob = p_resist)

observed_frequency <- n_positive / n_tested
sim<- calculate(p_resist,h, w, nsim=1)
sim

sims <- calculate(gammax, p_resist,
                  n_positive,
                  observed_frequency,
                  nsim = 1,
                  values = list(p0 = 0.9))

plot(sims$p_resist[1, , 1],
     ylim = c(0, 1),
     type = "b")


# the proportion of phenotype with homozygous allele
# the phenotype resistant is expressed 
resistant_phenotype_hom <- p_resist^2

# the phenotype susceptible is expressed 
susceptible_phenotype_hom <- p_susc^2

# the proportion of the phenotype with heterozygous allele
resistant_phenotype_heter <- 2*p_resist*p_susc 

# here we got the proportion of resistant phenotype 
proportion_resistant_phenotype <- resistant_phenotype_hom * 1 +
  susceptible_phenotype_hom* 0 +
  resistant_phenotype_heter * h


sims <- calculate(proportion_resistant_phenotype, p_resist,
                  n_positive,
                  observed_frequency,
                  nsim = 1,
                  values = list(p0 = 0.9))
sims


# we consider we have two alleles resistant with the susceptible allele
# p+q+r2 =1
# p represent the allele frequency of the first allele resistant R
# r represent the allele frequency of the second allele resistant Q in the locus
# q represent the allele frequency of the allele susceptible
# for each, the selection and fitness is different
# Assuming beta is a vector of coefficients for each allele's fitness

# number of alleles 
n_allele <- 3

beta <- normal(0,1, dim = n_alleles)
# replace with your actual beta values

# calculate relative fitness for each allele
gammax <- exp(beta)  # This will be a vector of length 4

# ITN coverage covariates - assuming you have one for each allele
itns_covariates <- c(0.1, 0.7, 0.5) 

# calculate selection pressure for each allele
st <- gammax * itns_covariates  


# relative fitness of each allele, we make assumption here
# Parameters
w <- 1 + st 
time_series <- 100


# Initial frequency of susceptible allele (S)\
p0 <- uniform(0,1)

# Initial frequency of resistant allele 1 (R1)

q0 <- 1 - p0

# Initial frequency of resistant allele 1 (R2)
r0 <- 1-p0-q0

# greta array to get the initial proportion
p <- zeros(time_series)  
r <- zeros(time_series) 


# Set first column to p0 (greta handles this as a graph operation)
p[1] <- p0
r[1] <- p0

# iterate through time series

mfunctiom<- for (t in 1:(time_series - 1)) {
  w_bar <-  q[t]*w[2] + r[t]*w[3]
  
  # Next generation frequencies

  p[t+1] <- (q[t] * w[2]) / w_bar
  r[t+1] <- (r[t] * w[3]) / w_bar
}

resist <- p+r
q <- 1 - resist



# Plot results
n_tested <- rep(100, time_series)
n_positive <- binomial(size = n_tested, prob =q)
observed_frequency <- n_positive / n_tested
sim1 <- calculate(n_positive, q, nsim=1)
sim1


library(greta)

# Parameters
n_alleles <- 3
time_series <- 100
n_tested <- rep(100, time_series)

# Initial allele frequencies
p0 <- uniform(0, 1) 
q0 <- 1 - p0  
r0 <- 1 - p0 - q0  

# Selection coefficients
beta <- normal(0, 1, dim = n_alleles)
gammax <- exp(beta)

# ITN coverage (example values)
itns_covariates <- matrix(
  c(0.1, 0.7, 0.5)
)

# Calculate selection pressure
st <- gammax  * itns_covariates 
w <- 1 + st  

# Initialize frequency 
p <- zeros(time_series)
q <- zeros(time_series)
r <- zeros(time_series)

p[1] <- p0
q[1] <- q0
r[1] <- r0

# Define Time Dynamics
for (t in 1:(time_series - 1)) {
  w_bar <- p[t]*w + q[t]*w + r[t]*w

}


p <- (p * w) / w_bar
q<- (q * w) / w_bar
r <- (r * w) / w_bar
p
q
r

# simulation
n_tested <- rep(100, time_series)
n_positive <- binomial(size = n_tested, prob = p)
observed_frequency <- n_positive / n_tested
sim<- calculate(p, nsim=1)


sims <- calculate(gammax, p_resist,
                  n_positive,
                  observed_frequency,
                  nsim = 1,
                  values = list(p0 = 0.9))


  














# selection pressure 
calculate_resistance_frequency <- function() {
  # Initialize vector for susceptible allele frequencies
  p_susc <- numeric(time_series)
  p_susc[1] <- p0
  
  # Calculate fitness based on h
  w <- case_when(
    h  ==  0, w <- w,            
    h <= 0.5 &  h >= 0.75, w <-  1 - (h * s),
    h < 0.75, w <- 1 - s,          
    TRUE ~ 1                                  
  )


binomial_likelihood <- function(x, n, p) {
    choose(n, x) * p^x * (1-p)^(n-x)
}

