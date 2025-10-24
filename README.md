# Spatiotemporal modelling of the spread of insecticide resistance in Africa intergrating genotype, allele frequency and phenotype
## Description
This repository maps insecticide resistance in Anopheles mosquitoes across Africa. It combines genotype data, phenotype (bioassay) outcomes, and allele-frequency measurements in a Bayesian hierarchical, semi-mechanistic framework to estimate resistance and susceptibility through space and time with quantified uncertainty.
The single locus version of this model could be found in https://github.com/idem-lab/ir_cube designed ([Golding](https://github.com/goldingn)).

<img width="720" height="1600" alt="image" src="https://github.com/user-attachments/assets/40b4dc8b-2a20-41aa-9ed8-8f4d74de6fc8" />

## Data inputs
- WHO discriminating concentration results from **Vector Atlas**, incorporating **IR Mapper** and **WHO Malaria Threats Map**
- **Genotypes**: Vector Atlas genotype frequency data (single-locus and multi-locus, where available) if avalaible
- **Covariates** (per location–time):
  - ITN use (modelled), (~timeframe to identify) (MAP)
  - IRS coverage, (~timeframe to identify) (WHO/MAP)
  - Population density, (~timeframe to identify)(WorldPop/MAP; UN-adjusted 2021–2030)
  - pesticide use (~timeframe and source to identify)


> Exact file paths and schemas will be documented as the code stabilizes

## Notation
We will develop and add notation through times of coding
 $l$: location; $t$: time; $c$: insecticide class
- Bioassay: $N$ tested, $D$ died
- Genotypes (multi-locus): for record $z$, $M^{z}$ tested, $N^{z}$ counts across genotype categories.
- Genotypes (single locus): allele-1 freq $p_{l,t}$ so $$q_{l,t}= 1 - p_{l,t}$$; HWE genotype mix $U=(p^2, 2pq, q^2)$
- $Z_{g,l,t}$: population fraction of genotype $g$ at $l,t$
- $U_{g,c} \in[0,1]$: susceptibility (prob. death) of genotype $g$ under insecticide class $c$
- Average population susceptibility:

$$ Q_{l,t,c}^{*}=\sum_{g} U_{g,c} Z_{g,l,t} $$

> Let $B$ be the number of loci, $G = 3^{B}$ the number of multilocus genotypes retained ($SS$, $SR$, $RR$ per locus; RS duplicates removed). Index genotypes by $g ∈ {1,…,G}$, loci by $lo ∈ {1,…,B}$, time by $t$.
# Road Map
## Defining the mathematical process
The mathematic path describing the model is in https://github.com/MiharisoaSylviane/phd_ir_mapping/tree/main/Model_description, IR_draft.Rmd
> Next step: Still need to define the math to determine the phenotype for each genotype $p_z$ 

## Simulation code matching the math by using others packages 
Checking if the code is producing the expected results
#### Minimal reproducible code
````
prob_left <- L*(1-P_t) + (1-L)*P_t
prob_right <- R*(1-P_t) + (1-R)*P_t
F <- prob_left * prob_right * (1 + L - R)
Z_prior <- apply(F, 1, prod); z_t <- Z_prior / sum(Z_prior)
````

````
#### Fitness (per locus):
SS <- (L==1) & (R==1); RR <- (L==0) & (R==0); SR <- (L==1) & (R==0)
G_loc <- 1*SS + w*RR + (h*w + (1-h))*SR
r_t <- apply(G_loc, 1, prod)

````
````
#### Update genotype distribution:
Z_star <- z_t * r_t; Z_next <- Z_star / sum(Z_star)
#### Update allele frequency per locus:
allele_from_Z <- function(Z_next, L, R) {
  G <- nrow(L); B <- ncol(L)
  stopifnot(length(Z_next) == G)
  a_counts <- (1 - L) + (1 - R)             # G x B ( resistant alleles)
  as.numeric(0.5 * (t(Z_next) %*% a_counts))   # length B
}


````
### Estimating the parameters used by running the model by fake real data
## Greta version of the code
Produce a greta version syntax of the simulation code

> Check: The goal is to assess the similarity of the results in a given parameters, assumptions and scenarios and data
## Fit greta version on the real data
### Simulation estimation study 
Minimal reproducting code


### Confirmation with real data



## Make prediction map



Sylviane Miharisoa · miharisoasylviane@gmail.com
Issues and ideas? please open a GitHub issue 
