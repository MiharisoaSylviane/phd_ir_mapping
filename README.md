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
 \(l\): location; \(t\): time; \(c\): insecticide class
- Bioassay: \(N\) tested, \(D\) died
- Genotypes (multi-locus): for record \(z\), \(M^{z}\) tested, \(N^{z}\) counts across genotype categories.
- Genotypes (single locus): allele-1 freq \(p_{l,t}\) (so \(q_{l,t}=1-p_{l,t}\)); HWE genotype mix \(U=(p^2,2pq,q^2)\)
- \(Z_{g,l,t}\): population fraction of genotype \(g\) at \((l,t)\)
- \(U_{g,c}\in[0,1]\): susceptibility (prob. death) of genotype \(g\) under insecticide class \(c\)
- Population susceptibility: $$\(Q^*_{l,t,c}=\sum_{g} U_{g,c}\, Z_{g,l,t}\)$$

# Building models
The mathematic path describing the model is in https://github.com/MiharisoaSylviane/phd_ir_mapping/tree/main/Model_description, IR_draft.Rmd
> Next step: Still need to define the mathematical process for the average multilocus genotype frequency Z_{t+1} at the next time point,and to determine the polygenic fitness r, which will be use the compute the next allele frequency at t+1. And also to define the phenotype resistant P* and Q* in function of insecticide resistance by using genotype
## Road Map
This road map contains the process donein the coding part of this model, we summarized the process explaining the why of each step:
### Step 1 : Per-locus genotype from allele frequency
For each locus lo and time/location, we define a prior p for the allele frequency of resistant Anopheles p_lo, initially from a prior uniform(0.1)
Map to genotype fractions under HWE: (p_lo^2,"  " 2p_lo (1-p_lo),"  "(1-p_lo )^2)
> *First draft in progress and only focus on the locus and genotype, time and insecticide class to consider after*
### Step 2 : Per locus genotype frequency composition
	Construct the multilocus genotype space g across loci
	Start with locus independence (linkage equilibrium); later, add LD/haplotypes if needed
 > First draft in progress
	Combine per-locus maps to get multilocus U_(g,c)
### Step 3 : Per-locus fitness, dominance and the polygenic fitness
Specify per-locus fitness (or susceptibility) parameters: w_lo and dominance h_lo for each insecticide class
From the silocus fitness, and dominance h, we define the polygenic fitness Gw (mention equation number)
>  * First draft in progress and only focus on the locus and genotype, time and insecticide class to consider after*, we have to mention in each step, line of code working but need to find a way to get Z
### Step 4 : Mutlticlocus genotype frequency for the next time point
	Infer the population mixture Z_(g,lo,t)using genotype data (Dirichlet–Multinomial)
	When only allele data are available, infer per-locus p_land back out implied genotype fractions
 > 
### Step 5 : Allele-frequency update 
	From Z_(g,lo,t)(and selection via U_(g,c)), update p_(lo,t+1)using standard selection recursions
	Optionally include migration or temporal smoothing to share information across space–time
### Step 6 : Population susceptibility / phenotype aggregation
	Compute population-level susceptibility (phenotype) for each insecticide class
>
## Simulation
### Step 7: Simulating data to test the model
	Creating simulated data to test based on the likelihood described in the mathematic description of the model
 >
## Application on real data
### Step 8: Visualisation and cleaning
	Normalized *Anopheles* data (detecting missing data, errors and others)
	Visualize the covariates
 >
### Step 9 : Checking and validation
	Assess the difference between the expected and the observed data for each step (genotype, allele frequency phenotype), by the distribution of residuals
 >
### Step 10 : Mapping and outputs
	- Generate maps and time series of p_lo, genotype fractions, and Q^* with uncertainty
	- discriminating concentration bioassay results against each class of insecticide
 >



Sylviane Miharisoa · miharisoasylviane@gmail.com
Issues and ideas? please open a GitHub issue 
