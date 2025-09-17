# phd_ir_mapping
## Description
This repository maps insecticide resistance in Anopheles mosquitoes across Africa. It combines genotype data, phenotype (bioassay) outcomes, and allele-frequency measurements in a Bayesian hierarchical, semi-mechanistic framework to estimate resistance and susceptibility through space and time with quantified uncertainty.
The single locus version of this model could be found in https://github.com/idem-lab/ir_cube designed (Golding).

## Data inputs
- WHO discriminating concentration results from **Vector Atlas**, incorporating **IR Mapper** and **WHO Malaria Threats Map**.
- **Genotypes**: Vector Atlas genotype frequency data (single-locus and multi-locus, where available) if avalaible.
- **Covariates** (per location–time):
  - ITN use (modelled), (~timeframe to identify) (MAP).
  - IRS coverage, (~timeframe to identify) (WHO/MAP).
  - Population density, (~timeframe to identify)(WorldPop/MAP; UN-adjusted 2021–2030)
  - pesticide use (~timeframe and source to identify)
 
> Exact file paths and schemas will be documented in `configs/` as the code stabilizes.

## Notation
We will develop and add notation through times of coding
 \(l\): location; \(t\): time; \(c\): insecticide class.
- Bioassay: \(N\) tested, \(D\) died.
- Genotypes (multi-locus): for record \(z\), \(M^{z}\) tested, \(N^{z}\) counts across genotype categories.
- Genotypes (single locus): allele-1 freq \(p_{l,t}\) (so \(q_{l,t}=1-p_{l,t}\)); HWE genotype mix \(U=(p^2,2pq,q^2)\).
- \(Z_{g,l,t}\): population fraction of genotype \(g\) at \((l,t)\).
- \(U_{g,c}\in[0,1]\): susceptibility (prob. death) of genotype \(g\) under insecticide class \(c\).
- Population susceptibility: $$\(Q^*_{l,t,c}=\sum_{g} U_{g,c}\, Z_{g,l,t}\)$$

## Models
The mathematic path describing the model is in https://github.com/MiharisoaSylviane/phd_ir_mapping/tree/main/Model_description, IR_draft.Rmd
However to summarize Roadmap of this model:

## Road Map
### Step 1 : Per-locus genotype from allele frequency
	For each locus lo and time/location, specify allele frequency p_lo  initially from data or a prior).
	Map to genotype fractions under HWE: (p_lo^2,"  " 2p_lo (1-p_lo),"  "(1-p_lo )^2).
> First draft
### Step 2 : Single-locus fitness / phenotype mapping
	Specify per-locus fitness (or susceptibility) parameters: w_lo and dominance h_lo for each insecticide class.
	Derive per-genotype susceptibility U_(g,c)at a single locus and validate against bioassay data.
### Step 3 : Multilocus genotype composition
	Construct the multilocus genotype space g across loci.
	Start with locus independence (linkage equilibrium); later, add LD/haplotypes if needed.
	Combine per-locus maps to get multilocus U_(g,c)(semi-mechanistic linkage from genotype to phenotype).
### Step 4 : Mutlticlocus genotyope for the next time point
	Infer the population mixture Z_(g,lo,t)using genotype data (Dirichlet–Multinomial).
	When only allele data are available, infer per-locus p_land back out implied genotype fractions.
### Step 5 : Allele-frequency update 
	From Z_(g,lo,t)(and selection via U_(g,c)), update p_(lo,t+1)using standard selection recursions.
	Optionally include migration or temporal smoothing to share information across space–time.
### Step 6 : Population susceptibility / phenotype aggregation
	Compute population-level susceptibility (phenotype) for each insecticide class:
$$Q_(lo,t,c)^*=∑_g U_(g,c) " " Z_(g,lo,tStep 1 : Per-locus genotype from allele frequency
	For each locus lo and time/location, specify allele frequency p_lo  initially from data or a prior).
	Map to genotype fractions under HWE: (p_lo^2,"  " 2p_lo (1-p_lo),"  "(1-p_lo )^2).
>First draft
Step 2 : Single-locus fitness / phenotype mapping
	Specify per-locus fitness (or susceptibility) parameters: w_loand dominance h_lofor each insecticide class.
	Derive per-genotype susceptibility U_(g,c)at a single locus and validate against bioassay data.
Step 3 : Multilocus genotype composition
	Construct the multilocus genotype space g across loci.
	Start with locus independence (linkage equilibrium); later, add LD/haplotypes if needed.
	Combine per-locus maps to get multilocus U_(g,c)(semi-mechanistic linkage from genotype to phenotype).
Step 4 : Genotype mixture inference
	Infer the population mixture Z_(g,lo,t)using genotype data (Dirichlet–Multinomial).
	When only allele data are available, infer per-locus p_land back out implied genotype fractions.
Step 5 : Allele-frequency update (dynamics)
	From Z_(g,lo,t)(and selection via U_(g,c)), update p_(lo,t+1)using standard selection recursions.
	Optionally include migration or temporal smoothing to share information across space–time.
Step 6 : Population susceptibility / phenotype aggregation
	Compute population-level susceptibility (phenotype) for each insecticide class:
$$Q_(lo,t,c)^*=∑_g U_(g,c) " " Z_(g,lo,t)$$
	Report both per-locus contributions and the final multilocus Q^*.
Step 7: Simulating data to test the model
	Creating simulated data to test based on the likelihood described in the mathematic description of the model
Step 8: Visualisation and cleaning
	Normalized Anopheles data (detecting missing data, errors and others)
	Visualize the covariates
Step 9 : Checking & validation
	Assess the difference between the expected and the observed data for each step (genotype, allele frequency phenotype)
Step 10 : Mapping and outputs
	Generate maps and time series of p_lo, genotype fractions, and Q^*with uncertainty.
	discriminating concentration bioassay results against each class of insecticide

	Report both per-locus contributions and the final multilocus Q^*.
### Step 7: Simulating data to test the model
	Creating simulated data to test based on the likelihood described in the mathematic description of the model
### Step 8: Visualisation and cleaning
	Normalized Anopheles data (detecting missing data, errors and others)
	Visualize the covariates
### Step 9 : Checking & validation
	Assess the difference between the expected and the observed data for each step (genotype, allele frequency phenotype)
### Step 10 : Mapping and outputs
	- Generate maps and time series of p_lo, genotype fractions, and Q^*with uncertainty.
	- discriminating concentration bioassay results against each class of insecticide

