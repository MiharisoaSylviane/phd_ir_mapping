# Title

# Introduction

Malaria remains one of the most pressing public health challenges
globally, particularly in sub-Saharan Africa. In 2022, approximately 94%
of global malaria cases and 95% of malaria-related deaths occurred in
this region, highlighting the persistent burden of this disease (WHO,
2023). Insecticide-treated nets (ITNs) and indoor residual spraying
(IRS) have been central to malaria control across the region,
consistently reducing human–vector contact and suppressing transmission
(Hobbs *et al.*, 2023). ITNs have been reported to contribute to
reductions of up to 50% in malaria incidence and 20% in malaria
mortality (Zikayo et *al.*, 2025), while IRS has been shown to reduce
malaria incidence by approximately 50% (Harrison, 2025). Regarding the
mode of application, ITNs are long-acting control tools that provide
protection for up to three years, whereas IRS is typically applied
annually (Willson *et al.*, 2015). Over the past two decades, ITNs have
become the dominant malaria prevention tool due to their scalability,
cost-effectiveness, and favorable safety profile (Batt *et al.*, 2015).
During most of this period, pyrethroids were the only WHO-recommended
active ingredient for ITNs and the main insecticide class used in IRS
(WHO, 2020), a constraint that simplified procurement and deployment but
imposed strong selective pressure on a single insecticide class (Moyes
*et al.*, 2021). Notably, pyrethroids have also demonstrated the highest
overall performance in malaria control compared to other insecticide
classes (Zhou *et al.*, 2022).

However, as ITN coverage expanded and usage intensified, resistance to
pyrethroids emerged rapidly and unevenly across Africa, generating
pronounced spatial and temporal heterogeneity in susceptibility among
major Anopheles vectors (Hancock *et al.*, 2020). Contributing to this
heterogeneity, resistance in African *Anopheles* arises through
multiple, often co-occurring mechanisms (Kinsinza *et al.*, 2017).
Target-site mutations in the voltage-gated sodium channel (Vgsc),
commonly referred to as knockdown resistance (kdr), reduce pyrethroid
binding affinity and impair knockdown efficacy. In parallel, metabolic
resistance—mediated by elevated cytochrome P450 monooxygenases,
glutathione S-transferases, and esterases—enhances detoxification
capacity and survival following insecticide exposure (Khan *et al.*,
2020). The prevalence of these mechanisms varies by species, region, and
time, and their frequencies can shift rapidly under changing
intervention pressures (Hancock, 2020; Kinsinza *et al.*, 2017).

Continental-scale syntheses that integrate WHO bioassay results with
molecular surveillance data (e.g., Vgsc genotype frequencies) have
revealed species-specific regional trends and heterogeneous resistance
trajectories, emphasizing the need for resistance-informed strategies
rather than uniform, one-size-fits-all approaches (Hancock, 2020).
Indoor insecticide-based interventions, particularly ITNs and IRS, tend
to suppress endophilic and anthropophilic vectors such as *Anopheles
gambiae* and *An. coluzzii* more strongly than the more exophagic *An.
arabiensis*, leading to measurable shifts in species composition that
can undermine the effectiveness of exclusively indoor tools and
complicate routine entomological monitoring (Sinka *et al.*, 2016).

Because resistance alters both survival and behavioural responses after
insecticide contact, and interventions are designed to directly or
indirectly reduce the mosquito community, they are expected to reshape
local vector abundance and community composition. Modelling and field
studies indicate that ITNs can elicit divergent behavioural outcomes
depending on their deterrent, inert, or attractive properties. Deterrent
ITNs reduce mosquito entry and enhance personal protection but may
decrease the likelihood of mosquitoes contacting treated surfaces, thus
limiting their community-level impact. Conversely, inert or attractive
ITNs may increase mosquito entry and contact rates, enhancing mortality
and reducing transmission when population-level coverage is high.
Moreover, the study by Sinka demonstrated the effect of interventions
(IRS and ITNs) on the relative abundance of three African Anopheles
species: Anopheles funestus, Anopheles gambiae/coluzzii, and Anopheles
arabiensis. Their findings revealed a diverse impact of interventions,
with the relative abundance of Anopheles funestus and Anopheles gambiae
tending to decrease, while that of Anopheles arabiensis tended to
increase.

## PROBLEM STATEMENT

Data on the true relative abundance of mosquito populations in the wild
remain limited and do not cover all areas where such information is
critical. The impact of interventions, such as insecticide-treated nets
(ITNs) and indoor residual spraying (IRS), on the relative abundance and
species composition of mosquito populations—comprising both resistant
and susceptible individuals—remains poorly understood. This gap has
prompted the use of modeling approaches to gain insights into these
dynamics. Current models of *Anopheles* relative abundance remain
limited. For example, Sinka *et al.*, 2010 modeled the relative
abundance of primary African vectors before and after indoor insecticide
interventions, and Killeen *et al.,* 2000 developed a simulation model
of African Anopheles ecology for malaria transmission analysis. While
these studies advanced understanding of vector population dynamics, they
often assume homogeneous susceptibility within species and do not fully
integrate insecticide resistance, spatiotemporal variation in
intervention coverage, or resulting shifts in species composition.
Consequently, such models may underestimate residual transmission and
fail to accurately identify areas requiring supplementary or alternative
vector control strategies, limiting their utility for evidence-based
planning. Despite growing evidence of widespread insecticide resistance
among African malaria vectors, the development of new insecticidal
compounds is neither a sustainable nor an immediate solution. The
discovery, formulation, and regulatory validation of novel insecticides
are time-consuming and costly, often requiring several years of
intensive research. Consequently, robust spatiotemporal mapping of
insecticide resistance is essential to guide decision-making and
identify geographic areas where the deployment of new products would
yield the greatest impact. Although insecticide resistance in African
malaria vectors is well documented and various modeling frameworks have
been developed, resistance mapping remains uncertain. Existing
approaches often have significant limitations: - They focus on limited
regions or species, failing to capture broader geographic and taxonomic
heterogeneity. - Current phenotypic resistance maps do not adequately
reflect resistance within species complexes and do not capture the
mechanisms under the resistance - Models based solely on allele
frequencies oversimplify resistance evolution, ignoring key factors such
as dominance effects, fitness trade-offs, and interactions among
multiple loci that together shape the phenotypic expression and spread
of resistance. Previous studies, such as those analyzing Vgsc-995S and
Vgsc-995F mutations, focused on a single target-site gene, revealing
geographically structured patterns and the influence of ITN coverage on
allele frequencies (Hancock, 2022). However, resistance is often
polygenic, involving multiple target-site and metabolic mechanisms that
interact to influence species composition and relative abundance under
realistic intervention pressures. Consequently, a multilocus genotype
approach, which accounts for allele combinations arising from
recombination, is essential to capture dominance patterns, resistance
expression, and to predict the speed and extent of resistance spread
under different intervention pressures. The effect of intervention on
the relative abundance and species composition in a population composed
by resistant and susceptible are not well understood.

## Research questions

To address gaps in understanding the interaction between resistance
evolution and vector control, this study explores three key questions: -
How has insecticide resistance been modeled over time? - How do
interventions and resistance influence changes in *Anopheles*
abundance? - How does resistance vary over time among different
Anopheles species and regions in Africa?

## Research objectives

-   Main objective The main objective of this study is to estimate and
    develop a spatiotemporal trends in insecticide resistance and the
    relative abundance of Anopheles species after intervention.

The specific objectives are: O1: Review of current insecticide
resistance models: To review existing models of insecticide resistance,
including the evolutionary processes considered, types of data used, and
model structures 02: Spatiotemporal model of the relative abundance of
Anopheles species: To investigate how insecticide resistance affects
changes in *Anopheles* relative abundance following ITN and IRS
deployment O3: Spatiotemporal model of insecticide resistance in
*Anopheles*: To integrate phenotypic (bioassay mortality) and genotypic
(multi-locus allele frequencies such as kdr and ace-1) data to
characterize spatiotemporal trends in resistance

## Scope of the study

Sustaining the effectiveness of malaria vector control interventions in
Africa requires a deep understanding of the genetic basis and dynamics
of insecticide resistance. Current strategies are constrained by limited
knowledge of how multi-locus resistance mechanisms interact and
influence Anopheles populations under real-world coverage of
insecticide-treated nets (ITNs) and indoor residual spraying (IRS).
Integrating entomological surveillance with phenotypic and genotypic
resistance data in a spatiotemporal modeling framework allows us to
capture patterns of species replacement, monitor the evolution of
resistance, and identify areas where interventions may be less
effective. This approach provides essential insights for evidence-based
decision-making, enabling the optimization and targeted deployment of
malaria control measures. Optimal deployment of interventions depends on
local resistance profiles, coverage histories, and vector ecology.
Without an explicit framework linking genotype-informed resistance to
vector abundance responses under realistic deployment scenarios, it
remains difficult to: - Prioritize products for deployment. - Schedule
rotations of insecticides effectively. - Forecast downstream effects on
malaria transmission.

## Literature review

### 1- Insecticide resistance mechanisms

Insecticide resistance in malaria vectors arises through multiple
adaptive mechanisms that enable mosquitoes to survive exposure to lethal
insecticide doses. These mechanisms are broadly classified into four
categories: target-site, metabolic, behavioral, and cuticular
resistance. Specifically, target-site resistance results from point
mutations in genes encoding insecticide-binding proteins, such as the
kdr mutation in voltage-gated sodium channels or the ace-1 mutation in
acetylcholinesterase, which reduce insecticide efficacy. And metabolic
resistance occurs through overexpression or amplification of detoxifying
enzymes—including cytochrome P450 monooxygenases, esterases, and
glutathione S-transferases—allowing rapid degradation of insecticides
and often conferring cross-resistance to multiple classes. Behavioral
resistance manifests as altered mosquito activity patterns, such as
feeding or resting outdoors, or avoiding oviposition on treated
surfaces, providing low to moderate levels of resistance and potentially
delaying the onset of physiological resistance. Cuticular resistance
arises from thickening or chemical modification of the exoskeleton,
which reduces insecticide penetration and offers moderate, sometimes
interclass, protection. These mechanisms frequently act in combination,
resulting in complex resistance phenotypes that compromise the efficacy
of standard vector control interventions and underscore the importance
of integrated resistance management strategies.

``` r
library(readxl)
data <- read_excel("C:/Users/Sylviane/Downloads/ir_mechanisms.xlsx")
# Display the data frame
knitr::kable(data, caption = "Data from CSV File")
```

| Mechanism | Key genes/families | Insecticide classes affected |
|:----------------------------------|:-----------------------|:-------------|
| Target-site, protein/target : sodium channel, AChE enzyme, GABA receptor | vgsc (voltage-gated sodium channel), ace-1 (acetylcholinesterase), rdl (GABA-gated chloride channel) | Pyrethroids/DDT, Organophosphates/carbamates, cyclodienes (dieldrin) |
| Metabolic, role in detoxicification : oxidation / breakdown of insecticides, conjugation of insecticide or metabolite to glutathione, hydrolysis or sequestration of insecticide molecules ) | Cytochrome P450 monooxygenases (CYPs (CYP6P3, CYP6M2)), Glutathione S- Transferases (GSTs (GSTe2)), Carboxylesterases, esterases | Multiple classes, cross-resistance likely |
| Cuticular / penetration (reducing insecticide uptake) | Cuticular protein genes, chitin-binding proteins | Slows uptake → broad effect |
| Behavioural/regulatory (enabling or facilitating behavioural avoidance) | Serine proteases, kinases, transcription factors | NA |

Data from CSV File

\<\<how insecticide resistance appear?\>\> what they mean by the
Molecular characterizations have revealed that various mutations in the
S1-S6 transmembrane segments of domain II of the sodium ion channel gene
give rise to resistance to these insecticides.

DDT and pyerthroids share the same mode of action on the insects nervous
system targeting the neuronal voltage gated channel sodium ion channel.

### 2- Relative abundance and species composition

### 3- Measurement of Insecticide resistance

### 4- Modeling insecticide resistance

### 5- Modeling the effect of intervention on the relative abundance

## 8- Research Methodology

The study will be divided into 3 chapters that will respond each
research question:

## 9- Research plan

## 10- schedule activities

-   research budget

## 11- References

Hobbs NP, Weetman D, Hastings IM. Insecticide resistance management
strategies for public health control of mosquitoes exhibiting polygenic
resistance: A comparison of sequences, rotations, and mixtures. Evol
Appl. 2023 Apr 5;16(4):936-959. doi: 10.1111/eva.13546. PMID: 37124088;
PMCID: PMC10130562.

Bhatt, S., Weiss, D., Cameron, E. et al. The effect of malaria control
on Plasmodium falciparum in Africa between 2000 and 2015. Nature 526,
207–211 (2015). <https://doi.org/10.1038/nature15535>

Moyes, C.L., Lees, R.S., Yunta, C. et al. Assessing cross-resistance
within the pyrethroids in terms of their interactions with key
cytochrome P450 enzymes and resistance in vector populations. Parasites
Vectors 14, 115 (2021). <https://doi.org/10.1186/s13071-021-04609-5>
Khan S, Uddin MN, Rizwan M, et al. Mechanism of Insecticide Resistance
in Insects/Pests. Polish Journal of Environmental Studies.
2020;29(3):2023–2030. <doi:10.15244/pjoes/108513>.

Kisinza, W.N., Nkya, T.E., Kabula, B. et al. Multiple insecticide
resistance in Anopheles gambiae from Tanzania: a major concern for
malaria vector control. Malar J 16, 439 (2017).
<https://doi.org/10.1186/s12936-017-2087-2>

Nwane, P., Etang, J., Chouaїbou, M. et al. Multiple insecticide
resistance mechanisms in Anopheles gambiae s.l. populations from
Cameroon, Central Africa. Parasites Vectors 6, 41 (2013).
<https://doi.org/10.1186/1756-3305-6-41> \<\<\< Declaration \>\>\>

Zikayo Amulaga R. (2025). A Review of the Efficacy of
Insecticide-Treated Bed Nets in Reducing Malaria Incidence among
Children under Five in Rural Sub-Saharan Africa. INOSR Scientific
Research 12(1):1-5. <https://doi.org/10.59298/INOSRSR/2025/12.1.1500>

Harrison Cole, 2025. stopmalarianow.org,
<https://stopmalarianow.org/evaluating-the-effectiveness-of-indoor-residual-spraying-in-malaria-endemic-regions/#what-is-indoor-residual-spraying-and-its-role-in-malaria-control>

Zhou Y, Zhang WX, Tembo E, Xie MZ, Zhang SS, Wang XR, Wei TT, Feng X,
Zhang YL, Du J, Liu YQ, Zhang X, Cui F, Lu QB. Effectiveness of indoor
residual spraying on malaria control: a systematic review and
meta-analysis. Infect Dis Poverty. 2022 Jul 23;11(1):83. doi:
10.1186/s40249-022-01005-8. PMID: 35870946; PMCID: PMC9308352

Depinay, JM.O., Mbogo, C.M., Killeen, G. et al. A simulation model of
African Anopheles ecology and population dynamics for the analysis of
malaria transmission. Malar J 3, 29 (2004).
<https://doi.org/10.1186/1475-2875-3-29>

Sinka, M.E., Golding, N., Massey, N.C. et al. Modelling the relative
abundance of the primary African vectors of malaria before and after the
implementation of indoor, insecticide-based vector control. Malar J 15,
142 (2016). <https://doi.org/10.1186/s12936-016-1187-8>

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](research_proposal_files/figure-markdown_github/pressure-1.png)

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
