# Create dummy matrices for Sylviane's multilocus diploid selection model

# To compute the change in allele frequencies across multiple loci in Sylviane's
# model, we need to create dummy matrices that encode the allele across each
# loci for all possible genotypes. We do this with two dummy matrices, one for
# the left hand allele and one for the right hand allele. In these matrices, a 1
# indicates the S allele at that locus and hand, and a 0 indicates the R allele.
# For heterozygotes, we assign the susceptible allele to the left and the
# resistant to the right, to reduce the number of genotypes that need to be
# considered in the calculations.

# # Example code:
# 
# # For the single locus case, these encode (1 = S, 0 = R):
# # SS
# # SR
# # RR
# create_dummy_matrices(n_loci = 1)
# 
# # For the two locus case, these encode (1 = S, 0 = R):
# # SS SS
# # SR SS
# # RR SS
# # SS SR
# # SR SR
# # RR SR
# # SS RR
# # SR RR
# # RR RR
# create_dummy_matrices(n_loci = 2)


create_dummy_matrices <- function(n_loci) {
  
  # Steps:
  # 1. Create all unique combinations of alleles across all loci.
  # 2. Remove heterozygotes with the resistant allele on the left (since these
  # are duplicates of the heterozygote genotypes with the resistant allele on
  # the right)
  # 3. Create dummy matrices for the left and right alleles
  
  
  # this function uses tidyr and dplyr code, so load tidyverse to access both
  require(tidyverse)
  
  # 1. Create all unique combinations of alleles across all loci.

  # The following code automates using expand_grid to create all possible
  # combinations of alleles (with duplicate heterozygotes) across the specified
  # number of loci
  
  # The actual code below is equivalent to the following example code (for
  # n_loci = 3), but enables passing the number of loci as a function argument,
  # rather than hard-coding it
  
  # possible_genotypes <- expand_grid(
  #   locus_3_left = c("susceptible", "resistant"),
  #   locus_3_right = c("susceptible", "resistant"),
  #   locus_2_left = c("susceptible", "resistant"),
  #   locus_2_right = c("susceptible", "resistant"),
  #   locus_1_left = c("susceptible", "resistant"),
  #   locus_1_right = c("susceptible", "resistant"),
  # )
  
  # create a list of these loci, left and right hand copies, and the two allele
  # types
  full_list <- list(
    left = c("susceptible", "resistant"),
    right = c("susceptible", "resistant")
  ) %>%
    replicate(
      n_loci,
      .,
      simplify = FALSE
    ) %>%
    unlist(
      recursive = FALSE
    )
  
  # set the names according to the locus and hand
  names(full_list) <- sprintf(
    "locus_%i_%s",
    rep(rev(seq_len(n_loci)), each = 2),
    rep(c("left", "right"), n_loci)
  )
  
  # apply expand_grid to get all possible combinations of these - the full set
  # of genotypes with duplicate heterozygotes
  possible_genotypes <- do.call(expand_grid,
                                full_list)
  
  # reformat this so that it is easier to work with, with numbered genotypes,
  # but still with duplicate heterozygotes
  genotypes_all <- possible_genotypes %>%
    mutate(
      genotype = row_number(),
      .before = everything()
    ) %>%
    pivot_longer(
      cols = starts_with("locus"),
      names_to = c("locus", "hand"),
      names_pattern = "(.*)_(.*)",
      values_to = "allele"
    ) %>%
    arrange(
      genotype,
      locus,
      hand
    )
  
  # 2. Remove heterozygotes with the resistant allele on the left (since these
  # are duplicates of the heterozygote genotypes with the resistant allele on
  # the right)
  
  # remove the duplicate heterozygotes
  genotypes <- genotypes_all %>%
    # assess each locus within each genotype
    group_by(
      genotype,
      locus
    ) %>%
    mutate(
      # label the heterozygotes
      heterozygote = n_distinct(allele) == 2,
      # label those with the left allele resistant
      left_resistant = hand == "left" & allele == "resistant",
      # label the RS heterozygote loci (redundant as we have the SR heterozygotes)
      redundant = heterozygote & left_resistant
    ) %>%
    # remove any *genotype* with a redundant heterozygote, since the whole
    # genotype is a duplicate
    group_by(
      genotype
    ) %>%
    filter(
      !any(redundant)
    ) %>%
    ungroup() %>%
    # remove the columns we used to compute these
    select(
      -heterozygote,
      -left_resistant,
      -redundant
    ) %>%
    # relabel the genotypes from 1:n_loci^3
    mutate(
      genotype = match(genotype, unique(genotype))
    )
  
  # 3. Create dummy matrices for the left and right alleles
  
  # now we have all the genotypes, compute the two matrices
  left_dummy <- genotypes %>%
    # only need the left hand alleles for this matrix
    filter(
      hand == "left"
    ) %>%
    # encode susceptibility as a dummy variable
    mutate(
      susceptible = case_when(
        allele == "susceptible" ~ 1,
        allele == "resistant" ~ 0,
      )
    ) %>%
    # remove the redundant columns
    select(
      -hand,
      -allele
    ) %>%
    # convert into a matrix
    pivot_wider(
      names_from = locus,
      values_from = susceptible
    ) %>%
    select(
      -genotype
    ) %>%
    as.matrix()
  
  right_dummy <- genotypes %>%
    # only need the right hand alleles for this matrix
    filter(
      hand == "right"
    ) %>%
    # encode susceptibility as a dummy variable
    mutate(
      susceptible = case_when(
        allele == "susceptible" ~ 1,
        allele == "resistant" ~ 0,
      )
    ) %>%
    # remove the redundant columns
    select(
      -hand,
      -allele
    ) %>%
    # convert into a matrix
    pivot_wider(
      names_from = locus,
      values_from = susceptible
    ) %>%
    select(
      -genotype
    ) %>%
    as.matrix()
  
  # return these in a list
  list(
    left = left_dummy,
    right = right_dummy
  )
  
}





