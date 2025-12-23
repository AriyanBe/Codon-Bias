library(tidyverse)

project_root <- "C:/Users/ariya/OneDrive - University Of Houston/Codon Bias"
results_dir  <- file.path(project_root, "Results_RSCU Results")

input_csv  <- file.path(results_dir, "Ecoli_RSCU.csv")
output_csv <- file.path(results_dir, "ChiSquare_CodonBias_by_AA.csv")

stopifnot(file.exists(input_csv))

# Read RSCU data
RSCU <- read.csv(input_csv)

# ---- FIX COLUMN NAME ----
# If your column is named differently, rename it here
# (adjust if needed after checking colnames(RSCU))
RSCU <- RSCU %>% 
  rename(AminoAcid = AminoAcid)

# Remove stop codons
RSCU <- RSCU %>% filter(AminoAcid != "*")

# ---- REMOVE SINGLE-CODON AMINO ACIDS BEFORE TESTING ----
RSCU_multi <- RSCU %>%
  group_by(AminoAcid) %>%
  filter(n() > 1) %>%
  ungroup()

# ---- Chi-square analysis ----
chi_results <- RSCU_multi %>%
  group_by(AminoAcid) %>%
  group_modify(~ {
    observed <- .x$Count
    expected <- rep(sum(observed) / length(observed), length(observed))

    chisq <- suppressWarnings(
      chisq.test(observed, p = rep(1 / length(observed), length(observed)))
    )

    tibble(
      Codon = .x$Codon,
      AminoAcid = unique(.x$AminoAcid),
      Observed = observed,
      Expected = expected,
      StdResidual = (observed - expected) / sqrt(expected),
      ChiSquare_p = chisq$p.value
    )
  }) %>%
  ungroup()

write.csv(chi_results, output_csv, row.names = FALSE)

cat("✅ Chi-square codon bias analysis complete\n")
