# =========================================================
# Goodness-of-Fit Visualization for Codon Usage Bias
# =========================================================

library(tidyverse)

# --- Paths ---
project_root <- "C:/Users/ariya/OneDrive - University Of Houston/Codon Bias"
results_dir  <- file.path(project_root, "Results_RSCU Results")

input_csv <- file.path(results_dir, "ChiSquare_CodonBias_by_AA.csv")
stopifnot(file.exists(input_csv))

# --- Read chi-square results ---
chi_results <- read.csv(input_csv)

# --- Explicitly remove stop codons (clarity + safety) ---
chi_results <- chi_results %>%
  filter(AminoAcid != "*")

# =========================================================
# SUMMARY TABLE (Main Results Table)
# =========================================================

summary_table <- chi_results %>%
  group_by(AminoAcid) %>%
  summarise(
    n_codons = n(),
    p_value = first(ChiSquare_p),
    preferred_codon = Codon[which.max(StdResidual)],
    max_residual = max(StdResidual),
    .groups = "drop"
  ) %>%
  arrange(p_value)

write.csv(
  summary_table,
  file.path(results_dir, "ChiSquare_summary_table.csv"),
  row.names = FALSE
)

# =========================================================
# RESCALE RESIDUALS (FOR VISUALIZATION ONLY)
# =========================================================
# NOTE: This does NOT affect statistics, only plotting

max_abs_resid <- max(abs(chi_results$StdResidual), na.rm = TRUE)

chi_results <- chi_results %>%
  mutate(StdResidual_scaled = StdResidual / max_abs_resid)

# =========================================================
# FACET PLOT (OVERVIEW FIGURE)
# =========================================================

facet_plot <- ggplot(chi_results, aes(Codon, StdResidual_scaled)) +
  geom_col(fill = "steelblue") +
  geom_hline(
    yintercept = c(-0.1, 0.1),
    linetype = "dashed",
    color = "red",
    linewidth = 0.6
  ) +
  facet_wrap(~ AminoAcid, scales = "free_x") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    title = "Standardized residuals of synonymous codons across amino acids",
    subtitle = "Residuals scaled for visualization; dashed lines indicate strong deviation",
    x = "Codon",
    y = "Scaled standardized residual"
  )

ggsave(
  filename = file.path(results_dir, "CodonBias_residuals_facet.png"),
  plot = facet_plot,
  width = 12,
  height = 8,
  dpi = 300
)

cat("GOF summary table and facet plot generated successfully\n")
