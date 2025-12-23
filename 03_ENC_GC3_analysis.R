# =========================================================
# ENC–GC3 Analysis (FINAL FIXED VERSION)
# =========================================================

library(tidyverse)
library(seqinr)

# --- Paths ---
project_root <- "C:/Users/ariya/OneDrive - University Of Houston/Research Thesis"

fasta_file  <- file.path(project_root, "Codons", "Ecoli_functional_CDS.fasta")
results_dir <- file.path(project_root, "Results_RSCU Results")

dir.create(results_dir, showWarnings = FALSE)

output_csv <- file.path(results_dir, "ENC_GC3_results.csv")
output_png <- file.path(results_dir, "ENC_GC3_plot.png")

stopifnot(file.exists(fasta_file))

# --- Read CDS as DNA ---
cds <- read.fasta(fasta_file, seqtype = "DNA", as.string = TRUE)

# --- GC3 calculation ---
calc_GC3 <- function(seq) {
  seq <- toupper(seq)
  n <- nchar(seq)
  if (n < 3) return(NA)

  codons <- substring(seq, seq(1, n - 2, 3), seq(3, n, 3))
  third  <- substring(codons, 3, 3)

  mean(third %in% c("G", "C"))
}

# --- ENC calculation (seqinr-safe) ---
calc_ENC <- function(seq) {
  seq <- toupper(seq)
  n <- nchar(seq)
  if (n < 3) return(NA)

  # Extract codons safely
  codons <- substring(seq, seq(1, n - 2, 3), seq(3, n, 3))

  # Translate using s2c (CRITICAL FIX)
  aa <- translate(s2c(paste(codons, collapse = "")))

  df <- data.frame(codon = codons, aa = aa, stringsAsFactors = FALSE)

  # Remove stop codons
  df <- df[df$aa != "*", ]

  families <- split(df$codon, df$aa)

  F_vals <- sapply(families, function(fam) {
    k <- length(fam)
    if (k <= 1) return(NA)
    p <- table(fam) / k
    sum(p^2)
  })

  F_vals <- F_vals[!is.na(F_vals)]
  if (length(F_vals) == 0) return(NA)

  ENC <- 1 / mean(F_vals)
  min(ENC, 61)
}

# --- Compute per-gene ENC and GC3 ---
results <- tibble(
  Gene = names(cds),
  GC3  = sapply(cds, calc_GC3),
  ENC  = sapply(cds, calc_ENC)
) %>%
  filter(!is.na(ENC), GC3 > 0, GC3 < 1)

# --- Save results ---
write.csv(results, output_csv, row.names = FALSE)

# --- Expected ENC curve ---
gc3_vals <- seq(0.01, 0.99, by = 0.01)
enc_expected <- 2 + gc3_vals + 29 / (gc3_vals^2 + (1 - gc3_vals)^2)

# --- Plot ---
p <- ggplot(results, aes(GC3, ENC)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_line(
    data = data.frame(GC3 = gc3_vals, ENC = enc_expected),
    aes(GC3, ENC),
    color = "red",
    linewidth = 1
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "ENC–GC3 Plot of E. coli Coding Sequences",
    x = "GC3",
    y = "ENC"
  )

ggsave(output_png, p, width = 7, height = 6, dpi = 300)

cat("ENC–GC3 analysis complete\n")
