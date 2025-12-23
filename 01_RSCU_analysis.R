# =========================================================
# Codon Usage Bias (RSCU) Analysis – E. coli
# =========================================================

# --- Load packages ---
library(tidyverse)
library(seqinr)

# --- Define paths ---
project_root <- "C:/Users/ariya/OneDrive - University Of Houston/Research Thesis"

codon_dir   <- file.path(project_root, "Codons")
results_dir <- file.path(project_root, "Results_RSCU Results")

dir.create(results_dir, showWarnings = FALSE)

input_csv  <- file.path(codon_dir, "Ecoli_codon_usage.csv")
output_csv <- file.path(results_dir, "Ecoli_RSCU.csv")
output_png <- file.path(results_dir, "RSCU_plot.png")

# --- Safety check ---
stopifnot(file.exists(input_csv))

# --- Read codon counts ---
codon_data <- read.csv(input_csv)

# --- Standard DNA genetic code ---
GENETIC_CODE <- c(
  "TTT"="F","TTC"="F","TTA"="L","TTG"="L",
  "CTT"="L","CTC"="L","CTA"="L","CTG"="L",
  "ATT"="I","ATC"="I","ATA"="I","ATG"="M",
  "GTT"="V","GTC"="V","GTA"="V","GTG"="V",
  "TCT"="S","TCC"="S","TCA"="S","TCG"="S",
  "CCT"="P","CCC"="P","CCA"="P","CCG"="P",
  "ACT"="T","ACC"="T","ACA"="T","ACG"="T",
  "GCT"="A","GCC"="A","GCA"="A","GCG"="A",
  "TAT"="Y","TAC"="Y","TAA"="*","TAG"="*",
  "CAT"="H","CAC"="H","CAA"="Q","CAG"="Q",
  "AAT"="N","AAC"="N","AAA"="K","AAG"="K",
  "GAT"="D","GAC"="D","GAA"="E","GAG"="E",
  "TGT"="C","TGC"="C","TGA"="*","TGG"="W",
  "CGT"="R","CGC"="R","CGA"="R","CGG"="R",
  "AGT"="S","AGC"="S","AGA"="R","AGG"="R",
  "GGT"="G","GGC"="G","GGA"="G","GGG"="G"
)

# --- Build codon → amino acid table ---
codon_table <- data.frame(
  Codon = names(GENETIC_CODE),
  AminoAcid = unname(GENETIC_CODE),
  stringsAsFactors = FALSE
)

# --- Merge codon counts with amino acids ---
merged <- merge(codon_data, codon_table, by = "Codon")

# --- Compute RSCU ---
RSCU <- merged %>%
  group_by(AminoAcid) %>%
  mutate(RSCU = Count / mean(Count)) %>%
  ungroup()

# --- Remove stop codons (CRITICAL FIX) ---
RSCU <- RSCU %>% filter(AminoAcid != "*")

# --- Order codons by amino acid, then codon ---
RSCU$Codon <- factor(
  RSCU$Codon,
  levels = RSCU %>%
    arrange(AminoAcid, Codon) %>%
    pull(Codon)
)

# --- Save cleaned RSCU table ---
write.csv(RSCU, output_csv, row.names = FALSE)

# --- Publication-quality plot (clean & interpretable) ---
p <- ggplot(RSCU, aes(x = Codon, y = RSCU)) +
  geom_col(fill = "steelblue") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Codon Usage Bias in E. coli",
    x = "Codon",
    y = "RSCU"
  )

ggsave(output_png, p, width = 14, height = 6, dpi = 300)

cat("RSCU analysis complete. Results saved.\n")
