# --- Install and load required packages ---
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("seqinr")) install.packages("seqinr")
library(tidyverse)
library(seqinr)

# --- Define folder paths ---
base_dir <- "Tests"
codon_dir <- file.path(base_dir, "Codons")
results_dir <- file.path(base_dir, "RSCU Results")  # New folder for results
dir.create(results_dir, showWarnings = FALSE)       # Create folder if it doesn't exist

input_csv <- file.path(codon_dir, "Ecoli_codon_usage.csv")
output_csv <- file.path(results_dir, "Ecoli_RSCU.csv")
output_plot <- file.path(results_dir, "RSCU_plot.png")

# --- Read codon usage data ---
codon_data <- read.csv(input_csv)

# --- Define the standard genetic code ---
GENETIC_CODE <- c(
  "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L",
  "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L",
  "ATT"="I", "ATC"="I", "ATA"="I", "ATG"="M",
  "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V",
  "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S",
  "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P",
  "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T",
  "GCT"="A", "GCC"="A", "GCA"="A", "GCG"="A",
  "TAT"="Y", "TAC"="Y", "TAA"="*", "TAG"="*",
  "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q",
  "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K",
  "GAT"="D", "GAC"="D", "GAA"="E", "GAG"="E",
  "TGT"="C", "TGC"="C", "TGA"="*", "TGG"="W",
  "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R",
  "AGT"="S", "AGC"="S", "AGA"="R", "AGG"="R",
  "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G"
)

# --- Build codon table ---
codon_table <- data.frame(
  Codon = names(GENETIC_CODE),
  AminoAcid = unlist(GENETIC_CODE),
  stringsAsFactors = FALSE
)

# --- Merge codon counts with amino acids ---
merged <- merge(codon_data, codon_table, by = "Codon", all.x = TRUE)

# --- Compute RSCU values ---
RSCU <- merged %>%
  group_by(AminoAcid) %>%
  mutate(RSCU = Count / (mean(Count[!is.na(AminoAcid)]))) %>%
  ungroup()

# --- Save RSCU values ---
write.csv(RSCU, output_csv, row.names = FALSE)
print(paste0("RSCU values calculated and saved to ", output_csv))

# --- Plot RSCU ---
library(ggplot2)
p <- ggplot(RSCU, aes(x = Codon, y = RSCU, fill = AminoAcid)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Codon Usage Bias in E. coli", x = "Codon", y = "RSCU Value") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# --- Save plot ---
ggsave(filename = output_plot, plot = p, width = 12, height = 6, dpi = 300)
print(paste0("RSCU plot saved to ", output_plot))
