from Bio import SeqIO
from collections import Counter
import os

# Input FASTA (all functional CDS)
fasta_file = "Tests/Codons/Ecoli_functional_CDS.fasta"

# Output CSV file
output_csv = "Tests/Codons/Ecoli_codon_usage.csv"

# Initialize codon counter
codons = [a+b+c for a in "ATGC" for b in "ATGC" for c in "ATGC"]
total_codon_counts = Counter({codon: 0 for codon in codons})

# Count codons across all CDS
for record in SeqIO.parse(fasta_file, "fasta"):
    seq = str(record.seq).upper()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in total_codon_counts:
            total_codon_counts[codon] += 1

# Print codon usage counts
print("Codon Usage Counts:")
for codon, count in total_codon_counts.items():
    print(f"{codon}: {count}")

# Save to CSV in Tests/Codons/
os.makedirs(os.path.dirname(output_csv), exist_ok=True)
with open(output_csv, "w") as f:
    f.write("Codon,Count\n")
    for codon, count in total_codon_counts.items():
        f.write(f"{codon},{count}\n")

print(f"\nCodon usage data saved to {output_csv}")
