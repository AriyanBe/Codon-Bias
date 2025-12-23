from Bio import SeqIO
import gzip
import os

# Input folder containing downloaded .gbff.gz files
input_folder = "Codons"

# Output FASTA file (all functional CDS)
output_file = os.path.join(input_folder, "Ecoli_functional_CDS.fasta")

count = 0

with open(output_file, "w") as out_f:
    for file in os.listdir(input_folder):
        if file.endswith(".gbff.gz"):
            print(f"Processing {file}...")
            gbff_path = os.path.join(input_folder, file)
            with gzip.open(gbff_path, "rt") as handle:
                for record in SeqIO.parse(handle, "genbank"):
                    for feature in record.features:
                        if feature.type == "CDS":
                            qualifiers = feature.qualifiers
                            product = qualifiers.get("product", [""])[0].lower()
                            gene = qualifiers.get("gene", ["unknown"])[0]

                            # Exclude hypothetical or pseudogenes
                            if "hypothetical" in product or "pseudogene" in product:
                                continue

                            if "translation" in qualifiers:
                                seq = qualifiers["translation"][0]
                                header = f">{gene}_{record.id}_{file}\n"
                                out_f.write(header)
                                out_f.write(seq + "\n")
                                count += 1

print(f"\nExtraction complete. {count} functional CDS saved to {output_file}")
