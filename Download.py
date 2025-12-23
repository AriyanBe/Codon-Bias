from Bio import Entrez
import urllib.request
import os

# Email for NCBI
Entrez.email = ""

# Create folder structure
base_dir = "."
codon_dir = os.path.join(base_dir, "Codons")

os.makedirs(codon_dir, exist_ok=True)

# Search for E. coli assemblies
handle = Entrez.esearch(
    db="assembly",
    term="Escherichia coli[Organism]",
    retmax=250  # Increase this if you want more genomes
)

record = Entrez.read(handle)
assembly_ids = record["IdList"]
print(f"Found {len(assembly_ids)} assemblies: {assembly_ids}")

# Download each assembly file
for asm_id in assembly_ids:
    handle = Entrez.esummary(db="assembly", id=asm_id, report="full")
    record = Entrez.read(handle)
    docsum = record['DocumentSummarySet']['DocumentSummary'][0]
    ftp_path = docsum["FtpPath_RefSeq"]

    if ftp_path:
        gbff_file = ftp_path + "/" + ftp_path.split("/")[-1] + "_genomic.gbff.gz"

        # Save into Tests/Codons/
        local_file = os.path.join(codon_dir, os.path.basename(gbff_file))

        print(f"Downloading: {gbff_file}")
        try:
            urllib.request.urlretrieve(gbff_file, local_file)
            print(f"Saved as: {local_file}\n")
        except Exception as e:
            print(f"Failed to download {gbff_file}: {e}\n")

print("Download complete.")
