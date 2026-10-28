import pandas as pd
from Bio import SeqIO

# 1. Load CSV
df = pd.read_csv("gene_presence_absence.csv", low_memory=False)
csv_genes = set(df['Gene'].astype(str))

# 2. Parse FASTA
fasta_genes = set()
fasta_records = list(SeqIO.parse("pan_genome_reference.fasta", "fasta"))

for record in fasta_records:
    # Get the ID exactly as it appears before the first space
    fasta_genes.add(record.id)

print(f"Total genes in CSV: {len(csv_genes)}")
print(f"Total sequences in FASTA: {len(fasta_records)}")

# 3. Find the missing ones
missing_in_fasta = csv_genes - fasta_genes
print(f"Genes in CSV but NOT matched in FASTA: {len(missing_in_fasta)}")

if len(missing_in_fasta) > 0:
    print("\nExamples of missing genes (first 10):")
    print(list(missing_in_fasta)[:10])