import pandas as pd
from Bio import SeqIO

# 1. Load the gene presence/absence matrix
df = pd.read_csv("gene_presence_absence.csv", low_memory=False)

# 2. Identify the isolate columns
# We exclude the known metadata columns to isolate just the sample columns
metadata_cols = ['Gene', 'Non-unique Gene name', 'Annotation']
isolate_columns = [col for col in df.columns if col not in metadata_cols]
n_samples = len(isolate_columns)  # Should be 33

# 3. Calculate 'No. isolates' (how many samples have a non-null entry for each gene)
# In Roary/Panaroo matrices, empty cells mean absence, filled cells mean presence.
# df[isolate_columns].notna().sum(axis=1) counts the non-empty cells across the sample columns.
df['No. isolates'] = df[isolate_columns].notna().sum(axis=1)

# 4. Define pangenome partitions based on the newly calculated column
core_genes = set(df[df['No. isolates'] >= 0.99 * n_samples]['Gene'])
soft_core_genes = set(df[(df['No. isolates'] >= 0.95 * n_samples) & (df['No. isolates'] < 0.99 * n_samples)]['Gene'])
shell_genes = set(df[(df['No. isolates'] >= 0.15 * n_samples) & (df['No. isolates'] < 0.95 * n_samples)]['Gene'])
cloud_genes = set(df[df['No. isolates'] < 0.15 * n_samples]['Gene'])

# 5. Parse the reference fasta and write to separate files
core_seqs, soft_core_seqs, shell_seqs, cloud_seqs = [], [], [], []

for record in SeqIO.parse("pan_genome_reference.fasta", "fasta"):
    gene_id = record.id.split()[0]
    if gene_id in core_genes:
        core_seqs.append(record)
    elif gene_id in soft_core_genes:
        soft_core_seqs.append(record)
    elif gene_id in shell_genes:
        shell_seqs.append(record)
    elif gene_id in cloud_genes:
        cloud_seqs.append(record)

# 6. Write the split FASTA files
SeqIO.write(core_seqs, "core_genes.fasta", "fasta")
SeqIO.write(soft_core_seqs, "soft_core_genes.fasta", "fasta")
SeqIO.write(shell_seqs, "shell_genes.fasta", "fasta")
SeqIO.write(cloud_seqs, "cloud_genes.fasta", "fasta")

print(f"Total genomes analyzed: {n_samples}")
print(f"Core: {len(core_seqs)} | Soft-Core: {len(soft_core_seqs)} | Shell: {len(shell_seqs)} | Cloud: {len(cloud_seqs)}")