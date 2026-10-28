import glob
import re
import io
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

for file_path in glob.glob("*.gbk"):
    if file_path.startswith("compliant_"):
        continue
        
    with open(file_path, "r") as f:
        lines = f.readlines()
        
    # STEP 1: Fix broken LOCUS lines caused by Prokka name overflows
    for i, line in enumerate(lines):
        if line.startswith("LOCUS"):
            # Scan ahead to find the true contig length from the 'source' feature
            length = "1000" # fallback
            for j in range(i+1, min(i+100, len(lines))):
                match = re.search(r'source\s+1\.\.(\d+)', lines[j])
                if match:
                    length = match.group(1)
                    break
                    
            # Extract the date
            date_match = re.search(r'(\d{2}-[A-Z]{3}-\d{4})', line)
            date = date_match.group(1) if date_match else "01-JAN-2000"
            
            # Extract the name and forcefully truncate to GenBank's 16-character limit
            raw_name = line[12:].strip().split()[0] if len(line) > 12 else "contig"
            safe_name = raw_name[:16].ljust(16)
            
            # Reconstruct a strictly compliant LOCUS line
            lines[i] = f"LOCUS       {safe_name} {length.rjust(11)} bp    DNA     linear       {date}\n"
            
    fixed_content = "".join(lines)
    
    # STEP 2: Use Biopython to add the missing 'gene' parent features
    records = []
    try:
        # Feed the repaired text directly into Biopython
        for record in SeqIO.parse(io.StringIO(fixed_content), "genbank"):
            new_features = []
            for feat in record.features:
                # Insert a 'gene' feature right before every 'CDS'
                if feat.type == "CDS":
                    gene_qualifiers = {}
                    if "locus_tag" in feat.qualifiers:
                        gene_qualifiers["locus_tag"] = feat.qualifiers["locus_tag"]
                    if "gene" in feat.qualifiers:
                        gene_qualifiers["gene"] = feat.qualifiers["gene"]

                    gene_feat = SeqFeature(
                        location=feat.location,
                        type="gene",
                        qualifiers=gene_qualifiers
                    )
                    new_features.append(gene_feat)
                new_features.append(feat)
                
            record.features = new_features
            records.append(record)
            
        out_name = f"compliant_{file_path}"
        SeqIO.write(records, out_name, "genbank")
        print(f"Successfully converted: {file_path} -> {out_name}")
    except Exception as e:
        print(f"Failed to process {file_path}: {e}")