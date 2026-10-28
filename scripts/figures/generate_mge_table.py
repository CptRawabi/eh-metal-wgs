import pandas as pd
import glob
import os

# --- FILE PATHS ---
BACMET_FILE = "/mnt/d/Desktop/Project_Deinococcus/BacMet/intial data/bacmet_with_genomic_coordinates.tsv"
MEF_FILE = "/mnt/d/Desktop/Project_Deinococcus/MEF/MGE_Outputs/ZPb-04/ZPb-04.csv"
ISESCAN_PATTERN = "/mnt/d/Desktop/Project_Deinococcus/ISEScan/*tabular_results_*"
VIRSORTER_PATTERN = "/mnt/d/Desktop/Project_Deinococcus/Virsorter/*boundary_infor*"
OUTPUT_FILE = "/mnt/d/Desktop/Project_Deinococcus/ZPb-04_MGE_BacMet_Table.csv"

def check_overlap(start1, end1, start2, end2):
    return max(start1, start2) <= min(end1, end2)

# --- 1. LOAD BACMET DATA ---
print("Loading BacMet data...")
bacmet_df = pd.read_csv(BACMET_FILE, sep='\t')
bacmet_df = bacmet_df.dropna(subset=['Start', 'End'])
bacmet_df = bacmet_df.sort_values(by=['Contig', 'Start'])

mges = [] 

# --- 2. LOAD MGE DATA ---
print("Loading MGE data...")

# Load MEF
if os.path.exists(MEF_FILE):
    mef_df = pd.read_csv(MEF_FILE, comment='#')
    for _, row in mef_df.iterrows():
        clean_contig = str(row['contig']).split()[0] 
        mges.append({
            'contig': clean_contig, 
            'start': int(row['start']), 
            'end': int(row['end']), 
            'mge': f"{row['name']} ({row['type']})"
        })

# Load ISEscan
isescan_files = glob.glob(ISESCAN_PATTERN)
if isescan_files:
    try:
        ise_df = pd.read_csv(isescan_files[0], sep='\t')
        i_contig = [c for c in ise_df.columns if 'seq' in c.lower() or 'contig' in c.lower()][0]
        i_start = [c for c in ise_df.columns if 'begin' in c.lower() or 'start' in c.lower()][0]
        i_end = [c for c in ise_df.columns if 'end' in c.lower()][0]
        i_type = [c for c in ise_df.columns if 'family' in c.lower() or 'cluster' in c.lower()][0]
        for _, row in ise_df.iterrows():
            mges.append({'contig': str(row[i_contig]).split()[0], 'start': int(row[i_start]), 'end': int(row[i_end]), 'mge': f"{row[i_type]} (ISEscan)"})
    except Exception as e:
        print(f"Skipping ISEscan: {e}")

# Load VirSorter
virsorter_files = glob.glob(VIRSORTER_PATTERN)
if virsorter_files:
    try:
        vir_df = pd.read_csv(virsorter_files[0], sep='\t')
        v_contig = [c for c in vir_df.columns if 'seq' in c.lower() or 'contig' in c.lower() or 'name' in c.lower()][0]
        v_start = [c for c in vir_df.columns if 'start' in c.lower() or 'begin' in c.lower()][0]
        v_end = [c for c in vir_df.columns if 'end' in c.lower()][0]
        for _, row in vir_df.iterrows():
            mges.append({'contig': str(row[v_contig]).split('_')[0], 'start': int(row[v_start]), 'end': int(row[v_end]), 'mge': "Prophage (VirSorter)"})
    except Exception as e:
        print(f"Skipping VirSorter: {e}")

# --- 3. CLUSTER & MATCH ---
print("Clustering genes and finding overlaps...")
results = []
current_cluster = []
CLUSTER_MAX_DISTANCE = 10000

for _, row in bacmet_df.iterrows():
    contig = str(row['Contig']).strip()
    start = int(row['Start'])
    end = int(row['End'])
    
    raw_gene = str(row['BacMet_ID'])
    gene = raw_gene.split('|')[1] if '|' in raw_gene else raw_gene

    if not current_cluster:
        current_cluster = [{'contig': contig, 'start': start, 'end': end, 'gene': gene}]
    else:
        last_gene = current_cluster[-1]
        if contig == last_gene['contig'] and (start - last_gene['end']) <= CLUSTER_MAX_DISTANCE:
            current_cluster.append({'contig': contig, 'start': start, 'end': end, 'gene': gene})
        else:
            c_contig = current_cluster[0]['contig']
            c_start = min(g['start'] for g in current_cluster)
            c_end = max(g['end'] for g in current_cluster)
            c_genes = ", ".join(g['gene'] for g in current_cluster)
            
            overlap_mges = [m['mge'] for m in mges if (c_contig in m['contig'] or m['contig'] in c_contig) and check_overlap(c_start, c_end, m['start'], m['end'])]
            overlap_str = ", ".join(set(overlap_mges)) if overlap_mges else "None detected"
            
            results.append({'Contig': c_contig, 'Region': f"{c_start}-{c_end}", 'BacMet genes': c_genes, 'MGE overlap': overlap_str})
            current_cluster = [{'contig': contig, 'start': start, 'end': end, 'gene': gene}]

if current_cluster:
    c_contig = current_cluster[0]['contig']
    c_start = min(g['start'] for g in current_cluster)
    c_end = max(g['end'] for g in current_cluster)
    c_genes = ", ".join(g['gene'] for g in current_cluster)
    overlap_mges = [m['mge'] for m in mges if (c_contig in m['contig'] or m['contig'] in c_contig) and check_overlap(c_start, c_end, m['start'], m['end'])]
    overlap_str = ", ".join(set(overlap_mges)) if overlap_mges else "None detected"
    results.append({'Contig': c_contig, 'Region': f"{c_start}-{c_end}", 'BacMet genes': c_genes, 'MGE overlap': overlap_str})

# --- 4. EXPORT ---
final_df = pd.DataFrame(results)
final_df.to_csv(OUTPUT_FILE, index=False)
print(f"Success! Table saved to: {OUTPUT_FILE}")