import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

# 1. Configure Arial font and publication-quality typography
mpl.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# 2. Load the specific full pairwise Galaxy dataset
file_name = 'Galaxy40-[FastANI on dataset 1-33 Output].tabular'
df = pd.read_csv(file_name, sep='\t', header=None, engine='python')

# 3. Helper function to strip directory paths and clean strain names
def clean_strain_name(val):
    s = str(val).strip()
    s = os.path.basename(s)
    for suffix in ['.fna', '.fasta', '.fa', '_fna_query', '_fna_ref', '_query', '_ref']:
        s = s.replace(suffix, '')
    s = s.replace('ZPb_04', 'ZPb4').replace('ZPb-04', 'ZPb4')
    return s

# Clean query, reference, and ANI columns
df['Query'] = df.iloc[:, 0].apply(clean_strain_name)
df['Ref'] = df.iloc[:, 1].apply(clean_strain_name)
df['ANI'] = pd.to_numeric(df.iloc[:, 2], errors='coerce')

# 4. Pivot into NxN matrix
matrix = df.pivot(index='Query', columns='Ref', values='ANI')

# Symmetrize matrix and enforce self-comparison = 100%
matrix = matrix.combine_first(matrix.T)
matrix = (matrix + matrix.T) / 2.0
for strain in matrix.columns:
    matrix.loc[strain, strain] = 100.0

# Fill missing non-pairwise cells with global minimum to avoid NaNs
fill_val = matrix.min().min() if not pd.isna(matrix.min().min()) else 95.0
matrix = matrix.fillna(fill_val)

# 5. Compute shared linkage tree on an explicit writable NumPy array copy
dist_matrix = 100.0 - matrix
dist_array = dist_matrix.to_numpy(copy=True)
np.fill_diagonal(dist_array, 0.0)

# Ensure valid distance bounds
dist_array = np.maximum(dist_array, 0.0)
condensed_dist = squareform(dist_array, checks=False)
shared_linkage = linkage(condensed_dist, method='average')

# 6. Build Strictly Green Colormap
green_colors = ['#e0f2f1', '#a5d6a7', '#4caf50', '#2e7d32', '#0d3b1e']
strict_green_cmap = LinearSegmentedColormap.from_list("StrictGreens", green_colors)

vmin_val = max(95.0, float(matrix.min().min()))

# 7. Generate Clustermap
g = sns.clustermap(
    matrix,
    row_linkage=shared_linkage,
    col_linkage=shared_linkage,
    cmap=strict_green_cmap,
    vmin=vmin_val,
    vmax=100.0,
    dendrogram_ratio=(0.10, 0.10),
    linewidths=0.3,
    linecolor='white',
    figsize=(11, 10),
    cbar_pos=(0.96, 0.35, 0.025, 0.3)
)

# 8. Unbolded axis tick labels and formatting
g.ax_heatmap.set_xlabel('', fontweight='normal')
g.ax_heatmap.set_ylabel('', fontweight='normal')

plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=9, fontweight='normal')
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=9, fontweight='normal')

# Format Colorbar / Legend
g.cax.set_ylabel('Average Nucleotide Identity (%)', fontsize=10, fontweight='normal', labelpad=10)
g.cax.tick_params(labelsize=9)

# Save figure
plt.savefig('FastANI_Clustered_Heatmap_Final.png', dpi=300, bbox_inches='tight')
print("Heatmap generated and saved as 'FastANI_Clustered_Heatmap_Final.png'.")