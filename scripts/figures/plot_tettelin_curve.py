import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit

# Configure Arial font and publication-quality typography
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# 1. Load Panaroo presence/absence matrix
df = pd.read_csv("gene_presence_absence.csv")

metadata_cols = [
    'Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 
    'Avg group size in seqs', 'Min group size in seqs', 'Max group size in seqs', 
    'Average sequence identity in seqs', 'Pan-genome reference', 'Core accessory', 
    'Conserved accessory', 'Cloud'
]

strain_cols = [c for c in df.columns if c not in metadata_cols]

matrix = df[strain_cols].notna().astype(int).values
n_strains = matrix.shape[1]
n_genes = matrix.shape[0]

print(f"Loaded {n_genes} gene clusters across {n_strains} strains.")

# 2. Permutation simulation for Core and Pan-genome
n_permutations = 100
num_genomes = np.arange(1, n_strains + 1)

core_curves = np.zeros((n_permutations, n_strains))
pan_curves = np.zeros((n_permutations, n_strains))

np.random.seed(42)
for i in range(n_permutations):
    perm_indices = np.random.permutation(n_strains)
    shuffled_matrix = matrix[:, perm_indices]

    core_sizes = []
    pan_sizes = []

    for n in range(1, n_strains + 1):
        subset = shuffled_matrix[:, :n]
        core_sizes.append(np.sum(np.all(subset == 1, axis=1)))
        pan_sizes.append(np.sum(np.any(subset == 1, axis=1)))

    core_curves[i, :] = core_sizes
    pan_curves[i, :] = pan_sizes

# Calculate average counts per genome step
mean_pan = np.mean(pan_curves, axis=0)
mean_core = np.mean(core_curves, axis=0)

# 3. Curve Fitting for Smooth Fitted Lines
def heaps_pan_fit(N, A, B):
    return A * (N ** B)

def exp_core_fit(N, A, B, C):
    return A * np.exp(-B * N) + C

# Smooth x-values for continuous fitted curves
x_smooth = np.linspace(1, n_strains, 300)

popt_pan, _ = curve_fit(heaps_pan_fit, num_genomes, mean_pan, p0=[mean_pan[0], 0.4])
print(f"Pan-genome fit: A = {popt_pan[0]:.2f}, B = {popt_pan[1]:.3f}")
pan_smooth_y = heaps_pan_fit(x_smooth, *popt_pan)

try:
    popt_core, _ = curve_fit(exp_core_fit, num_genomes, mean_core, p0=[mean_core[0], 0.3, mean_core[-1]])
    core_smooth_y = exp_core_fit(x_smooth, *popt_core)
except Exception:
    core_smooth_y = np.interp(x_smooth, num_genomes, mean_core)

# 4. Plotting (Scatter Points + Smooth Curves)
fig, ax = plt.subplots(figsize=(8, 5.5), dpi=300)

ax.set_facecolor('white')
fig.patch.set_facecolor('white')

# Bright high-contrast colors (matching reference style)
# Contrasting Light Green & Bright Orange colors
color_pan_line = '#FF7F00'      # Bright Orange (Line)
color_pan_points = '#FFB366'    # Soft Orange (Points)
color_core_line = '#2E7D32'     # Forest Green (Line)
color_core_points = '#81C784'   # Light Green (Points)

# 1. Plot Fitted Trend Lines
ax.plot(x_smooth, pan_smooth_y, color=color_pan_line, lw=1.8, label='Pan genome', zorder=2)
ax.plot(x_smooth, core_smooth_y, color=color_core_line, lw=1.8, label='Core genome', zorder=2)

# 2. Plot Discrete Data Points (Triangles)
ax.scatter(num_genomes, mean_pan, color=color_pan_points, marker='^', s=45, 
           edgecolor='none', label='Total gene families', zorder=3)
ax.scatter(num_genomes, mean_core, color=color_core_points, marker='^', s=45, 
           edgecolor='none', label='Core gene families', zorder=3)

# Axis labels & Limits
ax.set_xlabel("Number of Genomes", fontsize=12, fontfamily='Arial', labelpad=10)
ax.set_ylabel("Number of Gene Families", fontsize=12, fontfamily='Arial', labelpad=10)

ax.set_xlim(0.8, n_strains + 0.5)
ax.tick_params(axis='both', which='major', labelsize=10, direction='in', length=5, width=0.8)

# Full clean box frame around plot
for spine in ['top', 'bottom', 'left', 'right']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(0.8)
    ax.spines[spine].set_color('#333333')

# Legend placed outside to the right of the plot
legend = ax.legend(
    loc='center left', 
    bbox_to_anchor=(1.02, 0.5), 
    frameon=True, 
    facecolor='white', 
    edgecolor='#CCCCCC', 
    framealpha=1.0, 
    fontsize=9.5,
    borderpad=0.8,        # Increases padding inside the box
    handletextpad=0.8,    # Fixes text touching the line/marker handles
    labelspacing=0.6      # Adds consistent vertical space between rows
)

plt.tight_layout()
plt.savefig('tettelin_pangenome_curve.png', dpi=300, bbox_inches='tight')
print("Saved scatter + fitted curve plot to tettelin_pangenome_curve.png!")