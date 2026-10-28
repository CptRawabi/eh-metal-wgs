import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

# 1. Clean, Publication-Quality Settings
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = 1.0
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['ytick.major.width'] = 1.0
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False

def parse_eggnog_annotations(file_path):
    """Parses eggNOG-mapper tabular files and counts COG categories."""
    header = None
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#query'):
                header = line.strip().lstrip('#').split('\t')
                break
                
    df = pd.read_csv(file_path, sep='\t', comment='#', names=header)
    cog_col = 'COG_category' if 'COG_category' in df.columns else df.columns[6]
    
    cogs = df[cog_col].dropna()
    cogs = cogs[~cogs.isin(['-', ''])]
    split_cogs = cogs.apply(lambda x: list(str(x).replace(' ', ''))).explode()
    
    counts = split_cogs.value_counts().reset_index()
    counts.columns = ['COG_Letter', 'Count']
    return counts

# 2. Map the 4 files to the 3 Biological Groups
files_and_groups = {
    'core_annotations.tabular': 'Core',
    'softcore_annotations.tabular': 'Accessory', 
    'shell_annotations.tabular': 'Accessory',    
    'cloud_annotations.tabular': 'Unique'        
}

df_list = []
for file_path, group_name in files_and_groups.items():
    if os.path.exists(file_path):
        df_part = parse_eggnog_annotations(file_path)
        df_part['Group'] = group_name
        df_list.append(df_part)
    else:
        print(f"Warning: {file_path} not found.")

final_df = pd.concat(df_list, ignore_index=True)

# 3. Aggregate counts by COG Letter and Group
grouped_df = final_df.groupby(['COG_Letter', 'Group'])['Count'].sum().reset_index()

# 4. Calculate True Pangenome Relative Abundance (%)
total_pangenome_cogs = grouped_df['Count'].sum()
grouped_df['Percentage'] = (grouped_df['Count'] / total_pangenome_cogs) * 100

grouped_df['COG_Letter'] = '[' + grouped_df['COG_Letter'] + ']'

# 5. Pivot data for stacking
stacked_data = grouped_df.pivot(index='COG_Letter', columns='Group', values='Percentage').fillna(0)

# Enforce correct drawing order from bottom to top
group_order = ['Core', 'Accessory', 'Unique']
existing_cols = [col for col in group_order if col in stacked_data.columns]
stacked_data = stacked_data[existing_cols]

# 6. Generate the Figure
fig, ax = plt.subplots(figsize=(10, 6))

# Custom Palette: Core (Orange), Accessory (Green), Unique (Yellow)
custom_palette = ['#ff7f0e', '#2ca02c', '#ffcf20'] 

stacked_data.plot(
    kind='bar', 
    stacked=True, 
    ax=ax, 
    color=custom_palette, 
    edgecolor='black', # Added thin borders back so the bright yellow doesn't bleed into the white background
    linewidth=0.5,
    width=0.65,
    zorder=3 
)

# 7. Customize text and layout
ax.set_title('COG Distribution', fontsize=16, pad=15)
ax.set_ylabel('% COGs', fontsize=14, labelpad=10)
ax.set_xlabel('') 

# Add horizontal gridlines behind the bars
ax.yaxis.grid(True, color='#e0e0e0', linestyle='-', linewidth=1.0, zorder=0)
ax.set_axisbelow(True)

plt.xticks(rotation=0, fontsize=11)
plt.yticks(fontsize=11)

# Clean legend positioned at the bottom center
plt.legend(
    loc='upper center',
    bbox_to_anchor=(0.5, -0.12),
    ncol=3, 
    frameon=False,
    fontsize=12,
    handletextpad=0.4
)

plt.tight_layout()

# 8. Export files
plt.savefig('COG_Final_Distribution.png', dpi=300, bbox_inches='tight')
plt.savefig('COG_Final_Distribution.pdf', dpi=300, bbox_inches='tight')
plt.savefig('COG_Final_Distribution.svg', format='svg', bbox_inches='tight')

print("Chart generated successfully: COG_Final_Distribution.png, .pdf, and .svg")