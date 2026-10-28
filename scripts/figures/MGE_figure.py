import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load Strain Mapping
mapping_df = pd.read_csv("dataset_mapping.tsv", sep="\t")
mapping_dict = {
    f"dataset_{re.search(r'\d+', str(k)).group()}": str(v).strip()
    for k, v in zip(mapping_df["Dataset_ID"], mapping_df["Strain"])
    if re.search(r"\d+", str(k))
}

def normalize_strain(val):
    s_str = str(val).strip()
    match = re.search(r"dataset_(\d+)", s_str, re.IGNORECASE)
    if match:
        key = f"dataset_{match.group(1)}"
        if key in mapping_dict:
            return mapping_dict[key]
    clean_val = re.sub(r"[^a-zA-Z0-9]", "", s_str).lower()
    for std_strain in mapping_dict.values():
        if clean_val == re.sub(r"[^a-zA-Z0-9]", "", std_strain).lower():
            return std_strain
    return s_str

# 2. Load Raw Datasets
ise_df = pd.read_csv("/mnt/d/Desktop/Project_Deinococcus/ISEScan/isescan_master_merged_mapped.tsv", sep="\t")
vir_df = pd.read_csv("/mnt/d/Desktop/Project_Deinococcus/Virsorter/virsorter_master_merged_mapped.tsv", sep="\t")
mef_df = pd.read_csv("/mnt/d/Desktop/Project_Deinococcus/MEF/mef_master_merged.csv")
plasmid_df = pd.read_csv("/mnt/d/Desktop/Project_Deinococcus/PlasmidFinder/plasmid_summary_with_coordinates.tsv", sep="\t")

# 3. Normalize Strains and enforce 'ZPb4' cleanly in-memory (leaving raw files untouched)
for df in [ise_df, vir_df, mef_df, plasmid_df]:
    target_col = "Strain" if "Strain" in df.columns else df.columns[0]
    df["Strain"] = df[target_col].apply(normalize_strain)
    # Safely replace ZPb-04 with ZPb4 for outputs/figures only
    df["Strain"] = df["Strain"].replace({r"ZPb-04": "ZPb4", r"zpb-04": "ZPb4"}, regex=True)

# Filter simple IS elements out of MEF using the 'type' column
if "type" in mef_df.columns:
    mef_filtered = mef_df[~mef_df["type"].astype(str).str.contains(r'\bIS\b|insertion', case=False, na=False, regex=True)]
else:
    print("Warning: 'type' column not found in MEF data. Skipping IS filtering.")
    mef_filtered = mef_df

# 4. Count, Merge, and Clean Column Names
counts = [
    ise_df.groupby("Strain").size().rename("IS_Elements"),
    vir_df.groupby("Strain").size().rename("Prophages"),
    mef_filtered.groupby("Strain").size().rename("Transposons_ICEs"),
    plasmid_df.groupby("Strain").size().rename("Plasmids")
]

summary_df = pd.concat(counts, axis=1).fillna(0).astype(int)
summary_df["Total_MGEs"] = summary_df.sum(axis=1)

# Rename columns for clean publication labels
summary_df = summary_df.rename(columns={
    "IS_Elements": "IS elements",
    "Transposons_ICEs": "Transposons/ICEs",
    "Total_MGEs": "Total MGEs"
})

summary_df = summary_df.reset_index().sort_values("Total MGEs", ascending=True)

# 5. Export Updated Excel Matrix
summary_df.to_excel("/mnt/d/Desktop/Project_Deinococcus/MGE_Summary_Master_Updated.xlsx", index=False)

# 6. Generate Stacked Bar Chart
plt.figure(figsize=(12, 12))
p1 = plt.barh(summary_df["Strain"], summary_df["IS elements"], label="IS elements", color="#2b5c8f")
p2 = plt.barh(summary_df["Strain"], summary_df["Prophages"], left=summary_df["IS elements"], label="Prophages", color="#d95f02")
p3 = plt.barh(summary_df["Strain"], summary_df["Transposons/ICEs"], left=summary_df["IS elements"] + summary_df["Prophages"], label="Transposons/ICEs", color="#7570b3")
p4 = plt.barh(summary_df["Strain"], summary_df["Plasmids"], left=summary_df["IS elements"] + summary_df["Prophages"] + summary_df["Transposons/ICEs"], label="Plasmids", color="#e7298a")

plt.xlabel("Total Count", fontsize=12, fontweight="bold")
plt.ylabel("Strains", fontsize=12, fontweight="bold")
plt.title("Macro-Landscape of Mobile Genetic Elements", fontsize=14, fontweight="bold", pad=15)
plt.legend(title="MGE Type", frameon=True, loc="lower right")
plt.tight_layout()
plt.savefig("Figure_1_MGE_StackedBar.png", dpi=300)
plt.close()

# 7. Generate Heatmap
heatmap_data = summary_df.set_index("Strain")[["IS elements", "Prophages", "Transposons/ICEs", "Plasmids"]]
plt.figure(figsize=(9, 12))
sns.heatmap(heatmap_data, annot=True, fmt="d", cmap="YlGnBu", cbar_kws={'label': 'Count'}, linewidths=0.5)
plt.title("MGE Category Abundance", fontsize=14, fontweight="bold", pad=15)
plt.ylabel("Strains", fontsize=12, fontweight="bold")
plt.tight_layout()
plt.savefig("Figure_1_MGE_Heatmap.png", dpi=300)
plt.close()

print("Success: Excel matrix updated and Figure 1 generated with ZPb4!")