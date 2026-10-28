import pandas as pd

# Load mapped outputs
ise_df = pd.read_csv("/mnt/d/Desktop/Project_Deinococcus/ISEScan/isescan_master_merged_mapped.tsv", sep="\t")
vir_df = pd.read_csv("/mnt/d/Desktop/Project_Deinococcus/Virsorter/virsorter_master_merged_mapped.tsv", sep="\t")
mef_df = pd.read_csv("/mnt/d/Desktop/Project_Deinococcus/MEF/mef_master_merged.csv")

# Count elements per strain
ise_counts = ise_df.groupby("Strain").size().rename("IS_Elements")
vir_counts = vir_df.groupby("Strain").size().rename("Prophages")
mef_counts = mef_df.groupby("Strain").size().rename("Mobile_Elements")

# Combine into master summary matrix
summary_df = pd.concat([ise_counts, vir_counts, mef_counts], axis=1).fillna(0).astype(int)
summary_df["Total_MGEs"] = summary_df.sum(axis=1)
summary_df = summary_df.reset_index()

# Export everything into a single multi-sheet Excel file
output_excel = "/mnt/d/Desktop/Project_Deinococcus/MGE_Summary_Master.xlsx"
with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="Summary_Matrix", index=False)
    ise_df.to_excel(writer, sheet_name="ISEScan_Data", index=False)
    vir_df.to_excel(writer, sheet_name="VirSorter_Data", index=False)
    mef_df.to_excel(writer, sheet_name="MEF_Data", index=False)

print(f"Master Excel generated successfully at: {output_excel}")