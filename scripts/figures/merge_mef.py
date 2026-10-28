import os
import glob
import pandas as pd

# Automatically navigate to your MEF directory in WSL
target_dir = "/mnt/d/Desktop/Project_Deinococcus/MEF/MGE_Outputs"

if os.path.exists(target_dir):
    os.chdir(target_dir)
print(f"Working Directory: {os.getcwd()}")

# Search for all CSV files recursively
search_pattern = "**/*.csv"
files = glob.glob(search_pattern, recursive=True)

# Exclude master output file if re-running
files = [f for f in files if "mef_master_merged" not in f]
print(f"Found {len(files)} matching MEF CSV files.\n")

dfs = []
for file in files:
    print(f"Processing: {file}")
    strain_name = os.path.basename(file).replace(".csv", "")
    
    try:
        # Open file to find where the actual table starts (skipping text headers)
        start_row = 0
        with open(file, 'r') as f:
            for i, line in enumerate(f):
                # MEF tables have multiple columns. We look for a line with commas.
                if line.count(',') > 5: 
                    start_row = i
                    break
                    
        # Load the CSV starting exactly at the correct table header
        df = pd.read_csv(file, skiprows=start_row)
        df.insert(0, "Strain", strain_name)
        dfs.append(df)
        
    except Exception as e:
        print(f"  --> Skipping due to error: {e}")

if dfs:
    master_df = pd.concat(dfs, ignore_index=True)
    output_file = "mef_master_merged.csv"
    master_df.to_csv(output_file, index=False)
    print(f"\nSUCCESS: Merged {len(dfs)} files into '{output_file}'")
else:
    print("\nERROR: No CSV files successfully processed.")