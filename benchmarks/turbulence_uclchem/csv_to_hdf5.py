import pandas as pd
import glob
import os
import numpy as np

folder_path = "/scratch/09338/carlos9/turbulence_chemistry/chemistry_tracers_csv"

csv_files = glob.glob(os.path.join(folder_path, "*.csv"))

dfs = []

for file in csv_files:
    print(f"Loading {file}...")

    # Skip if file is empty
    if os.path.getsize(file) == 0:
        print(f"Skipping {file} (empty file)")
        continue

    try:
        df = pd.read_csv(file)
    except pd.errors.EmptyDataError:
        print(f"Skipping {file} (no valid CSV data)")
        continue
    
    unique_vals = df['tracer'].dropna().unique()
    if len(unique_vals) == 1:
        fill_value = unique_vals[0]
        df['tracer'] = df['tracer'].fillna(fill_value)
    else:
        print(f"Warning: more than one non-NaN tracer value in {file} → {unique_vals}")

    print(df['tracer'].unique())
    dfs.append(df)


combined_df = pd.concat(dfs, ignore_index=True)

output_file = "/scratch/09338/carlos9/turbulence_chemistry/turbulence_benchmark.h5"
combined_df.to_hdf(output_file, key="data", mode="w")

print(f"Saved combined DataFrame to {output_file}")
