"""Script to concatenate metrics_summary.csv files from multiple output directories."""

from pathlib import Path

import pandas as pd

# Define variables


# Define the base directory
base_dir = Path(
    "/Volumes/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/data/xenium_raw"
)

output_dir = Path(
    "/Volumes/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/data/out/summary_metrics/csv"
)

output_dir.mkdir(parents=True, exist_ok=True)

# List to store dataframes
dfs = []

# Iterate through run directories, then look for "output-*" folders within each
for run_folder in base_dir.iterdir():
    if run_folder.is_dir() and "RUN_" in run_folder.name:
        print(f"Processing run folder: {run_folder.name}")

        # Look for "output-*" folders within this run folder
        for output_folder in run_folder.glob("output-*"):
            if output_folder.is_dir():
                csv_path = output_folder / "metrics_summary.csv"

                if csv_path.exists():
                    print(f"Reading: {csv_path}")
                    df = pd.read_csv(csv_path)

                    # Add columns to track which folders the data came from
                    df["run"] = run_folder.name.split("_")[-1]
                    df["run_folder"] = run_folder.name
                    df["source_folder"] = output_folder.name

                    dfs.append(df)
                else:
                    print(f"Warning: metrics_summary.csv not found in {output_folder}")
            else:
                print(f"Skipping non-directory: {output_folder}")
    else:
        print(f"Skipping folder (no RUN_ in name): {run_folder.name}")

# Concatenate all dataframes
if dfs:
    # Concatenate by shared columns (intersection of all columns)
    combined_df = pd.concat(dfs, axis=0, ignore_index=True, join="inner")

    # Extract Condition from region_name
    if "region_name" in combined_df.columns:
        combined_df["Condition"] = combined_df["region_name"].str.split("_").str[0]
        print(
            f"\n Condition column added: {combined_df['Condition'].unique().tolist()}"
        )
    else:
        print("\nWarning: 'region_name' column not found")

    print(f"\nSuccessfully concatenated {len(dfs)} files")
    print(f"Combined dataframe shape: {combined_df.shape}")
    print(f"\nColumns in combined dataframe:\n{combined_df.columns.tolist()}")

    # Save the combined dataframe
    output_path = output_dir / "combined_metrics_summary_all.csv"
    combined_df.to_csv(output_path, index=False)
    print(f"\nCombined data saved to: {output_path}")

    # Display first few rows
    print("\nFirst few rows:")
    print(combined_df.head())
else:
    print("No CSV files found to concatenate")
