"""Script to read and visualize the combined metrics summary CSV file."""

from pathlib import Path

import pandas as pd

# Define the base directory
base_dir = Path("/Volumes/sep22/home/wet_lab/_Experiments/009_ST_Xenium/out/data")

# Save the combined dataframe
df = pd.read_csv(base_dir / "combined_metrics_summary.csv")
print(df.head())
