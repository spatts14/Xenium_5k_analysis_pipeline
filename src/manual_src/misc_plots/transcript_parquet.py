"""Combine all transcripts.parquet files into a single CSV file."""

import os
import re
from pathlib import Path

import pandas as pd

# HPC
BASE_DIR = Path(
    "/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1"
)
# LOCAL
# BASE_DIR = Path("/Volumes/phenotypingsputumasthmaticsaurorawellcomea1/")

XENIUM_DIR = BASE_DIR / "live/Sara_Patti/009_ST_Xenium/data/xenium_raw"
OUTPATH_DIR = BASE_DIR / "live/Sara_Patti/009_ST_Xenium/data/out/transcripts/"

GENE_LIST = ["16S", "MRC1", "KRT5"]

dfs = []

# Iterate through top-level directories in xenium/
for batch_dir in os.listdir(XENIUM_DIR):
    batch_path = os.path.join(XENIUM_DIR, batch_dir)
    if not os.path.isdir(batch_path):
        continue

    # Extract Batch: e.g. RUN_2 from 20251028__162706__SP25164_SARA_PATTI_RUN_2
    batch_match = re.search(r"(RUN_\d+)", batch_dir)
    batch = batch_match.group(1) if batch_match else batch_dir

    # Iterate through subdirectories
    for sub_dir in os.listdir(batch_path):
        sub_path = os.path.join(batch_path, sub_dir)
        if not os.path.isdir(sub_path):
            continue

        # Extract ROI: starts with COPD, IPF, PM08, or MICA
        roi_match = re.search(
            r"((?:COPD|IPF|PM08|MICA)[^_]*(?:_[^_]+)*?)(?=__\d{8})", sub_dir
        )
        roi = roi_match.group(1) if roi_match else sub_dir

        parquet_path = os.path.join(sub_path, "transcripts.parquet")
        if not os.path.exists(parquet_path):
            print(f"WARNING: No transcripts.parquet found in {sub_path}")
            continue

        print(f"Loading: {parquet_path}")
        print(f"  Batch: {batch} | ROI: {roi}")

        df = pd.read_parquet(parquet_path)
        df["Batch"] = batch
        df["ROI"] = roi

        # Filter to only keep transcripts in GENE_LIST
        n_before = len(df)
        df = df[df["feature_name"].isin(GENE_LIST)]
        n_after = len(df)
        print(
            f"  Transcripts: {n_before:,} → {n_after:,} (kept {n_after / n_before * 100:.1f}% matching gene list)"
        )

        if n_after == 0:
            print(f"  WARNING: No transcripts matched GENE_LIST in {sub_path}")
            continue

        dfs.append(df)

# Combine all dataframes
combined = pd.concat(dfs, ignore_index=True)

print(f"\nCombined shape: {combined.shape}")
print(f"Genes found: {combined['feature_name'].unique().tolist()}")
print(f"Batches found: {combined['Batch'].unique().tolist()}")
print(f"ROIs found: {combined['ROI'].unique().tolist()}")

# Save to CSV
output_path = OUTPATH_DIR / "xenium_transcripts_combined_gene_list.csv"
combined.to_csv(output_path, index=False)
print(f"\nSaved to: {output_path}")
