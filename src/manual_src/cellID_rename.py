"""Rename ROIs in adata.obs["ROI"] based on mapping from per-tissue CSV files."""

import glob
import os
import re
from pathlib import Path

import pandas as pd
import scanpy as sc

output = Path(
    "/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/"
)

# Folder containing the per-tissue cell stats CSVs
CSV_FOLDER = output / "docs/xenium_explorer_cell_IDs"
ADATA_PATH = output / "data/out/adata/all_samples.h5ad"
RENAME_ROIS = ["MICA_III_319_315_311", "MICA_III_325_337_379"]

## STEP 1:  Build mapping of cell_ID to ROI label from the per-tissue CSV files
# Set path for output

# Load csv files and build map of cell_ID to ROI label
cell_to_roi = {}
csv_files = sorted(glob.glob(os.path.join(CSV_FOLDER, "*.csv")))

if not csv_files:
    raise FileNotFoundError(
        f"No CSV files found in '{CSV_FOLDER}'. "
        "Check the folder path and make sure it contains *_cells_stats.csv files."
    )

print(f"Found {len(csv_files)} CSV file(s) in '{CSV_FOLDER}'. Processing...")

for csv_path in csv_files:
    # Read the Selection name from the comment header (line 1)
    with open(csv_path) as f:
        first_line = f.readline().strip()

    # Extract ROI name from "#Selection name: MICA_III_311"
    match = re.search(r"Selection name:\s*(\S+)", first_line)
    if match:
        # Normalise: MICAIII_319 -> MICA_III_319
        raw_name = match.group(1)
        roi_name = re.sub(r"MICAIII", "MICA_III", raw_name)
    else:
        # Fallback: derive from filename  e.g. MICA_III_311_cells_stats.csv
        basename = os.path.basename(csv_path)
        match2 = re.search(r"(MICA_III_\d+)", basename)
        if match2:
            roi_name = match2.group(1)
        else:
            roi_name = basename.replace(".csv", "").replace("_cells_stats", "")
            print(
                f"WARNING: Could not extract ROI name from '{basename}' using '{roi_name}'"
            )

    # Read the cell data (skip comment lines starting with #)
    df = pd.read_csv(csv_path, comment="#")
    df.columns = df.columns.str.strip()

    n_cells = len(df)
    print(f"{os.path.basename(csv_path)}")
    print(f"ROI label : {roi_name}")
    print(f"Cells     : {n_cells}\n")

    if "Cell ID" not in df.columns:
        print(
            f"WARNING: 'Cell ID' column not found in {csv_path}. Available columns: {list(df.columns)}"
        )
        continue

    for cell_id in df["Cell ID"]:
        cell_to_roi[cell_id.strip()] = roi_name

print(f"Total cell ROI mappings built: {len(cell_to_roi)}\n")


## STEP 2: Load adata and subset to only the RENAME_ROIS combined ROIs, then apply mapping to adata.obs["ROI"]

# Load adata
adata_path = ADATA_PATH
if not adata_path.exists():
    raise FileNotFoundError(
        f"adata file not found at '{adata_path}'. Check the path and filename."
    )
adata = sc.read_h5ad(adata_path)

before = adata.obs["ROI"].value_counts()
print("ROI value counts BEFORE remapping:")
print(before, "\n")

# Subset to only the RENAME_ROIS combined ROIs
rename_rois = adata[adata.obs["ROI"].isin(RENAME_ROIS)].copy()
print(f"Subset of ROIs to rename shape: {rename_rois.shape}")
print(f"Example obs_names: {rename_rois.obs_names[:5].tolist()}")

# STEP 3: Strip ROI name from cell ID
new_obs_names = list(adata.obs_names)  # start with full adata obs_names
new_roi_values = adata.obs["ROI"].tolist()  # start with full adata ROI values
unmatched = []

for i, obs_name in enumerate(adata.obs_names):
    current_roi = new_roi_values[i]

    # Only process MICA III cells
    if current_roi not in RENAME_ROIS:
        continue

    # Strip wrong prefix to recover bare cell ID
    if obs_name.startswith(current_roi + "_"):
        bare_cell_id = obs_name[len(current_roi) + 1 :]
    else:
        match = re.search(r"([a-z]+-\d+)$", obs_name)
        bare_cell_id = match.group(1) if match else None

    if bare_cell_id and bare_cell_id in cell_to_roi:
        correct_roi = cell_to_roi[bare_cell_id]
        new_obs_names[i] = f"{correct_roi}_{bare_cell_id}"  # update in place
        new_roi_values[i] = correct_roi  # update in place
    else:
        unmatched.append(obs_name)

# Apply corrections back onto original adata
adata.obs_names = new_obs_names
adata.obs["ROI"] = new_roi_values


# STEP 4: Validate and save corrected complete adata

if not adata.obs_names.is_unique:
    dupes = pd.Series(adata.obs_names)
    dupes = dupes[dupes.duplicated()]
    print(f"WARNING: {len(dupes)} duplicate obs_names:")
    print(dupes.values[:10])
else:
    print("obs_names are unique")

print(f"\nROI counts after fix:\n{adata.obs['ROI'].value_counts()}\n")

if unmatched:
    print(f"WARNING: {len(unmatched)} MICA III cells not matched in any CSV:")
    print(unmatched[:10])
else:
    print("All ROIs renamed cells successfully")

adata.write_h5ad(ADATA_PATH.with_name("all_samples_fixed.h5ad"))
print("Saved corrected complete adata.")
