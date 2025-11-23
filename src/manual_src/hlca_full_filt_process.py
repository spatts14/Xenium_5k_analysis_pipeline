"""Filter and process the HLCA full reference dataset - Memory Optimized."""

import logging
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

log_file = "hlca_analysis.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_file),
    ],
)

logging.info("Starting HLCA plotting script...")

# Set directories and file paths
output_dir = Path("/rds/general/user/sep22/home/Projects/_Public_datasets/HLCA/data/")
adata_path = output_dir / "hlca_full_unprocessed.h5ad"
filtered_path = output_dir / "hlca_full_filtered.h5ad"
processed_path = output_dir / "hlca_full_processed.h5ad"

# Set figure directory for this module
sc.settings.figdir = output_dir / "figs" / "hlca_full_filt_process"
sc.settings.figdir.mkdir(parents=True, exist_ok=True)

# Check if adata file exists
if not adata_path.exists():
    raise FileNotFoundError(f"{adata_path} not found.")

logging.info(f"Loading adata file from {adata_path}...")
# Use backed mode to avoid loading entire dataset into memory
adata = sc.read_h5ad(adata_path, backed="r")
logging.info(f"Original shape: {adata.shape}")

logging.info("Filter for disease and tissue_level_2...")
for key in ["disease", "tissue_level_2", "cell_type"]:
    if key not in adata.obs:
        raise KeyError(f"{key} not found in adata.obs columns.")

# Create filter mask
keep_diseases = [
    "normal",
    "pulmonary fibrosis",
    "interstitial lung disease",
    "chronic obstructive pulmonary disease",
]
remove_tissues = ["nose", "inferior turbinate", "trachea"]

disease_mask = adata.obs["disease"].isin(keep_diseases)
tissue_mask = ~adata.obs["tissue_level_2"].isin(remove_tissues)
combined_mask = disease_mask & tissue_mask

logging.info(f"Cells passing filters: {combined_mask.sum()} / {len(combined_mask)}")

# Load only filtered subset into memory
adata = adata[combined_mask, :].to_memory()
logging.info(f"Filtered shape: {adata.shape}")

logging.info(f"Remaining diseases: {adata.obs['disease'].unique()}")
logging.info(f"Remaining tissue types: {adata.obs['tissue_level_2'].unique()}")

# Confirm filtering
logging.info(f"Disease value counts:\n{adata.obs['disease'].value_counts()}")
logging.info(f"Tissue value counts:\n{adata.obs['tissue_level_2'].value_counts()}")

# Save filtered dataset
logging.info(f"Saving filtered adata to {filtered_path}...")
adata.write_h5ad(filtered_path)

logging.info("Process filtered dataset...")

# Basic filtering
logging.info("Filtering cells and genes...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
logging.info(f"Shape after filtering: {adata.shape}")

# Store raw counts BEFORE normalization
if "counts" not in adata.layers:
    if adata.raw is not None:
        adata.layers["counts"] = adata.raw[:, adata.var_names].X.copy()
        logging.info("Stored raw counts in adata.layers['counts'].")
    else:
        # If no raw, store current X if it's counts
        totals = adata.X.sum(axis=1)
        totals = totals.A1 if hasattr(totals, "A1") else np.array(totals).ravel()
        if np.median(totals) > 1e3:  # Looks like raw counts
            adata.layers["counts"] = adata.X.copy()
            logging.info("Stored current X as raw counts.")
        else:
            logging.warning("No raw counts available!")
else:
    logging.info("Raw counts layer already exists.")

# Check if data is already normalized
logging.info("Checking normalization...")
totals = adata.X.sum(axis=1)
totals = totals.A1 if hasattr(totals, "A1") else np.array(totals).ravel()

if np.median(totals) > 2e4:
    logging.info("Normalizing data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
else:
    logging.info("Data already normalized.")

# Check if data is log-transformed
logging.info("Checking log transformation...")
X_sample = adata.X[:100, :100]
arr = X_sample.A if hasattr(X_sample, "A") else X_sample

if np.allclose(arr, np.round(arr)):
    logging.info("Applying log transformation...")
    sc.pp.log1p(adata)
else:
    logging.info("Data already log-transformed.")

# Dimensionality reduction
logging.info("Calculating PCA and UMAP...")
sc.tl.pca(adata, n_comps=75, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=75)
sc.tl.umap(adata)

# Save processed data
adata.write_h5ad(processed_path)
logging.info(f"Saved processed adata to {processed_path}...")

# Plot UMAP
obs_lists = ["disease", "tissue_level_2", "cell_type"]
for obs in obs_lists:
    sc.pl.umap(
        adata,
        color=obs,
        title=f"HLCA full reference data UMAP colored by {obs}",
        save=f"_hlca_full_ref_{obs}.png",
    )

logging.info("Filtering and processing script completed.")
