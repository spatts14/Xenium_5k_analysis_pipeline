"""Filter and process the HLCA full reference dataset."""

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
        logging.StreamHandler(sys.stdout),  # prints to console (stdout)
        logging.FileHandler(log_file),  # also saves to file
    ],
)

logging.info("Starting HLCA plotting script...")

# Set directories and file paths
output_dir = Path("/rds/general/user/sep22/home/Projects/_Public_datasets/HLCA/data/")
adata_path = output_dir / "hlca_full_unprocessed.h5ad"
filtered_path = output_dir / "hlca_full_filtered.h5ad"
processed_path = output_dir / "hlca_full_processed.h5ad"

# Set figure directory for this module (overrides global setting)
sc.settings.figdir = output_dir / "figs" / "hlca_full_filt_process"
sc.settings.figdir.mkdir(parents=True, exist_ok=True)

# Check if adata file exists
if not adata_path.exists():
    raise FileNotFoundError(f"{adata_path} not found.")

logging.info(f"Loading adata file from {adata_path}...")
adata = sc.read_h5ad(adata_path)

logging.info("Filter for disease and tissue_level_2...")
for key in ["disease", "tissue_level_2", "cell_type"]:
    if key not in adata.obs:
        raise KeyError(f"{key} not found in adata.obs columns.")

# Filter for specific diseases
keep_diseases = [
    "normal",
    "pulmonary fibrosis",
    "interstitial lung disease",
    "chronic obstructive pulmonary disease",
]
adata = adata[adata.obs["disease"].isin(keep_diseases), :].copy()
logging.info(f"Remaining disease: {adata.obs['disease'].unique()}")

# Filter for specific tissues
remove = ["nose", "inferior turbinate", "trachea"]
adata = adata[~adata.obs["tissue_level_2"].isin(remove), :].copy()
logging.info(
    f"Removed {remove}. Remaining tissue types: {adata.obs['tissue_level_2'].unique()}"
)

# Confirm filtering
logging.info(f"Disease value counts:\n{adata.obs['disease'].value_counts()}")
logging.info(f"Tissue value counts:\n{adata.obs['tissue_level_2'].value_counts()}")

adata.write_h5ad(filtered_path)
logging.info(f"Saved filtered adata to {filtered_path}...")

logging.info("Process filtered dataset...")

# Basic filtering and normalization
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)

# Store raw counts
if "counts" not in adata.layers:
    if adata.raw is not None:
        # Subset raw to match current adata var_names
        adata.layers["counts"] = adata.raw[:, adata.var_names].X.copy()
        logging.info("Stored raw counts in adata.layers['counts'].")
    else:
        logging.info("adata.raw is missing. Cannot store raw counts.")
else:
    logging.info("Raw counts layer already exists. Skipping.")

# Check if data is already normalized
logging.info("Checking normalization...")
totals = adata.X.sum(axis=1)
totals = totals.A1 if hasattr(totals, "A1") else np.array(totals).ravel()

if np.median(totals) > 2e4:  # crude but works for most scRNA-seq
    logging.info("Data does not appear to be normalized! Normalizing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    logging.info("Normalized total counts.")
else:
    logging.info("Data already appears normalized. Skipping.")

# Check if data is log-transformed
logging.info("Checking log transformation...")
X_sample = adata.X[:100, :100]
arr = X_sample.A if hasattr(X_sample, "A") else X_sample

if np.allclose(arr, np.round(arr)):  # looks like raw counts
    logging.info("Data does not appear to be log transformed! Transforming data...")
    sc.pp.log1p(adata)
    logging.info("Applied log1p transformation.")
else:
    logging.info("Data already appears log-transformed. Skipping.")

# Feature selection and scaling
sc.pp.highly_variable_genes(adata, n_top_genes=5000, flavor="seurat_v3")
sc.pp.scale(adata, max_value=10)

# Dimensionality reduction
sc.tl.pca(adata, n_comps=75, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=75)
sc.tl.umap(adata)

# Save processed data
adata.write_h5ad(processed_path)
logging.info(f"Saved processed adata to {processed_path}...")

# Plot UMAP

# Plot UMAP
obs_lists = ["disease", "tissue_coarse_unharmonized", "cell_type"]
for obs in obs_lists:
    sc.pl.umap(
        adata,
        color=obs,
        title=f"HLCA full reference data UMAP colored by {obs}",
        save=f"_hlca_full_ref_{obs}.png",
    )

logging.info("Filtering and processing script completed.")
