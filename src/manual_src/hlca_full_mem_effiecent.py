"""Memory efficient loading and filtering of HLCA dataset."""

import gc
import logging
import sys
from pathlib import Path

import scanpy as sc

# Set log files
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

if not adata_path.exists():
    raise FileNotFoundError(f"{adata_path} not found.")

logging.info(f"Reading AnnData file from {adata_path}")
adata = sc.read_h5ad(adata_path, backed="r")

# Filter before loading into memory
mask = adata.obs["disease"].isin(
    [
        "normal",
        "pulmonary fibrosis",
        "interstitial lung disease",
        "chronic obstructive pulmonary disease",
    ]
)

adata_subset = adata[mask, :].to_memory()

logging.info(f"Diseases: {adata_subset.obs['disease'].value_counts()}")
logging.info(
    f"Tissue types: {adata_subset.obs['tissue_coarse_unharmonized'].value_counts()}"
)

adata = None  # free memory
gc.collect()

logging.info("Saving filtered data...")
adata_subset.write_h5ad(filtered_path)

logging.info("HLCA plotting script completed successfully.")
