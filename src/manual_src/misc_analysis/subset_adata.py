"""Subset adata based on metadata columns."""

# Subset adata based on metadata columns
import logging
import os
from pathlib import Path

import scanpy as sc

# Set directories and file names
base_dir = Path(
    "/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/output/2026-02-22_analysis_run_HVG2000/"
)
module_dir = base_dir / "annotate"
h5ad_file = "adata.h5ad"

output_dir = module_dir / "subset_adata"
os.mkdir(output_dir, exist_ok=True)

# Set up logging
log_file = module_dir / "subsetting_by_meta.log"
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s - %(levelname)s] %(message)s",
    handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


# Load annotated adata
logger.info(f"Loading data from {module_dir / f'{h5ad_file}'}")
adata = sc.read_h5ad(module_dir / f"{h5ad_file}")
logger.info(f"Data loaded successfully. Shape: {adata.shape}")

# Meta data to subset by
# Column name
meta_col = "condition"
# Values to subset by
meta_values = ["IPF", "PM08", "COPD_V1"]

# Confirm that meta column exists
if meta_col not in adata.obs.columns:
    logger.error(f"Meta column '{meta_col}' not found in adata.obs.columns")
    raise ValueError(f"Meta column '{meta_col}' not found in adata.obs.columns")

# Confirm that meta values exist in the meta column
missing_values = [val for val in meta_values if val not in adata.obs[meta_col].unique()]
if missing_values:
    logger.error(f"Meta values {missing_values} not found in meta column '{meta_col}'")
    raise ValueError(
        f"Meta values {missing_values} not found in meta column '{meta_col}'"
    )

# Subset data based on meta column and values
logger.info(f"Subsetting data based on column '{meta_col}' and values {meta_values}")
adata_subset = adata[adata.obs[meta_col].isin(meta_values)].copy()
logger.info(f"Data subsetted successfully. New shape: {adata_subset.shape}")

# Name of new adata file
new_adata_file = f"adata_subset_by_{meta_col}.h5ad"

# Save new adata
logger.info(f"Saving subsetted data to {output_dir / new_adata_file}")
adata_subset.write_h5ad(output_dir / new_adata_file)
logger.info("Subsetted data saved successfully.")
