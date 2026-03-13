"""Subset adata based on metadata columns."""

# Subset adata based on metadata columns
import logging
import os
from pathlib import Path

import anndata as ad
import scanpy as sc


def subset_adata_by_meta(
    adata: ad.AnnData,
    meta_col: str,
    meta_values: list[str],
    logger=None,
) -> ad.AnnData:
    """Subset an AnnData object based on values in a metadata column.

    Args:
    adata : anndata.AnnData
        Input AnnData object.
    meta_col : str
        Column name in adata.obs to subset on.
    meta_values : list[str]
        Values in meta_col to retain.
    logger : logging.Logger, optional
        Logger for error reporting.

    Returns:
    anndata.AnnData
        Subsetted AnnData object containing only cells matching meta_values.
    """
    # Check column exists
    if meta_col not in adata.obs.columns:
        msg = f"Meta column '{meta_col}' not found in adata.obs.columns"
        if logger:
            logger.error(msg)
        raise ValueError(msg)

    # Check values exist
    available_values = set(adata.obs[meta_col].unique())
    missing_values = [val for val in meta_values if val not in available_values]

    if missing_values:
        msg = f"Meta values {missing_values} not found in meta column '{meta_col}'"
        if logger:
            logger.error(msg)
        raise ValueError(msg)

    # Subset
    mask = adata.obs[meta_col].isin(meta_values)
    adata_subset = adata[mask].copy()

    return adata_subset


# Set directories and file names
base_dir = Path(os.getenv("BASE_DIR"))
module_dir = base_dir / "annotate/level_1_annotation/"
h5ad_file = "adata_level_1_annotation1.h5ad"

output_dir = base_dir / "subset_adata"
os.makedirs(output_dir, exist_ok=True)

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

# Step 1: subset by condition
condition_col = "condition"
condition_values = ["IPF", "PM08", "COPD"]
logger.info(
    f"Subsetting data based on column '{condition_col}' and values {condition_values}"
)
adata_condition_subset = subset_adata_by_meta(
    adata,
    meta_col=condition_col,
    meta_values=condition_values,
    logger=logger,
)
logger.info(f"Condition subset complete. New shape: {adata_condition_subset.shape}")

# Step 2: subset the condition-filtered data by timepoint
timepoint_col = "timepoint"
timepoint_values = ["NA", "V1"]
logger.info(
    "Subsetting condition-filtered data based on column "
    f"'{timepoint_col}' and values {timepoint_values}"
)
adata_subset = subset_adata_by_meta(
    adata_condition_subset,
    meta_col=timepoint_col,
    meta_values=timepoint_values,
    logger=logger,
)
logger.info(f"Final subset complete. New shape: {adata_subset.shape}")


# Check output of subsetted adata
logger.info("Checking subsetted data...")
logger.info(f"Conditions: {adata_subset.obs['condition'].value_counts()}")
logger.info(f"Timepoints: {adata_subset.obs['timepoint'].value_counts()}")

# Save new adata
file_name = "adata_subset_COPD_IPF_PM08.h5ad"
logger.info(f"Saving subsetted data to {output_dir / file_name}")
adata_subset.write_h5ad(output_dir / file_name)
logger.info("Subsetted data saved successfully.")
