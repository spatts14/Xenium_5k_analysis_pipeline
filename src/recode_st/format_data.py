"""Import and format Xenium data into Zarr and AnnData formats."""

import re
import warnings
from logging import getLogger
from pathlib import Path

import anndata as ad
from spatialdata_io import xenium

from recode_st.config import IOConfig

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def extract_roi_name(folder_name: str) -> str:
    """Extract ROI name from folder name."""
    parts = folder_name.split("__")
    for part in parts:
        if part.startswith(("IPF", "COPD", "PM08")):
            return part
    raise ValueError(f"Could not parse ROI name from: {folder_name}")


def extract_run_number(folder_name: str) -> int:
    """Extract run number from folder name, e.g. RUN_1 → 1."""
    match = re.search(r"RUN_(\d+)", folder_name, re.IGNORECASE)
    if match:
        return int(match.group(1))
    raise ValueError(f"Could not extract run number from folder: {folder_name}")


def convert_all_xenium(io_config: IOConfig):
    """Convert all Xenium ROIs into Zarrs and one combined AnnData."""
    xenium_root = Path(io_config.xenium_dir)
    if not xenium_root.exists():
        raise FileNotFoundError(f"Xenium root directory not found: {xenium_root}")

    # Prepare output directories
    zarr_root = Path(io_config.zarr_dir)
    zarr_root.mkdir(parents=True, exist_ok=True)

    output_data_dir = Path(io_config.output_data_dir)
    output_data_dir.mkdir(parents=True, exist_ok=True)

    output_adata_dir = output_data_dir / "adata"
    output_adata_dir.mkdir(parents=True, exist_ok=True)

    all_adatas = []

    # Loop over each run folder directly inside xenium_raw/
    for run_dir in sorted(xenium_root.iterdir()):
        if not run_dir.is_dir():
            continue

        try:
            run_number = extract_run_number(run_dir.name)
        except ValueError:
            logger.warning(f"Skipping folder (no run number): {run_dir.name}")
            continue

        logger.info(f"Processing run {run_number}: {run_dir.name}")

        # Look for ROI folders inside this run directory
        for roi_folder in run_dir.iterdir():
            if not roi_folder.is_dir() or not roi_folder.name.startswith("output-"):
                continue

            try:
                roi_name = extract_roi_name(roi_folder.name)
            except ValueError as err:
                logger.warning(f"{err}")
                continue

            logger.info(f"Processing ROI: {roi_name}")

            try:
                sdata = xenium(roi_folder)
            except Exception as err:
                logger.error(f"Failed loading Xenium data for {roi_name}: {err}")
                continue

            # Save per-ROI Zarr
            roi_zarr_path = zarr_root / f"{roi_name}.zarr"
            try:
                sdata.write(roi_zarr_path, overwrite=True)
            except Exception as err:
                logger.error(f"Failed writing Zarr for {roi_name}: {err}")
                continue

            # Extract AnnData and annotate
            try:
                adata = sdata.tables["table"]
                adata.obs["ROI"] = roi_name
                adata.obs["batch"] = run_dir.name
                adata.obs["run"] = run_number
                adata.obs["condition"] = (
                    roi_name.split("_")[0]
                    if roi_name.split("_")[0] in ["IPF", "COPD", "PM08"]
                    else "Unknown"
                )
                all_adatas.append(adata)
            except KeyError:
                logger.warning(f"No table found in {roi_name}, skipping AnnData.")

    # Combine and save all AnnData objects
    if all_adatas:
        combined = ad.concat(all_adatas, join="outer", label="ROI", fill_value=0)
        combined_path = output_adata_dir / "all_samples.h5ad"
        combined.write(combined_path)
        logger.info(f"Combined AnnData written to {combined_path}")
    else:
        logger.warning("No AnnData objects were created. Combined AnnData not written.")


def run_format(io_config: IOConfig):
    """Run the formatting module across multiple ROIs."""
    convert_all_xenium(io_config)
    logger.info("Finished formatting all ROIs.")
