"""Import and format Xenium data into Zarr and AnnData formats."""

import warnings
from logging import getLogger
from pathlib import Path

import anndata as ad
from spatialdata_io import xenium

from recode_st.config import IOConfig
from recode_st.logging_config import configure_logging

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def extract_roi_name(folder_name: str) -> str:
    """Extract ROI name from folder name."""
    # Example: "output-...__IPF_RBH_19__20251001__141533" → "IPF_RBH_19"
    parts = folder_name.split("__")
    for part in parts:
        if part.startswith(("IPF", "COPD", "PM08")):
            return part
    raise ValueError(f"Could not parse ROI name from: {folder_name}")


def convert_all_xenium(io_config: IOConfig):
    """Convert all Xenium ROIs into Zarrs and one combined AnnData."""
    xenium_path = Path(io_config.xenium_dir)

    zarr_root = Path(io_config.zarr_dir)
    zarr_root.mkdir(parents=True, exist_ok=True)

    output_data_dir = Path(io_config.output_data_dir)
    output_data_dir.mkdir(parents=True, exist_ok=True)

    output_adata = output_data_dir / "adata"
    output_adata.mkdir(parents=True, exist_ok=True)

    all_adatas = []

    if not xenium_path.exists():
        raise FileNotFoundError(f"Xenium root directory not found: {xenium_path}")

    # Loop over run directories inside Xenium root
    for run_dir in xenium_path.iterdir():
        if not run_dir.is_dir() or not run_dir.name.lower().startswith("run_"):
            continue

        # Loop over date/run folders inside each run directory
        for date_dir in run_dir.iterdir():
            if not date_dir.is_dir():
                continue

            run_name = (
                date_dir.name
            )  # e.g. "20251001__141239__SP25164_SARA_PATTI_RUN_1"
            logger.info(f"Processing run: {run_name}")

            # Loop over ROI folders
            for roi_folder in date_dir.iterdir():
                if roi_folder.is_dir() and roi_folder.name.startswith("output-"):
                    roi_name = extract_roi_name(roi_folder.name)
                    logger.info(f"Processing ROI: {roi_name}")

                    try:
                        sdata = xenium(roi_folder)
                    except Exception as err:
                        logger.error(
                            f"Failed loading Xenium data for {roi_name}: {err}"
                        )
                        continue

                    # Save Zarr per ROI
                    roi_zarr_path = zarr_root / f"{roi_name}.zarr"
                    try:
                        sdata.write(roi_zarr_path, overwrite=True)
                    except Exception as err:
                        logger.error(f"Failed writing Zarr for {roi_name}: {err}")
                        continue

                    # Extract AnnData and add labels
                    try:
                        adata = sdata.tables["table"]
                        adata.obs["ROI"] = roi_name
                        adata.obs["batch"] = run_name
                        all_adatas.append(adata)
                    except KeyError:
                        logger.warning(
                            f"No table found in {roi_name}, skipping AnnData."
                        )

    # Concatenate all adatas
    if all_adatas:
        combined = ad.concat(all_adatas, join="outer", label="ROI", fill_value=0)
        combined_path = output_adata / "all_samples.h5ad"
        combined.write(combined_path)
        logger.info(f"Combined AnnData written to {combined_path}")
    else:
        logger.warning("No AnnData objects were created. Combined AnnData not written.")


def run_format(io_config: IOConfig):
    """Run the formatting module across multiple ROIs."""
    convert_all_xenium(io_config)
    logger.info("Finished formatting all ROIs.")


if __name__ == "__main__":
    configure_logging()
    logger = getLogger("recode_st.0_format")
    run_format(IOConfig())
