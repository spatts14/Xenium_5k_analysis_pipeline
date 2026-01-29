"""Import and format Xenium data into Zarr and AnnData formats."""

import warnings
from logging import getLogger
from pathlib import Path

import anndata as ad
from spatialdata_io import xenium

from recode_st.config import IOConfig

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)


def convert_all_xenium(io_config: IOConfig):
    """Convert all Xenium ROIs into Zarrs and one combined AnnData."""
    xenium_path = Path(io_config.xenium_dir)

    zarr_root = Path(io_config.zarr_dir)
    zarr_root.mkdir(parents=True, exist_ok=True)

    output_data_dir = Path(io_config.output_data_dir)
    output_data_dir.mkdir(parents=True, exist_ok=True)

    adata_dir = Path(io_config.adata_dir)
    adata_dir.mkdir(parents=True, exist_ok=True)

    all_adatas = []

    logger.info(f"Starting conversion of Xenium data in {xenium_path}")

    if not xenium_path.exists():
        raise FileNotFoundError(f"Xenium root directory not found: {xenium_path}")

    # Get all run directories and log them
    all_run_dirs = [d for d in xenium_path.iterdir() if d.is_dir()]
    logger.info(f"Found {len(all_run_dirs)} run directories:")
    for run_dir in sorted(all_run_dirs):
        logger.info(f"  - {run_dir.name}")

    # Loop over run/date folders (e.g., "20260122__165134__SP25164_SARA_PATTI_RUN_4")
    for run_dir in xenium_path.iterdir():
        if not run_dir.is_dir():
            continue

        run_name = run_dir.name  # e.g., "20260122__165134__SP25164_SARA_PATTI_RUN_4"
        logger.info(f"Processing run directory: {run_name}")

        # Count and log output-* folders within this run
        output_folders = [
            f for f in run_dir.iterdir() if f.is_dir() and f.name.startswith("output-")
        ]
        logger.info(f"Found {len(output_folders)} output-* folders in {run_name}")

        # Loop over each output-* folder
        for output_folder in output_folders:
            folder_name = output_folder.name
            # e.g., "output-XETG00431__0051416__COPD_46005_V1__20260122__165617"
            #   or  "output-XETG00431__0051416__IPF_RBH_06__20260122__165618"

            logger.info(f"Processing folder: {folder_name}")

            # Extract ROI information from folder name
            try:
                # Split by "__" and get the third part
                parts = folder_name.split("__")
                if len(parts) < 3:
                    logger.warning(f"Unexpected folder format: {folder_name}, skipping")
                    continue

                roi_info = parts[2]  # e.g., "COPD_46005_V1" or "IPF_RBH_06"
                roi_name = roi_info  # Use the full roi_info as roi_name

                roi_parts = roi_info.split("_")

                # Check if last part is a valid timepoint
                valid_timepoints = ["V1", "V2", "V3"]

                if roi_parts[-1] in valid_timepoints:
                    # Has timepoint: COPD_46005_V1
                    timepoint = roi_parts[-1]  # "V1"
                    sample_parts = roi_parts[:-1]  # ["COPD", "46005"]
                    condition = sample_parts[0]  # "COPD"
                    sample_ID = "_".join(sample_parts)  # "COPD_46005"
                else:
                    # No timepoint: IPF_RBH_06
                    timepoint = "NA"
                    condition = roi_parts[0]  # "IPF"
                    sample_ID = roi_info  # "IPF_RBH_06"

                logger.info(
                    f"Parsed - ROI: {roi_name}, Condition: {condition}, "
                    f"Sample ID: {sample_ID}, Timepoint: {timepoint}"
                )

            except Exception as err:
                logger.error(f"Failed to parse folder name {folder_name}: {err}")
                continue

            # Load Xenium data
            try:
                sdata = xenium(output_folder)
                logger.info(f"Loaded Xenium data for {roi_name}")
            except Exception as err:
                logger.error(f"Failed loading Xenium data for {roi_name}: {err}")
                continue

            # Save Zarr per ROI
            roi_zarr_path = zarr_root / f"{roi_name}.zarr"
            try:
                sdata.write(roi_zarr_path, overwrite=True)
                logger.info(f"Zarr for {roi_name} saved to {roi_zarr_path}")
            except Exception as err:
                logger.error(f"Failed writing Zarr for {roi_name}: {err}")
                continue

            # Extract AnnData and add labels
            try:
                # Extract AnnData from SpatialData
                logger.info(f"Creating AnnData for {roi_name}")
                adata = sdata.tables["table"]

                # Validate AnnData structure
                if adata.X is None:
                    logger.error(f"AnnData for {roi_name} has no .X matrix, skipping")
                    continue

                logger.info(f"AnnData for {roi_name} has shape {adata.shape}")

                # Add metadata to AnnData
                adata.obs["ROI"] = roi_name  # "COPD_46005_V1" or "IPF_RBH_06"
                adata.obs["batch"] = run_name.split("_")[
                    -1
                ]  # "4" (extracts from "RUN_4")
                adata.obs["condition"] = condition  # "COPD" or "IPF"
                adata.obs["sample_ID"] = sample_ID  # "COPD_46005" or "IPF_RBH_06"
                adata.obs["timepoint"] = timepoint  # "V1" or "NA"

                # Add timepoint label
                timepoint_mapping = {
                    "V1": "baseline",
                    "V2": "6_weeks",
                    "V3": "6_months",
                }
                adata.obs["timepoint_label"] = timepoint_mapping.get(timepoint, "NA")

                logger.info(
                    f"Added metadata to AnnData - ROI: {roi_name}, "
                    f"Batch: {adata.obs['batch'].iloc[0]}, "
                    f"Condition: {condition}, Sample ID: {sample_ID}, "
                    f"Timepoint: {timepoint} ({adata.obs['timepoint_label'].iloc[0]})"
                )

                # Rename obs_names to include ROI so they are unique across ROIs
                adata.obs_names = roi_name + "_" + adata.obs_names.astype(str)

                if not adata.obs_names.is_unique:
                    raise ValueError(f"Non-unique obs_names in AnnData for {roi_name}")

                logger.info(f"obs_names are unique in AnnData for {roi_name}")
                logger.info(f"Example obs_names: {adata.obs_names[:3].tolist()}")
                logger.info(
                    f"AnnData for {roi_name} has shape {adata.shape} "
                    f"with {adata.n_vars} genes"
                )

                # Final validation before saving
                if adata.X is None:
                    logger.error(f"AnnData for {roi_name} lost .X matrix before saving")
                    continue

                # Save individual AnnData
                adata_path = adata_dir / f"{roi_name}.h5ad"
                logger.info(
                    f"Saving AnnData for {roi_name} with .X shape: {adata.X.shape}"
                )
                adata.write(adata_path)
                logger.info(f"AnnData for {roi_name} saved to {adata_path}")

                all_adatas.append(adata)

            except KeyError:
                logger.warning(f"No table found in {roi_name}, skipping AnnData")
            except Exception as err:
                logger.error(f"Failed processing AnnData for {roi_name}: {err}")
                continue

    logger.info(f"Conversion complete. Processed {len(all_adatas)} ROIs successfully")

    # Concatenate all adatas (MOVED BEFORE return)
    if all_adatas:
        combined = ad.concat(all_adatas, join="outer", fill_value=0)
        combined_path = adata_dir / "all_samples.h5ad"
        combined.write(combined_path)
        logger.info(f"Combined AnnData written to {combined_path}")
    else:
        logger.warning("No AnnData objects were created. Combined AnnData not written.")

    return all_adatas


def concatenate_h5ad_files(adata_dir: SyntaxWarning):
    """Load all .h5ad files from a directory and concatenate them.

    Args:
        adata_dir : str
            Path to directory containing individual .h5ad files

    Returns:
        ad.AnnData
            Concatenated AnnData object
    """
    adata_dir = Path(adata_dir)

    if not adata_dir.exists():
        raise FileNotFoundError(f"Directory not found: {adata_dir}")

    # Find all .h5ad files (excluding any existing all_samples.h5ad)
    h5ad_files = [f for f in adata_dir.glob("*.h5ad") if f.name != "all_samples.h5ad"]

    if not h5ad_files:
        raise ValueError(f"No .h5ad files found in {adata_dir}")

    logger.info(f"Found {len(h5ad_files)} .h5ad files to concatenate")

    # Load all AnnData objects
    all_adatas = []
    failed_files = []

    for h5ad_file in sorted(h5ad_files):
        try:
            logger.info(f"Loading {h5ad_file.name}...")
            adata = ad.read_h5ad(h5ad_file)

            # Verify .X is not None
            if adata.X is None:
                logger.error(f"{h5ad_file.name} has None for .X matrix, skipping")
                failed_files.append(h5ad_file.name)
                continue

            logger.info(f"  Shape: {adata.shape}, .X shape: {adata.X.shape}")
            all_adatas.append(adata)

        except Exception as e:
            logger.error(f"Failed to load {h5ad_file.name}: {e}")
            failed_files.append(h5ad_file.name)
            continue

    if not all_adatas:
        raise ValueError("No valid AnnData objects were loaded")

    logger.info(f"\nSuccessfully loaded {len(all_adatas)} AnnData objects")
    if failed_files:
        logger.warning(f"Failed to load {len(failed_files)} files: {failed_files}")

    # Concatenate
    logger.info("\nConcatenating AnnData objects...")
    logger.info("Using join='outer' to keep all genes...")

    try:
        combined = ad.concat(all_adatas, join="outer", fill_value=0, merge="unique")

        # Verify .X is not None after concatenation
        if combined.X is None:
            logger.warning("Combined AnnData has None for .X with join='outer'")
            logger.info("Trying join='inner' instead...")
            combined = ad.concat(all_adatas, join="inner", merge="unique")

            if combined.X is None:
                raise ValueError(
                    "Concatenation failed: .X is None even with join='inner'"
                )

        logger.info(f"\nCombined AnnData shape: {combined.shape}")
        logger.info(f"Combined .X shape: {combined.X.shape}")
        logger.info(f"Number of cells: {combined.n_obs}")
        logger.info(f"Number of genes: {combined.n_vars}")

        # Check for unique obs_names
        if not combined.obs_names.is_unique:
            logger.warning("obs_names are not unique after concatenation!")
            logger.info("Making obs_names unique...")
            combined.obs_names_make_unique()

        # Save concatenated AnnData
        output_path = adata_dir / "all_samples.h5ad"

        logger.info(f"\nSaving concatenated AnnData to {output_path}...")
        combined.write(output_path)
        logger.info("Done!")

        return combined

    except Exception as e:
        logger.error(f"Concatenation failed: {e}")
        raise


def run_format(io_config: IOConfig):
    """Run the formatting module across multiple ROIs."""
    # Convert all Xenium data (this now also concatenates)
    convert_all_xenium(io_config)

    logger.info("Finished formatting all ROIs.")
