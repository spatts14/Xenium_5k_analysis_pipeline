"""Data loading module for multi-ROI Xenium data."""

import warnings
from logging import getLogger
from pathlib import Path
from typing import List, Union
import spatialdata as sd
from anndata import AnnData, concat
from zarr.errors import PathNotFoundError
from recode_st.config import IOConfig
from recode_st.helper_function import seed_everything
from recode_st.logging_config import configure_logging

warnings.filterwarnings("ignore")
logger = getLogger(__name__)


def load_xenium_data(
    zarr_paths: Union[str, Path, List[Union[str, Path]]],
    batch_key: str = "roi_run"
) -> AnnData:
    """
    Load Xenium data from one or multiple ROI zarr stores.

    Parameters
    ----------
    zarr_paths : str, Path, or List[str/Path]
        Path(s) to the Xenium zarr store(s). Can be a single path or list of paths.
    batch_key : str, default="roi_run"
        Key to store batch information in adata.obs.
        Options: "roi", "run", or "roi_run".

    Returns
    -------
    AnnData
        Combined AnnData object with all ROIs.
    """
    # Convert single path to list
    if isinstance(zarr_paths, (str, Path)):
        zarr_paths = [zarr_paths]

    logger.info(f"Loading {len(zarr_paths)} Xenium ROI(s)...")

    adata_list = []

    for zarr_path in zarr_paths:
        zarr_path = Path(zarr_path)

        # Extract run and ROI names
        run_folder = zarr_path.parent.name
        run_name = run_folder.split("_")[-1] if "RUN" in run_folder else "unknown_run"
        parts = zarr_path.name.replace("output-", "").split("__")
        roi_name = next(p for p in parts if p.startswith(("COPD", "IPF", "PM08"))) # this is hard coded into for Sara's data, should be changed for other datasets

        try:
            logger.info(f"Loading ROI {roi_name} (Run {run_name}) from {zarr_path}")
            sdata = sd.read_zarr(zarr_path)
            adata = sdata.tables["table"]

            # Add metadata to obs
            adata.obs["roi"] = roi_name
            adata.obs["run"] = run_name
            adata.obs["roi_run"] = f"{run_name}_{roi_name}"

            # Assign batch_key column explicitly
            adata.obs[batch_key] = adata.obs[batch_key]

            # Store original zarr path in uns for reference
            if "roi_metadata" not in adata.uns:
                adata.uns["roi_metadata"] = {}
            adata.uns["roi_metadata"][adata.obs["roi_run"][0]] = str(zarr_path)

            # Spatial coordinates
            if "spatial" in adata.obsm:
                logger.info(f"  Spatial coordinates found for ROI {roi_name}")

            adata_list.append(adata)
            logger.info(f"  Loaded {adata.n_obs} cells, {adata.n_vars} genes")

        except PathNotFoundError as err:
            logger.error(f"File not found: {zarr_path}")
            raise err
        except KeyError as err:
            logger.error(f"'table' not found in zarr store: {zarr_path}")
            logger.error("Available tables: " + str(list(sdata.tables.keys()) if 'sdata' in locals() else "None"))
            raise err
        except Exception as e:
            logger.error(f"Error loading {zarr_path}: {str(e)}")
            raise e

    # Concatenate all ROIs
    if len(adata_list) > 1:
        logger.info("Concatenating ROIs...")
        adata_combined = concat(
            adata_list,
            join="outer",
            merge="same",
            uns_merge="first",
            label=batch_key,
            keys=[adata.obs[batch_key][0] for adata in adata_list],
            index_unique="_",
        )

        # Merge roi_metadata
        roi_metadata = {}
        for adata in adata_list:
            if "roi_metadata" in adata.uns:
                roi_metadata.update(adata.uns["roi_metadata"])
        adata_combined.uns["roi_metadata"] = roi_metadata

        logger.info(
            f"Combined: {adata_combined.n_obs} cells, {adata_combined.n_vars} genes"
        )
        logger.info(
            f"Batches ({batch_key}): {', '.join(adata_combined.obs[batch_key].unique())}"
        )
    else:
        adata_combined = adata_list[0]
        logger.info(
            f"Single ROI loaded: {adata_combined.n_obs} cells, {adata_combined.n_vars} genes"
        )

    logger.info("Data loading completed")
    return adata_combined


if __name__ == "__main__":
    # Configure logging
    configure_logging()

    # Load configuration
    config = IOConfig()

    # Set random seed for reproducibility
    seed_everything(config.seed)

    # Load Xenium data
    adata = load_xenium_data(
        zarr_paths=config.input.zarr_paths,
        batch_key=config.input.batch_key
    )

    # Save the combined AnnData object if output path is provided
    if config.output.adata_path:
        adata.write_h5ad(config.output.adata_path)
        logger.info(f"Combined AnnData saved to {config.output.adata_path}")   
