"""Muspan module."""

import warnings
from logging import getLogger

import scanpy as sc
import spatialdata as sd

from recode_st.helper_function import seed_everything
from recode_st.logging_config import configure_logging
from recode_st.paths import output_path, zarr_path

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def run_muspan():
    """Run Muspan on Xenium data."""
    try:
        import muspan  # noqa: F401
    except ModuleNotFoundError as err:
        logger.error(
            "Could not load necessary MuSpAn package. You can obtain this with:\n"
            "    pip install 'recode_st[muspan] @ git+"
            "https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics.git"
        )
        raise err
    # Set variables
    module_name = "6_muspan"  # name of the module
    module_dir = output_path / module_name
    seed = 21122023  # seed for reproducibility

    # Set seed
    seed_everything(seed)

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(output_path / "5_spatial_stats" / "adata.h5ad")  # noqa: F841
    sdata = sd.read_zarr(zarr_path)  # noqa: F841


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.6_muspan")  # re-name the logger to match the module

    try:
        run_muspan()
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
