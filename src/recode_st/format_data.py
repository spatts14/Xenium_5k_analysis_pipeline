"""Module for formatting Xenium data into Zarr format."""

import warnings
from logging import getLogger
from pathlib import Path

from spatialdata_io import xenium

from recode_st.config import FormatDataModuleConfig
from recode_st.logging_config import configure_logging
from recode_st.paths import xenium_path, zarr_path

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def convert_xenium_to_zarr(xenium_path: Path, zarr_path: Path):
    """Convert Xenium data to Zarr format."""
    try:
        # Load Xenium data
        logger.info("Reading Xenium data...")
        sdata = xenium(xenium_path)
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
        raise err

    try:
        # Write to Zarr format
        logger.info("Writing to Zarr...")
        sdata.write(zarr_path, overwrite=True)
    except ValueError as err:
        logger.error(f"Failed writing to Zarr: {err}")
        raise err


def run_format(config: FormatDataModuleConfig):
    """Run the formatting module."""
    convert_xenium_to_zarr(xenium_path, zarr_path)

    logger.info("Finished formatting data.")


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.0_format")  # re-name the logger to match the module

    run_format(FormatDataModuleConfig(module_name="0_format"))
