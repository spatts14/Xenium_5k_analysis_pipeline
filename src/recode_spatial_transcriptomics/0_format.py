"""Module for formatting Xenium data into Zarr format."""

# Import packages
import logging
import warnings  # ? what is the best way to suppress warnings from package inputs?

from spatialdata_io import xenium

from .paths import logging_path, xenium_path, zarr_path

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    # Set up logging
    logging.basicConfig(
        filename=logging_path / "0_format.txt",  # output file
        filemode="w",  # overwrites the file each time
        format="%(asctime)s - %(levelname)s - %(message)s",  # log format
        level=logging.INFO,  # minimum level to log
    )

    # Load into spatialdata format
    logging.info("Reading Xenium data...")
    sdata = xenium(xenium_path)

    # Convert to zarr format
    logging.info("Writing to Zarr...")
    sdata.write(zarr_path)

    # Convert to zarr format
    logging.info("Finished formatting data.")
