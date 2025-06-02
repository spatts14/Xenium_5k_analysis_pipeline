# Import packages
import logging
import warnings  # ? what is the best way to suppress warnings from package inputs?
from pathlib import Path

from spatialdata_io import xenium

warnings.filterwarnings("ignore")

# Set directories
input_path = (
    "/Users/sarapatti/Desktop/PhD_projects/Llyod_lab/ReCoDe-spatial-transcriptomics"
)
output_path = "/Users/sarapatti/Desktop/PhD_projects/Llyod_lab/ReCoDe-spatial-transcriptomics/analysis"
xenium_path = Path(input_path) / "data/xenium"
zarr_path = Path(input_path) / "data/xenium.zarr"
logging_path = "/Users/sarapatti/Desktop/PhD_projects/Llyod_lab/ReCoDe-spatial-transcriptomics/analysis/logging"

# Set up logging
logging.basicConfig(
    filename=Path(logging_path) / "0_format.txt",  # output file
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
