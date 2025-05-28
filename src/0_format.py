# Import packages
import warnings  # ? what is the best way to suppress warnings from package inputs?

warnings.filterwarnings("ignore")

import spatialdata as sd
from spatialdata_io import xenium
import scanpy as sc
import squidpy as sq

from pathlib import Path

import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler("my_script.log"), logging.StreamHandler()],
)

# Set directories -
input_path = "/Users/sarapatti/Desktop/PhD_projects/Llyod_lab/ReCoDe-spatial-transcriptomics"
output_path = "/Users/sarapatti/Desktop/PhD_projects/Llyod_lab/ReCoDe-spatial-transcriptomics/analysis"
xenium_path = Path(input_path) / "data/xenium"
zarr_path = Path(input_path) / "data/xenium.zarr"

# Load into spatialdata format
logging.info("Reading Xenium data...")
sdata = xenium(xenium_path)

# Convert to zarr format
logging.info("Writing to Zarr...")
sdata.write(zarr_path)

# Convert to zarr format
logging.info("Finished formatting data.")
