# Import packages
import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import muspan as ms
import numpy as np
import scanpy as sc
import spatialdata as sd
import squidpy as sq

# Set variables
module_name = "6_muspan"  # name of the module

# Set directories
base_dir = "/Users/sarapatti/Desktop/PhD_projects/Llyod_lab/ReCoDe-spatial-transcriptomics"
input_path = base_dir
output_path = Path(base_dir) / "analysis"
logging_path = Path(output_path) / "logging"
zarr_path = Path(input_path) / "data/xenium.zarr"

# Confirm directories exist
if not Path(input_path).exists():
    raise FileNotFoundError(f"Input path {input_path} does not exist.")
if not Path(output_path).exists():
    raise FileNotFoundError(f"Output path {output_path} does not exist.")


# Create output directories if they do not exist
os.makedirs(Path(output_path) / module_name, exist_ok=True)

# Set up logging
os.makedirs(
    logging_path, exist_ok=True
)  # should set up all these directories at the start of the pipeline?
logging.basicConfig(
    filename=Path(logging_path) / f"{module_name}.txt",  # output file
    filemode="w",  # overwrites the file each time
    format="%(asctime)s - %(levelname)s - %(message)s",  # log format
    level=logging.INFO,  # minimum level to log
)

# change directory to output_path/module_name
os.chdir(
    Path(output_path) / module_name
)  # need to so plots save in the correct directory

# Import data
logging.info("Loading Xenium data...")
adata = sc.read_h5ad(Path(output_path) / "4_view_images/adata.h5ad")
sdata = sd.read_zarr(zarr_path)
