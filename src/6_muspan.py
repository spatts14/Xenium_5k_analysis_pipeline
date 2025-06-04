"""Muspan module."""

# Import packages
import logging
import os
from pathlib import Path

import scanpy as sc
import spatialdata as sd

from .paths import base_dir, logging_path, output_path, zarr_path

# Set variables
module_name = "6_muspan"  # name of the module

# Confirm directories exist
if not Path(base_dir).exists():
    raise FileNotFoundError(f"Input path {base_dir} does not exist.")
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
