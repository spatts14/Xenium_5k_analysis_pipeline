# Import packages
import numpy as np

import spatialdata as sd
from spatialdata_io import xenium

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

import random
import torch
import os

from helper_function.py import seed_everything

# Set seed
seed_everything(21122023)

# Set variables

# Set directories
input_path = "./"
output_path = "./"
xenium_path = f"{input_path}, /Xenium" # ^ Change to file path rather than f" string
zarr_path = f"{output_path}, /Xenium.zarr" # ^ Change to file path rather than f" string

# Import data
sdata = sd.read_zarr(zarr_path)
adata =  adata = sc.read_h5ad(f"{output_path}/data.h5ad")

# View plots
# ^ View and save each plot
sq.pl.spatial_scatter(
    adata,
    library_id="spatial",
    shape=None,
    color=[
        "leiden",
    ],
    wspace=0.4,
)