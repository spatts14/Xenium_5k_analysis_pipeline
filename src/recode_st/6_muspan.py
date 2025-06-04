"""Muspan module."""

# Import packages

import scanpy as sc
import spatialdata as sd

from .helper_function import seed_everything
from .paths import output_path, zarr_path

if __name__ == "__main__":
    # Set seed
    seed_everything(21122023)

    # Set variables

    # Import data
    sdata = sd.read_zarr(zarr_path)
    adata = adata = sc.read_h5ad(output_path / "data.h5ad")
