# Import packages

import scanpy as sc
import spatialdata as sd

# Set seed
seed_everything(21122023)

# Set variables

# Set directories
input_path = "./"
output_path = "./"
xenium_path = f"{input_path}, /Xenium"  # ^ Change to file path rather than f" string
zarr_path = (
    f"{output_path}, /Xenium.zarr"  # ^ Change to file path rather than f" string
)

# Import data
sdata = sd.read_zarr(zarr_path)
adata = adata = sc.read_h5ad(f"{output_path}/data.h5ad")
