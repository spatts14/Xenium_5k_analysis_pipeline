# Import packages

import scanpy as sc
import spatialdata as sd

from helper_function.py import seed_everything

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

# Preform dimension reduction analysis
sc.pp.pca(adata)  # compute prinical components
sc.pp.neighbors(adata)  # compute a neighborhood graph
sc.tl.umap(adata)  # calculate umap
sc.tl.leiden(adata)  # determine cell clusters


# plot UMAP
# ^ Save plot
sc.pl.umap(
    adata,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "leiden",
    ],
    wspace=0.4,
)
