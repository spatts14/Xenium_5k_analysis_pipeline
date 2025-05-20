# Import packages

import scanpy as sc
import spatialdata as sd
import squidpy as sq

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

# $ Calculate spatial statistics

# build spatial neighborhood graph
# compute connectivity
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)


# compute centrality scores

# closeness centrality - measure of how close the group is to other nodes (i.e cells)
# clustering coefficient - measure of the degree to which nodes (cells) cluster together.
# degree centrality - fraction of non-group members connected to group members

sq.gr.centrality_scores(adata, cluster_key="leiden")

# plot centrality scores
sq.pl.centrality_scores(adata, cluster_key="leiden", figsize=(16, 5))


# Compute co-occurrence probability

# Create subset table layer
sdata.tables["subsample"] = sc.pp.subsample(adata, fraction=0.5, copy=True)
adata_subsample = sdata.tables["subsample"]

# Visualize co-occurrence
sq.gr.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
)
sq.pl.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
    clusters="12",
    figsize=(10, 10),
)
sq.pl.spatial_scatter(
    adata_subsample,
    color="leiden",
    shape=None,
    size=2,
)
