"""Dimension reduction module."""

import warnings
from logging import getLogger

import geosketch
import numpy as np
import scanpy as sc
import squidpy as sq

from recode_st.config import DimensionReductionModuleConfig, IOConfig
from recode_st.logging_config import configure_logging

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def run_dimension_reduction(
    config: DimensionReductionModuleConfig, io_config: IOConfig
):
    """Run dimension reduction on Xenium data."""
    # Set variables
    module_dir = io_config.output_dir / config.module_name
    n_comps = config.n_comps
    n_neighbors = config.n_neighbors
    resolution = config.resolution
    cluster_name = config.cluster_name
    norm_approach = (
        config.norm_approach if hasattr(config, "norm_approach") else "cell_area"
    )

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set the directory where to save the ScanPy figures
    sc.settings.figdir = module_dir

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(
        io_config.output_dir / "1_quality_control" / f"adata_{norm_approach}.h5ad"
    )

    # Set your target size
    n_total = 5000  # total cells you want in your dev set

    # Calculate proportional samples per ROI
    roi_counts = adata.obs["ROI"].value_counts()
    sampled_indices = []

    for roi in adata.obs["ROI"].unique():
        roi_mask = adata.obs["ROI"] == roi
        roi_data = adata[roi_mask]

        # Proportional to original ROI size
        n_roi_sketch = int(n_total * (roi_counts[roi] / len(adata)))
        n_roi_sketch = max(50, n_roi_sketch)  # ensure minimum samples per ROI
        n_roi_sketch = min(n_roi_sketch, len(roi_data))  # don't exceed available cells

        # Geometric sketch within this ROI
        if n_roi_sketch >= len(roi_data):
            # Take all cells if sketch size >= available cells
            roi_sketch_local = np.arange(len(roi_data))
        else:
            roi_sketch_local = geosketch.gs(roi_data.X, n_roi_sketch, replace=False)

        # Convert local indices to global indices
        global_indices = np.where(roi_mask)[0][roi_sketch_local]
        sampled_indices.extend(global_indices)

    adata = adata[
        sampled_indices, :
    ].copy()  # replace adata with the subsampled version for dev

    # Verify distribution
    print(adata.obs["ROI"].value_counts().sort_index())
    print(f"Total cells: {len(adata)}")

    # Highly variable genes
    logger.info("Selecting highly variable genes...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,  # select top highly variable genes
        # batch_key='run'
    )

    # Perform dimension reduction analysis
    logger.info("Compute PCA...")
    sc.pp.pca(adata, n_comps=50)  # Number of PC calculated compute principal components
    sc.pl.pca_variance_ratio(
        adata,
        log=True,
        n_pcs=50,  # Number of PCs shown in the plot
        show=False,
        save=f"_{config.module_name}.png",
    )
    logger.info(f"PCA Variance plot saved to {sc.settings.figdir}")

    logger.info("Compute neighbors...")
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,  # compute a neighborhood graph
        n_pcs=n_comps,  # For 5K panel, 30-50 PCs is typical
    )

    logger.info("Create UMAPs and cluster cells..")
    sc.tl.umap(adata)  # calculate umap
    sc.tl.leiden(
        adata,
        resolution=resolution,  # choose resolution for clustering
        key_added=cluster_name,
    )  # name clusters

    # plot UMAP
    logger.info("Plotting UMAPs...")
    sc.pl.umap(
        adata,
        color=[
            "total_counts",
            "n_genes_by_counts",
            cluster_name,
        ],
        wspace=0.4,
        show=False,
        save=f"_{config.module_name}.png",  # save the figure with the module name
        frameon=False,
    )
    logger.info(f"UMAP plot saved to {sc.settings.figdir}")

    # plot visualization of leiden clusters
    logger.info(f"Plotting {cluster_name} clusters...")
    for roi in adata.obs["ROI"].unique():
        subset = adata[adata.obs["ROI"] == roi]
        sq.pl.spatial_scatter(
            subset,
            library_id="spatial",
            shape=None,
            color=[cluster_name],
            wspace=0.4,
            save=module_dir / f"{cluster_name}_{roi}_spatial.png",
        )
    logger.info(f"{cluster_name} spatial scatter plot saved to {module_dir}")

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.2_dimension_reduction")

    # Set seed
    # seed_everything(21122023)

    try:
        run_dimension_reduction(
            DimensionReductionModuleConfig(module_name="2_dimension_reduction"),
            IOConfig(),
        )
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
