"""Dimension reduction module."""

import warnings
from logging import getLogger

import scanpy as sc
import squidpy as sq

from recode_st.config import DimensionReductionModuleConfig, IOConfig
from recode_st.helper_function import seed_everything
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

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set the directory where to save the ScanPy figures
    sc.settings.figdir = module_dir

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "1_qc" / "adata.h5ad")

    # Perform dimension reduction analysis
    logger.info("Compute PCA...")
    sc.pp.pca(adata, n_comps=n_comps)  # compute principal components
    sc.pl.pca_variance_ratio(
        adata,
        log=True,
        n_pcs=50,
        show=False,
        save=f"_{config.module_name}.png",
    )
    logger.info(f"PCA Variance plot saved to {sc.settings.figdir}")

    logger.info("Compute neighbors...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)  # compute a neighborhood graph

    logger.info("Create UMAPs and cluster cells..")
    sc.tl.umap(adata)  # calculate umap
    sc.tl.leiden(
        adata,
        resolution=resolution,  # choose resolution for clustering
        key_added=cluster_name,
    )  # name leiden clusters

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
    sq.pl.spatial_scatter(
        adata,
        library_id="spatial",
        shape=None,
        color=[cluster_name],
        wspace=0.4,
        save=module_dir / f"{cluster_name}_spatial.png",
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
    seed_everything(21122023)

    try:
        run_dimension_reduction(
            DimensionReductionModuleConfig(module_name="2_dimension_reduction"),
            IOConfig(),
        )
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
