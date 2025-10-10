"""Dimension reduction module."""

import warnings
from logging import getLogger

import geosketch
import numpy as np
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
    norm_approach = config.norm_approach
    subsample_data = config.subsample_data

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set the directory where to save the ScanPy figures
    sc.settings.figdir = module_dir

    # Import data
    logger.info(f"Loading Xenium data normalized with {norm_approach}...")
    adata = sc.read_h5ad(
        io_config.output_dir / "1_quality_control" / f"adata_{norm_approach}.h5ad"
    )

    #! Not sure if I need this
    # # Highly variable genes
    # logger.info("Selecting highly variable genes...")
    # sc.pp.highly_variable_genes(
    #     adata,
    #     n_top_genes=2000,  # select top highly variable genes
    #     # batch_key='run'
    # )

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

    if subsample_data is True:
        logger.info("Subsample data for dev...")
        orig_size = len(adata)
        logger.info(f"Original size: {orig_size} cells")

        # Set your target size
        n_total = 10000  # total cells you want in your dev set

        # Compute PCA if not already present
        if "X_pca" not in adata.obsm:
            logger.info("PCA not computed. Computing PCA...")
            sc.pp.pca(adata, n_comps=50)  # adjust n_comps as needed

        # Calculate proportional samples per ROI
        roi_counts = adata.obs["ROI"].value_counts()
        sampled_indices = []

        logger.info(f"Sketching {n_total} cells across {len(roi_counts)} ROIs...")

        for roi in adata.obs["ROI"].unique():
            roi_mask = adata.obs["ROI"] == roi
            roi_data = adata[roi_mask]

            # Proportional to original ROI size
            n_roi_sketch = int(n_total * (roi_counts[roi] / len(adata)))
            n_roi_sketch = max(50, n_roi_sketch)  # ensure minimum samples per ROI
            n_roi_sketch = min(
                n_roi_sketch, len(roi_data)
            )  # don't exceed available cells

            logger.info(
                f"  ROI {roi}: sketching {n_roi_sketch} from {len(roi_data)} cells"
            )

            # Geometric sketch within this ROI
            if n_roi_sketch >= len(roi_data):
                # Take all cells if sketch size >= available cells
                roi_sketch_local = np.arange(len(roi_data))
            else:
                # Use PCA representation (already dense)
                roi_sketch_local = geosketch.gs(
                    roi_data.obsm["X_pca"], n_roi_sketch, replace=False
                )

            # Convert local indices to global indices
            global_indices = np.where(roi_mask)[0][roi_sketch_local]
            sampled_indices.extend(global_indices)

        adata = adata[
            sampled_indices, :
        ].copy()  # replace adata with the subsampled version for dev

        adata.write_h5ad(module_dir / f"adata_sketch_{norm_approach}.h5ad")

        # Verify distribution
        logger.info("\n Sketching Results")
        logger.info(f"Sketched size: {len(adata)} cells")
        logger.info(f"Reduction: {100 * (1 - len(adata) / orig_size):.1f}%")
        logger.info("\nCells per ROI:")
        logger.info(adata.obs["ROI"].value_counts().sort_index())
    else:
        logger.info("Skipping sub-sampling data for dev...")

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
            "ROI",
            cluster_name,
        ],
        wspace=0.4,
        show=False,
        save=f"_{config.module_name}_{norm_approach}.png",  # save figure
        frameon=False,
    )

    sc.pl.umap(
        adata,
        color=[
            "PTPRC",  # immune cells
            "CD3E",  # T cells
            "CD68",  # macrophages
            "EPCAM",  # epithelial cells
            "COL1A1",  # collagen
            "PDGFRA",  # fibroblasts
            "ACTA2",  # smooth muscle cells
            "VWF",  # endothelial cells
        ],
        wspace=0.4,
        show=False,
        save=f"_{config.module_name}_cell_markers_{norm_approach}.png",  # save figure
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
    seed_everything(6)

    try:
        run_dimension_reduction(
            DimensionReductionModuleConfig(module_name="2_dimension_reduction"),
            IOConfig(),
        )
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
