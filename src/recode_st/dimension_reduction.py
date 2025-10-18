"""Dimension reduction module."""

import warnings
from logging import getLogger
from pathlib import Path

import geosketch
import numpy as np
import scanpy as sc
import squidpy as sq

from recode_st.config import DimensionReductionModuleConfig, IOConfig

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def subsample_strategy_func(
    io_config: IOConfig,
    subsample_strategy: str,  # expected: "none", "compute", or "load"
    module_dir: Path,
    norm_approach: str,
    n_total: int = 10_000,
    min_cells_per_roi: int = 50,
    n_pca: int = 50,
) -> sc.AnnData:
    """Subsample or load a development dataset from an AnnData object.

    Args:
        io_config (IOConfig): IO configuration object.
        subsample_strategy (str): One of {"none", "compute", "load"}.
            - "compute": subsample and save the data
            - "load": load pre-subsampled data
            - "none": use full data
        module_dir (Path): Directory where the sketch file will be written or read.
        norm_approach (str): Label for the normalization approach, used in the filename.
        n_total (int, optional): Total num of cells to include in the subsampled data.
        Defaults to 10,000.
        min_cells_per_roi (int, optional): Minimum number of cells to sample per ROI.
        Defaults to 50.
        n_pca (int, optional): Number of PCA components to compute if missing.
        Defaults to 50.

    Returns:
        sc.AnnData: The subsampled (or loaded) AnnData object.
    """
    sketch_path = module_dir / f"adata_sketch_{norm_approach}.h5ad"

    if subsample_strategy == "load":
        logger.info("Loading subsampled data for dimension reduction...")
        if not sketch_path.exists():
            raise FileNotFoundError(
                f"Expected subsampled file not found: {sketch_path}\n"
                "Run with subsample_strategy='compute' first."
            )

        adata = sc.read_h5ad(sketch_path)
        logger.info(f"Loaded subsampled dataset with {len(adata)} cells.")

    elif subsample_strategy == "compute":
        logger.info("Computing subsampled data for dimension reduction...")
        logger.info(f"Loading full dataset normalized with {norm_approach}...")
        adata = sc.read_h5ad(
            io_config.output_dir / "1_quality_control" / f"adata_{norm_approach}.h5ad"
        )
        logger.info("Subsampling data for dev...")
        orig_size = len(adata)
        logger.info(f"Original size: {orig_size} cells")

        # Compute PCA if not already present
        if "X_pca" not in adata.obsm:
            logger.info("PCA not found. Computing PCA...")
            sc.pp.pca(adata, n_comps=60)  # compute 60 PCs

        roi_counts = adata.obs["ROI"].value_counts()
        sampled_indices = []
        logger.info(f"Sketching {n_total} cells across {len(roi_counts)} ROIs...")

        for roi in adata.obs["ROI"].unique():
            roi_mask = adata.obs["ROI"] == roi
            roi_data = adata[roi_mask]
            n_roi_sketch = int(n_total * (roi_counts[roi] / len(adata)))
            n_roi_sketch = max(min_cells_per_roi, n_roi_sketch)
            n_roi_sketch = min(n_roi_sketch, len(roi_data))

            logger.info(
                f"  ROI {roi}: sketching {n_roi_sketch} from {len(roi_data)} cells"
            )

            if n_roi_sketch >= len(roi_data):
                roi_sketch_local = np.arange(len(roi_data))
            else:
                roi_sketch_local = geosketch.gs(
                    roi_data.obsm["X_pca"], n_roi_sketch, replace=False
                )

            global_indices = np.where(roi_mask)[0][roi_sketch_local]
            sampled_indices.extend(global_indices)

        adata = adata[sampled_indices, :].copy()
        logger.info(f"Saved subsampled dataset: size {len(adata)} cells")
        adata.write_h5ad(sketch_path)

        logger.info("\nSketching Results")
        logger.info(f"Sketched size: {len(adata)} cells")
        logger.info(f"Reduction: {100 * (1 - len(adata) / orig_size):.1f}%")
        logger.info("\nCells per ROI:")
        logger.info(adata.obs["ROI"].value_counts().sort_index())

    elif subsample_strategy == "none":
        logger.info(f"Using full dataset normalized with {norm_approach}...")
        logger.info(f"Loading full dataset normalized with {norm_approach}...")
        adata = sc.read_h5ad(
            io_config.output_dir / "1_quality_control" / f"adata_{norm_approach}.h5ad"
        )
        if "X_pca" not in adata.obsm:
            logger.info("Computing PCA...")
            sc.pp.pca(adata, n_comps=60)  # compute 60 PCs
        logger.info("PCA computation complete.")

    else:
        raise ValueError(
            f"Invalid subsample_strategy: {subsample_strategy}. "
            "Choose from {'none', 'compute', 'load'}."
        )

    return adata


def run_dimension_reduction(
    config: DimensionReductionModuleConfig, io_config: IOConfig
):
    """Run dimension reduction on Xenium data."""
    # Set variables
    module_dir = io_config.output_dir / config.module_name
    n_pca = config.n_pca
    n_neighbors = config.n_neighbors
    resolution = config.resolution
    cluster_name = config.cluster_name
    norm_approach = config.norm_approach
    subsample_strategy = config.subsample_strategy

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set the directory where to save the ScanPy figures
    sc.settings.figdir = module_dir

    logger.info(f"Loading data using subsample strategy: {subsample_strategy}...")
    adata = subsample_strategy_func(
        io_config=io_config,
        subsample_strategy=subsample_strategy,
        module_dir=module_dir,
        norm_approach=norm_approach,
        n_total=10000,  # target ~10000 total cells
        min_cells_per_roi=100,  # ensure each ROI has at least 100
        n_pca=n_pca,  # compute 30 PCs if missing
    )

    logger.info("Compute neighbors...")
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,  # compute a neighborhood graph
        n_pcs=n_pca,  # For 5K panel, 30-50 PCs is typical
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
