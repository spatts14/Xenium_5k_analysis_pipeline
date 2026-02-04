"""Dimension reduction module."""

import warnings
from collections import Counter
from logging import getLogger
from pathlib import Path
from typing import Any, Literal

import geosketch
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq
from sklearn.metrics import silhouette_samples, silhouette_score

from recode_st.config import DimensionReductionModuleConfig, IOConfig
from recode_st.helper_function import configure_scanpy_figures

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)

SubsampleStrategy = Literal["none", "compute", "load"]


def create_subsample(
    adata: sc.AnnData,
    n_total: int = 10_000,
    min_cells_per_roi: int = 50,
    n_pca: int = 50,
) -> sc.AnnData:
    """Create a subsampled dataset using geosketch.

    Args:
        adata: Input AnnData object (full dataset)
        n_total: Total number of cells to subsample
        min_cells_per_roi: Minimum cells to keep per ROI
        n_pca: Number of PCA components for sketching

    Returns:
        Subsampled AnnData object
    """
    logger.info("Subsampling data...")
    orig_size = len(adata)
    logger.info(f"Original size: {orig_size} cells")

    # Ensure PCA is computed
    if "X_pca" not in adata.obsm:
        logger.info("Computing PCA for sketching...")
        sc.pp.pca(adata, n_comps=n_pca)

    roi_counts = adata.obs["ROI"].value_counts()
    sampled_indices = []
    logger.info(f"Sketching {n_total} cells across {len(roi_counts)} ROIs...")

    for roi in adata.obs["ROI"].unique():
        roi_mask = adata.obs["ROI"] == roi
        roi_data = adata[roi_mask]

        # Calculate proportional sample size for this ROI
        n_roi_sketch = int(n_total * (roi_counts[roi] / len(adata)))
        n_roi_sketch = max(min_cells_per_roi, n_roi_sketch)
        n_roi_sketch = min(n_roi_sketch, len(roi_data))

        logger.info(f"  ROI {roi}: sketching {n_roi_sketch} from {len(roi_data)} cells")

        # Sketch or take all
        if n_roi_sketch >= len(roi_data):
            roi_sketch_local = np.arange(len(roi_data))
        else:
            roi_sketch_local = geosketch.gs(
                roi_data.obsm["X_pca"], n_roi_sketch, replace=False
            )

        global_indices = np.where(roi_mask)[0][roi_sketch_local]
        sampled_indices.extend(global_indices)

    adata_sub = adata[sampled_indices, :].copy()

    # Log results
    logger.info("\nSubsampling Results:")
    logger.info(f"  Original: {orig_size} cells")
    logger.info(f"  Subsampled: {len(adata_sub)} cells")
    logger.info(f"  Reduction: {100 * (1 - len(adata_sub) / orig_size):.1f}%")
    logger.info("\nCells per ROI:")
    logger.info(adata_sub.obs["ROI"].value_counts().sort_index())

    return adata_sub


def load_data(
    io_config: IOConfig,
    config: DimensionReductionModuleConfig,
    norm_approach: str,
    subsample_strategy: SubsampleStrategy = "none",  # make default "none"
) -> sc.AnnData:
    """Load data - if subsample, load sub-sampled data.

    Args:
        io_config: IO configuration
        config: Dimension reduction module configuration
        norm_approach: Normalization approach label
        subsample_strategy: One of "none", "compute", or "load"
    Returns:
        AnnData object ready for analysis

    Raises:
        FileNotFoundError: If load strategy is used but file doesn't exist
        ValueError: If invalid strategy is provided
    """
    # Validate subsample_strategy
    valid_strategies = {"none", "compute", "load"}
    if subsample_strategy not in valid_strategies:
        raise ValueError(
            f"Invalid subsample_strategy '{subsample_strategy}'. "
            f"Must be one of {valid_strategies}"
        )

    # Define file paths
    data_path = (
        Path(io_config.output_dir) / "quality_control" / f"adata_{norm_approach}.h5ad"
    )
    module_dir = Path(io_config.output_dir) / config.module_name
    subsample_path = module_dir / f"adata_subsample_{config.norm_approach}.h5ad"

    if subsample_strategy == "none":
        logger.info("Loading full dataset (no subsampling)...")
        adata = sc.read_h5ad(data_path)

        # Compute PCA if needed
        if "X_pca" not in adata.obsm:
            logger.info("Computing PCA...")
            sc.pp.pca(adata, n_comps=config.n_pca)

        logger.info(f"Loaded {len(adata)} cells")
        return adata

    elif subsample_strategy == "compute":
        logger.info("Creating new subsample...")

        # Load full dataset
        adata = sc.read_h5ad(data_path)

        # Create subsample
        adata_sub = create_subsample(
            adata,
            n_total=config.subsample_n_total,
            min_cells_per_roi=config.subsample_min_cells_per_roi,
            n_pca=config.n_pca,
        )
        # Save subsample
        logger.info(f"Saving subsample to {subsample_path}")
        adata_sub.write_h5ad(subsample_path)

        return adata_sub

    elif subsample_strategy == "load":
        logger.info("Loading existing subsample...")

        if not subsample_path.exists():
            raise FileNotFoundError(
                f"Subsample file not found: {subsample_path}\n"
                f"Run with subsample_strategy='compute' first."
            )

        adata = sc.read_h5ad(subsample_path)
        logger.info(f"Loaded {len(adata)} cells")
        logger.info("Cells per ROI:")
        logger.info(adata.obs["ROI"].value_counts().sort_index())

        return adata

    else:
        raise ValueError(
            f"Invalid subsample_strategy: '{subsample_strategy}'. "
            f"Must be one of: 'none', 'compute', 'load'"
        )


def compute_dimensionality_reduction(
    adata: sc.AnnData,
    n_pca: int = 30,  # number of PCA components to use to compute neighbors
    n_neighbors: int = 10,  # number of neighbors for graph
    min_dist: float = 0.5,  # minimum distance for UMAP
) -> sc.AnnData:
    """Compute PCA, neighbors, UMAP, and clustering.

    Args:
        adata: Input AnnData object
        n_pca: Number of PCA components to use
        n_neighbors: Number of neighbors for graph
        resolution: Resolution for Leiden clustering
        cluster_name: Name for cluster annotation
        min_dist: Minimum distance for UMAP

    Returns:
        Nothing. Saves figures and AnnData with computed dimensionality reduction.
    """
    logger.info(
        f"Computing neighbors and UMAP with n_pca={n_pca}, n_neighbors={n_neighbors}..."
    )
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pca)

    logger.info("Computing UMAP...")
    sc.tl.umap(adata, min_dist=min_dist, spread=2.0)

    return adata


def calculate_clusters(
    adata: sc.AnnData,
    res_list: list[float] = [
        0.1,
        0.3,
        0.5,
        0.8,
        1.0,
    ],
) -> sc.AnnData:
    """Compute Leiden clustering.

    Args:
        adata: Input AnnData object
        res_list: Resolution for Leiden clustering

    Returns:
        AnnData with computed clustering.
    """
    flavor = "igraph"
    logger.info(
        f"Using the '{flavor}' implementation of Leiden for faster performance."
    )

    for res in res_list:
        logger.info(f"Clustering with resolution={res}...")

        # Cluster at this resolution
        key = f"leiden_res_{res}"
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=key,
            flavor=flavor,  # Faster than default leidenalg
            n_iterations=2,  # Fewer iterations, usually sufficient
        )
        logger.info(
            f"Leiden clustering completed successfully.Added '{key} to adata.obs."
        )

        # Verify the clustering was added
        if key not in adata.obs.columns:
            raise ValueError(f"Clustering column '{key}' was not created successfully.")

        logger.info(f"Number of clusters found: {len(adata.obs[key].unique())}")
    return adata


def evaluate_resolutions(
    adata: sc.AnnData,
    res_list: list[float] = [
        0.1,
        0.3,
        0.5,
        0.8,
        1.0,
    ],
    n_pcs: int = 30,
    use_rep: str = "X_pca",
) -> pd.DataFrame:
    """Evaluate pre-computed clustering with silhouette scores.

    Note: This function assumes clustering has already been computed using
    calculate_clusters(). It does NOT recompute clusters.

    Args:
        adata: AnnData object with 'leiden_res_{resolution}' columns
        res_list: List of resolution values to evaluate
        n_pcs: Number of PCs to use for silhouette calculation
        use_rep: Representation to use for silhouette score

    Returns:
        DataFrame with resolution metrics (n_clusters, silhouette scores, cluster sizes)
    """
    logger.info(f"Evaluating {len(res_list)} resolutions...")

    # Validate that clustering exists
    missing = [res for res in res_list if f"leiden_res_{res}" not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"Clustering not found for resolutions: {missing}. "
            f"Run calculate_clusters() first."
        )

    # Get PCA representation
    if use_rep not in adata.obsm:
        raise ValueError(f"Representation '{use_rep}' not found in adata.obsm.")
    X = adata.obsm[use_rep][:, :n_pcs]

    results = []

    for res in res_list:
        key = f"leiden_res_{res}"
        labels = adata.obs[key].astype(int).values
        n_clusters = len(np.unique(labels))

        logger.info(f"  Evaluating resolution {res} ({n_clusters} clusters)...")

        # Calculate silhouette scores (skip if only 1 cluster)
        if n_clusters > 1:
            sil_score = silhouette_score(X, labels)
            sil_samples = silhouette_samples(X, labels)

            # Per-cluster silhouette scores
            cluster_sil_scores = {
                cid: sil_samples[labels == cid].mean() for cid in np.unique(labels)
            }
            min_cluster_sil = min(cluster_sil_scores.values())
            max_cluster_sil = max(cluster_sil_scores.values())
        else:
            sil_score = np.nan
            min_cluster_sil = np.nan
            max_cluster_sil = np.nan

        # Cluster sizes
        cluster_sizes = Counter(labels)

        results.append(
            {
                "resolution": res,
                "n_clusters": n_clusters,
                "silhouette_score": sil_score,
                "min_cluster_silhouette": min_cluster_sil,
                "max_cluster_silhouette": max_cluster_sil,
                "min_cluster_size": min(cluster_sizes.values()),
                "max_cluster_size": max(cluster_sizes.values()),
                "cluster_key": key,
            }
        )

        logger.info(f"Silhouette={sil_score:.3f}")

    metrics_df = pd.DataFrame(results)

    # Log summary
    best_idx = metrics_df["silhouette_score"].idxmax()
    best = metrics_df.loc[best_idx]
    logger.info(
        f"Best resolution by silhouette: {best['resolution']} "
        f"({best['n_clusters']} clusters, silhouette={best['silhouette_score']:.3f})"
    )

    return metrics_df


def plot_dimensionality_reduction(
    adata: sc.AnnData,
    norm_approach: str,
    module_name: str,
    n_neighbors: int,
    figdir: Path,
    cmap: Any = sns.color_palette("crest", as_cmap=True),
    cluster_name: str = "leiden",
    config: DimensionReductionModuleConfig = None,
) -> None:
    """Create and save dimensionality reduction plots.

    Args:
        adata: AnnData with computed UMAP and clusters
        norm_approach: Normalization approach label
        module_name: Name of module for file naming
        n_neighbors: Number of neighbors used
        figdir: Directory to save figures
        cmap: Colormap for plots
        cluster_name: Name of cluster annotation
        config: Configuration object containing visualization settings
    """
    logger.info("Plotting UMAPs...")

    # QC and metadata plots
    logger.info("Plotting UMAPs with QC meterics...")

    # Check if cluster column exists before plotting
    color_list = ["total_counts", "n_genes_by_counts"]
    if cluster_name in adata.obs.columns:
        color_list.append(cluster_name)  # adding cluster_name to color list
        logger.info(f"Including '{cluster_name}' in UMAP plots.")
    else:
        logger.warning(
            f"Cluster column '{cluster_name}' not found in adata.obs."
            f"Skipping cluster visualization."
        )
        logger.info(f"Available columns in adata.obs: {list(adata.obs.columns)}")

    sc.pl.umap(
        adata,
        color=color_list,
        ncols=3,
        cmap=cmap,
        wspace=0.4,
        show=False,
        save=f"_{module_name}_{norm_approach}_neighbors_{n_neighbors}_{config.resolution}.pdf",
        frameon=False,
    )

    # Plot observation fields if available
    if config and config.obs_vis_list:
        logger.info("Plotting UMAPs with observation fields...")
        sc.pl.umap(
            adata,
            color=config.obs_vis_list,
            ncols=3,
            cmap=cmap,
            wspace=0.4,
            show=False,
            save=f"_{module_name}_{norm_approach}_neighbors_{n_neighbors}_{config.resolution}_obs_fields.pdf",
            frameon=False,
        )
    else:
        logger.info("Skipping observation fields plot (obs_vis_list not configured)")

    # Plot marker genes if available
    if config and config.marker_genes:
        logger.info("Plotting UMAPs with marker genes...")
        sc.pl.umap(
            adata,
            color=config.marker_genes,
            cmap=cmap,
            ncols=4,
            wspace=0.4,
            show=False,
            save=f"_{module_name}_cell_markers_{norm_approach}_neighbors_{n_neighbors}_{config.resolution}.pdf",
            frameon=False,
        )
    else:
        logger.info("Skipping marker genes plot (marker_genes not configured)")

    logger.info(f"UMAP plots saved to {figdir}")


def plot_spatial_distribution(
    adata: sc.AnnData,
    cluster_name: str,
    module_dir: Path,
) -> None:
    """Plot spatial distribution of clusters for each ROI.

    Args:
        adata: AnnData with cluster annotations
        cluster_name: Name of cluster annotation
        module_dir: Directory to save plots
    """
    # Check if cluster column exists before plotting
    if cluster_name not in adata.obs.columns:
        logger.warning(
            f"Cluster column '{cluster_name}' not found in adata.obs."
            f"Skipping spatial cluster plots."
        )
        logger.info(f"Available columns in adata.obs: {list(adata.obs.columns)}")
        return

    logger.info(f"Plotting {cluster_name} spatial distribution...")

    for roi in adata.obs["ROI"].unique():
        subset = adata[adata.obs["ROI"] == roi]
        sq.pl.spatial_scatter(
            subset,
            library_id="spatial",
            shape=None,
            color=[cluster_name],
            wspace=0.4,
            save=module_dir / f"{cluster_name}_{roi}_spatial.pdf",
        )

    logger.info(f"Spatial plots saved to {module_dir}")


def run_dimension_reduction(
    config: DimensionReductionModuleConfig, io_config: IOConfig
) -> sc.AnnData:
    """Run dimension reduction on Xenium data.

    Args:
        config: Configuration for dimension reduction
        io_config: IO configuration

    Returns:
        AnnData object with computed dimensionality reduction
    """
    # Setup
    module_dir = io_config.output_dir / config.module_name
    module_dir.mkdir(exist_ok=True, parents=True)
    spatial_plots_dir = module_dir / "spatial_plots"
    spatial_plots_dir.mkdir(exist_ok=True, parents=True)

    # Set figure directory for this module (overrides global setting)
    sc.settings.figdir = module_dir

    # Set figure settings to ensure consistency across all modules
    configure_scanpy_figures(str(io_config.output_dir))
    cmap = sns.color_palette("Spectral", as_cmap=True)
    palette = sns.color_palette("hls", 16)

    adata = load_data(
        io_config=io_config,
        config=config,
        norm_approach=config.norm_approach,
        subsample_strategy=config.subsample_strategy,
    )

    # Plot PCA variance
    logger.info("Computing and plotting PCA variance ratio...")
    sc.pp.pca(adata, n_comps=80, svd_solver="arpack", use_highly_variable=True)

    logger.info("Plotting PCA variance...")
    sc.pl.pca_variance_ratio(
        adata,
        log=True,
        n_pcs=80,
        show=False,
        save=f"_{config.module_name}.pdf",
    )

    logger.info("Plotting PCA scatter plots...")
    sc.pl.pca(
        adata,
        color=["total_counts", "n_genes_by_counts"],
        cmap=cmap,
        dimensions=[(0, 1)],
        ncols=2,
        size=2,
        show=False,
        save=f"_{config.module_name}.pdf",
    )

    # Plot PCA with observation fields if available
    if config.obs_vis_list:
        logger.info("Plotting PCA with observation fields...")
        for _obs_field in config.obs_vis_list:
            sc.pl.pca(
                adata,
                color=_obs_field,
                palette=palette,
                dimensions=[(0, 1)],
                ncols=3,
                size=2,
                alpha=0.5,
                show=False,
                save=f"_{config.module_name}_{_obs_field}_PC1_PC2.pdf",
            )
        for _obs_field in config.obs_vis_list:
            sc.pl.pca(
                adata,
                color=_obs_field,
                palette=palette,
                dimensions=[(2, 3)],
                ncols=3,
                size=2,
                alpha=0.5,
                show=False,
                save=f"_{config.module_name}_{_obs_field}_PC3_PC4.pdf",
            )
    else:
        logger.info(
            "Skipping PCA observation fields plot (obs_vis_list not configured)"
        )

    # Compute dimensionality reduction
    logger.info("Computing dimension reduction...")
    adata = compute_dimensionality_reduction(
        adata=adata,
        n_pca=config.n_pca,
        n_neighbors=config.n_neighbors,
        min_dist=0.1,
    )

    # Compute clustering
    logger.info("Computing clustering...")
    adata = calculate_clusters(adata=adata, res_list=[config.resolution])

    # Evaluate clustering quality (no re-clustering)
    metrics_df = evaluate_resolutions(
        adata, res_list=[config.resolution], n_pcs=config.n_pca
    )

    # Set best clustering as main cluster annotation
    best_res = metrics_df.loc[metrics_df["silhouette_score"].idxmax(), "resolution"]
    logger.info(f"Setting leiden_res_{best_res} as primary cluster annotation...")
    adata.obs["leiden_best"] = adata.obs[f"leiden_res_{best_res}"]

    logger.info("Saving adata with computed dimension reduction...")
    output_path = module_dir / "adata.h5ad"
    adata.write_h5ad(output_path)
    logger.info(f"Final results saved to {output_path}")

    logger.info("Plotting dimension reduction...")

    # Plot results
    plot_dimensionality_reduction(
        adata=adata,
        cluster_name="leiden_best",
        norm_approach=config.norm_approach,
        module_name=config.module_name,
        n_neighbors=config.n_neighbors,
        cmap=cmap,
        figdir=module_dir,
        config=config,
    )

    logger.info("Plotting spatial distribution of clusters...")
    plot_spatial_distribution(
        adata=adata,
        cluster_name="leiden_best",
        module_dir=spatial_plots_dir,
    )

    logger.info(f"Dimension reduction module '{config.module_name}' complete.")
