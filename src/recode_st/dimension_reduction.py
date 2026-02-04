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

# Default resolution list - single source of truth
DEFAULT_RESOLUTIONS: list[float] = [0.1, 0.3, 0.5, 0.8, 1.0]


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

        logger.info(f" ROI {roi}: sketching {n_roi_sketch} from {len(roi_data)} cells")

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
    logger.info(" Subsampling Results:")
    logger.info(f"Original: {orig_size} cells")
    logger.info(f"Subsampled: {len(adata_sub)} cells")
    logger.info(f"Reduction: {100 * (1 - len(adata_sub) / orig_size):.1f}%")

    return adata_sub


def load_data(
    io_config: IOConfig,
    config: DimensionReductionModuleConfig,
    subsample_strategy: SubsampleStrategy = "none",
) -> sc.AnnData:
    """Load data with optional subsampling.

    Args:
        io_config: IO configuration
        config: Dimension reduction module configuration
        subsample_strategy: One of "none", "compute", or "load"

    Returns:
        AnnData object ready for analysis

    Raises:
        FileNotFoundError: If load strategy is used but file doesn't exist
        ValueError: If invalid strategy is provided
    """
    valid_strategies = {"none", "compute", "load"}
    if subsample_strategy not in valid_strategies:
        raise ValueError(
            f"Invalid subsample_strategy '{subsample_strategy}'. "
            f"Must be one of {valid_strategies}"
        )

    # Define file paths
    data_path = (
        Path(io_config.output_dir)
        / "quality_control"
        / f"adata_{config.norm_approach}.h5ad"
    )
    module_dir = Path(io_config.output_dir) / config.module_name
    subsample_path = module_dir / f"adata_subsample_{config.norm_approach}.h5ad"

    if subsample_strategy == "none":
        logger.info("Loading full dataset (no subsampling)...")
        adata = sc.read_h5ad(data_path)
        logger.info(f"  Loaded {len(adata)} cells")
        return adata

    elif subsample_strategy == "compute":
        logger.info("Creating new subsample...")
        adata = sc.read_h5ad(data_path)
        adata_sub = create_subsample(
            adata,
            n_total=config.subsample_n_total,
            min_cells_per_roi=config.subsample_min_cells_per_roi,
            n_pca=config.n_pca,
        )
        logger.info(f"  Saving subsample to {subsample_path}")
        adata_sub.write_h5ad(subsample_path)
        return adata_sub

    else:  # "load"
        logger.info("Loading existing subsample...")
        if not subsample_path.exists():
            raise FileNotFoundError(
                f"Subsample file not found: {subsample_path}\n"
                f"Run with subsample_strategy='compute' first."
            )
        adata = sc.read_h5ad(subsample_path)
        logger.info(f"Loaded {len(adata)} cells")
        return adata


def compute_dimensionality_reduction(
    adata: sc.AnnData,
    n_pca: int = 30,
    n_neighbors: int = 10,
    min_dist: float = 0.1,
    spread: float = 2.0,
) -> sc.AnnData:
    """Compute PCA, neighbors, and UMAP.

    Args:
        adata: Input AnnData object
        n_pca: Number of PCA components to use for neighbors
        n_neighbors: Number of neighbors for graph construction
        min_dist: Minimum distance for UMAP
        spread: Spread parameter for UMAP

    Returns:
        AnnData with computed neighbors and UMAP embeddings.
    """
    logger.info(f"Computing neighbors (n_pcs={n_pca}, n_neighbors={n_neighbors})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pca)

    logger.info(f"Computing UMAP (min_dist={min_dist}, spread={spread})...")
    sc.tl.umap(adata, min_dist=min_dist, spread=spread)

    return adata


def calculate_clusters(
    adata: sc.AnnData,
    res_list: list[float] = DEFAULT_RESOLUTIONS,
    flavor: str = "igraph",
    n_iterations: int = 2,
) -> sc.AnnData:
    """Compute Leiden clustering at multiple resolutions.

    Args:
        adata: Input AnnData object (must have neighbors computed)
        res_list: List of resolution values for clustering
        flavor: Leiden implementation ('igraph' or 'leidenalg')
        n_iterations: Number of Leiden iterations
        random_state: Random state for reproducibility

    Returns:
        AnnData with 'leiden_res_{resolution}' columns added.

    Raises:
        ValueError: If neighbors have not been computed.
    """
    if res_list is None:
        res_list = DEFAULT_RESOLUTIONS

    if "neighbors" not in adata.uns:
        raise ValueError(
            "Neighbors not computed. Run compute_dimensionality_reduction() first."
        )

    logger.info(f"Computing Leiden clustering at {len(res_list)} resolutions...")
    logger.info(f"  Using '{flavor}' implementation")

    for res in res_list:
        key = f"leiden_res_{res}"
        logger.info(f"Resolution {res}...")

        try:
            sc.tl.leiden(
                adata,
                resolution=res,
                key_added=key,
                flavor=flavor,
                n_iterations=n_iterations,
            )
        except Exception as e:
            logger.warning(f"{flavor} failed: {e}. Falling back to leidenalg...")
            sc.tl.leiden(
                adata,
                resolution=res,
                key_added=key,
                flavor="leidenalg",
                n_iterations=n_iterations,
            )

        if key not in adata.obs.columns:
            raise ValueError(f"Clustering column '{key}' was not created.")

        n_clusters = adata.obs[key].nunique()
        logger.info(f"For {key} created {n_clusters} clusters")

    return adata


def evaluate_resolutions(
    adata: sc.AnnData,
    res_list: list[float] = DEFAULT_RESOLUTIONS,
    n_pcs: int = 30,
    use_rep: str = "X_pca",
) -> pd.DataFrame:
    """Evaluate pre-computed clustering with silhouette scores.

    Args:
        adata: AnnData with 'leiden_res_{resolution}' columns
        res_list: List of resolution values to evaluate
        n_pcs: Number of PCs for silhouette calculation
        use_rep: Representation to use for silhouette score

    Returns:
        DataFrame with metrics: resolution, n_clusters, silhouette scores, cluster sizes

    Raises:
        ValueError: If clustering columns are missing or representation not found.
    """
    if res_list is None:
        res_list = DEFAULT_RESOLUTIONS

    logger.info(f"Evaluating {len(res_list)} resolutions...")

    # Validate clustering exists
    missing = [res for res in res_list if f"leiden_res_{res}" not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"Clustering not found for resolutions: {missing}. "
            f"Run calculate_clusters() first."
        )

    # Validate representation exists
    if use_rep not in adata.obsm:
        raise ValueError(f"Representation '{use_rep}' not found in adata.obsm.")

    X = adata.obsm[use_rep][:, :n_pcs]
    results = []

    for res in res_list:
        key = f"leiden_res_{res}"
        labels = adata.obs[key].astype(int).values
        n_clusters = len(np.unique(labels))

        logger.info(f"Resolution {res} ({n_clusters} clusters)...")

        if n_clusters > 1:
            sil_score = silhouette_score(X, labels)
            sil_samples = silhouette_samples(X, labels)
            cluster_sil = {
                cid: sil_samples[labels == cid].mean() for cid in np.unique(labels)
            }
            min_cluster_sil = min(cluster_sil.values())
            max_cluster_sil = max(cluster_sil.values())
        else:
            sil_score = np.nan
            min_cluster_sil = np.nan
            max_cluster_sil = np.nan

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

        logger.info(f"Silhouette: {sil_score:.3f}")

    metrics_df = pd.DataFrame(results)

    # Log best resolution
    best_idx = metrics_df["silhouette_score"].idxmax()
    best = metrics_df.loc[best_idx]
    logger.info(
        f" Best resolution: {best['resolution']} "
        f"({best['n_clusters']} clusters, silhouette={best['silhouette_score']:.3f})"
    )

    return metrics_df


def plot_dimensionality_reduction(
    adata: sc.AnnData,
    module_dir: Path,
    norm_approach: str,
    n_neighbors: int,
    leiden_key: str = "leiden_best",
    cmap: Any = None,
    config: DimensionReductionModuleConfig | None = None,
) -> None:
    """Create and save dimensionality reduction plots.

    Args:
        adata: AnnData with computed UMAP and clusters
        module_dir: Directory to save figures
        norm_approach: Normalization approach label
        n_neighbors: Number of neighbors used
        leiden_key: Leiden clustering column name
        cmap: Colormap for continuous variables
        config: Configuration object with visualization settings
    """
    if cmap is None:
        cmap = sns.color_palette("crest", as_cmap=True)

    logger.info("Plotting UMAPs...")

    # Build color list
    color_list = ["total_counts", "n_genes_by_counts"]
    if leiden_key in adata.obs.columns:
        color_list.append(leiden_key)
        logger.info(f"  Including '{leiden_key}' in plots")
    else:
        logger.warning(f"  Cluster column '{leiden_key}' not found, skipping")
    # Main UMAP plot
    sc.pl.umap(
        adata,
        color=color_list,
        ncols=3,
        cmap=cmap,
        wspace=0.4,
        show=False,
        save=f"_{leiden_key}_{norm_approach}_n{n_neighbors}.pdf",
        frameon=False,
    )

    # Observation fields
    if config and config.obs_vis_list:
        logger.info("  Plotting observation fields...")
        sc.pl.umap(
            adata,
            color=config.obs_vis_list,
            ncols=3,
            cmap=cmap,
            wspace=0.4,
            show=False,
            save=f"_{leiden_key}_{norm_approach}_obs_fields.pdf",
            frameon=False,
        )

    # Marker genes
    if config and config.marker_genes:
        logger.info("  Plotting marker genes...")
        sc.pl.umap(
            adata,
            color=config.marker_genes,
            cmap=cmap,
            ncols=4,
            wspace=0.4,
            show=False,
            save=f"_{leiden_key}_markers_{norm_approach}.pdf",
            frameon=False,
        )

    logger.info(f"Plots saved to {module_dir}")


def plot_spatial_distribution(
    adata: sc.AnnData,
    module_dir: Path,
    leiden_key: str = "leiden_best",
) -> None:
    """Plot spatial distribution of clusters for each ROI.

    Args:
        adata: AnnData with cluster annotations
        module_dir: Directory to save plots
        leiden_key: Leiden clustering column name
    """
    if leiden_key not in adata.obs.columns:
        logger.warning(
            f"Cluster column '{leiden_key}' not found, skipping spatial plots"
        )
        return

    if "ROI" not in adata.obs.columns:
        logger.warning("'ROI' column not found, skipping spatial plots")
        return

    logger.info(f"Plotting spatial distribution of '{leiden_key}'...")

    for roi in adata.obs["ROI"].unique():
        subset = adata[adata.obs["ROI"] == roi]
        sq.pl.spatial_scatter(
            subset,
            library_id="spatial",
            shape=None,
            color=[leiden_key],
            wspace=0.4,
            save=module_dir / f"{leiden_key}_{roi}_spatial.pdf",
        )

    logger.info(f"  Spatial plots saved to {module_dir}")


def run_dimension_reduction(
    config: DimensionReductionModuleConfig,
    io_config: IOConfig,
) -> sc.AnnData:
    """Run complete dimension reduction pipeline.

    Args:
        config: Configuration for dimension reduction
        io_config: IO configuration

    Returns:
        AnnData with computed dimensionality reduction and clustering
    """
    # Setup directories
    module_dir = Path(io_config.output_dir) / config.module_name
    module_dir.mkdir(exist_ok=True, parents=True)
    spatial_plots_dir = module_dir / "spatial_plots"
    spatial_plots_dir.mkdir(exist_ok=True, parents=True)

    sc.settings.figdir = module_dir
    configure_scanpy_figures(str(io_config.output_dir))
    cmap = sns.color_palette("Spectral", as_cmap=True)

    # Load data
    adata = load_data(
        io_config=io_config,
        config=config,
        subsample_strategy=config.subsample_strategy,
    )

    # Ensure PCA is computed (only compute if not already present)
    if "X_pca" not in adata.obsm:
        logger.info(f"Computing PCA (n_comps={config.n_pca})...")
        sc.pp.pca(
            adata, n_comps=config.n_pca, svd_solver="arpack", use_highly_variable=True
        )
    else:
        logger.info("PCA already computed, skipping...")

    # Plot PCA variance
    sc.pl.pca_variance_ratio(
        adata,
        log=True,
        n_pcs=config.n_pca,
        show=False,
        save=f"_{config.module_name}.pdf",
    )

    # Compute neighbors and UMAP
    adata = compute_dimensionality_reduction(
        adata=adata,
        n_pca=config.n_pca,
        n_neighbors=config.n_neighbors,
        min_dist=0.1,
    )

    # Compute clustering - use single resolution from config
    adata = calculate_clusters(adata=adata, res_list=[config.resolution])

    # Evaluate the single resolution
    metrics_df = evaluate_resolutions(
        adata, res_list=[config.resolution], n_pcs=config.n_pca
    )
    metrics_df.to_csv(module_dir / "resolution_metrics.csv", index=False)

    # Set the clustering result as 'leiden_best'
    best_res = metrics_df.loc[metrics_df["silhouette_score"].idxmax(), "resolution"]
    logger.info(f"Setting 'leiden_best' to leiden_res_{best_res}")
    adata.obs["leiden_best"] = adata.obs[f"leiden_res_{best_res}"]

    # Save results
    output_path = module_dir / "adata.h5ad"
    adata.write_h5ad(output_path)
    logger.info(f"Results saved to {output_path}")

    # Plot results
    plot_dimensionality_reduction(
        adata=adata,
        module_dir=module_dir,
        norm_approach=config.norm_approach,
        leiden_key="leiden_best",
        n_neighbors=config.n_neighbors,
        cmap=cmap,
        config=config,
    )

    plot_spatial_distribution(
        adata=adata,
        module_dir=spatial_plots_dir,
        leiden_key="leiden_best",
    )

    logger.info(f"Dimension reduction module '{config.module_name}' complete.")
    return adata
