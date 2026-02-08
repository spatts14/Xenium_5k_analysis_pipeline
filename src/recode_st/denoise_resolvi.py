"""Denoising and resolution enhancement using ResolVI."""

import warnings
from logging import getLogger

import numpy as np
import scanpy as sc
import scipy.sparse as sp
import scvi

from recode_st.config import DenoiseResolVIModuleConfig, IOConfig
from recode_st.helper_function import configure_scanpy_figures

# Suppress specific warnings to reduce noise in logs
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)


def check_layers(adata):
    """Check that adata has the necessary layers and obsm for ResolVI.

    Args:
        adata: AnnData object to check.

    Raises:
        ValueError: If required layers or obsm keys are missing.

    Returns:
        adata: AnnData object to check.
    """
    if "counts" not in adata.layers:
        logger.info("adata does not have 'counts' layer!")
        raise ValueError("adata must have raw counts in adata.layers['counts']")
    if "X_spatial" not in adata.obsm:
        # Check common alternative keys
        if "spatial" in adata.obsm:
            logger.info("Copying 'spatial' to 'X_spatial' for ResolVI...")
            adata.obsm["X_spatial"] = adata.obsm["spatial"].copy()
        elif "X_centroid" in adata.obsm:
            logger.info("Copying 'X_centroid' to 'X_spatial' for ResolVI...")
            adata.obsm["X_spatial"] = adata.obsm["X_centroid"].copy()
        elif "x_centroid" in adata.obs and "y_centroid" in adata.obs:
            # Xenium often stores coordinates in obs
            logger.info("Creating 'X_spatial' from x_centroid/y_centroid columns...")
            adata.obsm["X_spatial"] = np.column_stack(
                [adata.obs["x_centroid"].values, adata.obs["y_centroid"].values]
            )
        else:
            logger.error(f"Available obsm keys: {list(adata.obsm.keys())}")
            logger.error(f"Available obs columns: {list(adata.obs.columns)}")
            raise ValueError(
                "Could not find spatial coordinates. ResolVI requires coordinates in "
                "adata.obsm['X_spatial']. Please add them manually."
            )
    logger.info(f"X_spatial shape: {adata.obsm['X_spatial'].shape}")

    return adata


def run_resolvi(
    adata,
    n_latent=30,
    n_hidden=128,
    n_layers=2,
    max_epochs=400,
    batch_size=256,
    early_stopping=True,
    use_gpu=True,
):
    """Run ResolVI denoising on QC-filtered Xenium data.

    Args:
        adata : AnnData
            QC-filtered spatial data
        n_latent : int
            Dimensionality of latent space
        n_hidden : int
            Number of hidden units in neural network layers
        n_layers : int
            Number of hidden layers
        max_epochs : int
            Maximum training epochs
        batch_size : int
            Training batch size
        early_stopping : bool
            Whether to use early stopping
        use_gpu : bool
            Whether to use GPU if available

    Returns:
        model : RESOLVI
            Trained ResolVI model
        adata : AnnData
            AnnData with denoised representations added
    """
    logger.info("Staring ResolVI denoise function...")

    logger.info("Set up adata for ResolVI...")
    scvi.external.RESOLVI.setup_anndata(adata, layer="counts", batch_key="ROI")

    logger.info("Initialize model...")
    model = scvi.external.RESOLVI(
        adata,
        dispersion="gene",
        n_latent=n_latent,
        n_hidden=n_hidden,
        n_layers=n_layers,
        downsample_counts=False,
    )

    logger.info("Model architecture:")
    logger.info(f"Latent dimensions: {n_latent}")
    logger.info(f"Hidden units: {n_hidden}")
    logger.info(f"Layers: {n_layers}")

    # Train model
    logger.info(f"Training ResolVI (max {max_epochs} epochs)...")

    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        early_stopping=False,  # Disabled since no validation set is configured
        check_val_every_n_epoch=None,  # No validation monitoring
    )

    logger.info(f"Training completed at epoch {model.history['elbo_train'].shape[0]}")

    logger.info("Extracting denoised representations...")

    # 1. Latent representation (for clustering, UMAP, etc.)
    adata.obsm["X_resolvi"] = model.get_latent_representation()

    # 2. Denoised normalized expression
    adata.layers["resolvi_denoised"] = model.get_normalized_expression(
        library_size=1e4  # Normalize to 10,000 counts
    )

    print(f"Added adata.obsm['X_resolvi']: latent representation ({n_latent} dims)")
    print("Added adata.layers['resolvi_denoised']: denoised expression")

    return model, adata


def qc_normalize_data(adata):
    """Quality control and normalization of data.

    Args:
        adata: AnnData object to process.

    Returns:
        adata: Processed AnnData object.
    """
    logger.info("Filtering low-quality cells and genes...")
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)

    logger.info("Normalizing by cell area...")
    if "cell_area" in adata.obs.columns:
        cell_area_inv = 1 / adata.obs["cell_area"].values  # shape (n_cells,)

        if sp.issparse(adata.X):
            # Sparse-safe multiplication
            scaling = sp.diags(cell_area_inv)
            adata.X = scaling.dot(adata.X)
        else:
            # Dense case
            adata.X = adata.X * cell_area_inv[:, None]

        logger.info("Log transforming...")
        sc.pp.log1p(adata)  # Log transform
        # Identify highly variable genes

        logger.info("Identifying highly variable genes...")
        sc.pp.highly_variable_genes(
            adata, flavor="seurat", n_top_genes=2000, inplace=True
        )

        logger.info("Scaling data...")
        sc.pp.scale(adata, max_value=10)


def post_resolvi_analysis(adata, resolution=0.5, save_path=None):
    """Perform downstream analysis on ResolVI-denoised data.

    Args:
        adata : AnnData
            Data with ResolVI representations
        resolution : float
            Leiden clustering resolution
        save_path : str, optional
            Path to save figures

    Returns:
        adata : AnnData
            Data with clustering and UMAP added
    """
    logger.info("Post ResolVI analysis to visualize results...")

    logger.info("Calculating clusters for non-ResolVI data for comparison...")
    sc.pp.neighbors(adata, use_rep="X_pca", key_added="neighbors_standard")
    sc.tl.umap(adata, neighbors_key="neighbors_standard")
    adata.obsm["X_umap_standard"] = adata.obsm["X_umap"].copy()
    sc.tl.leiden(
        adata,
        resolution=resolution,
        neighbors_key="neighbors_standard",
        key_added="leiden_standard",
    )

    logger.info("Calculating clusters for non-ResolVI data for comparison...")
    sc.pp.neighbors(
        adata, use_rep="X_resolvi", n_neighbors=40, key_added="neighbors_resolvi"
    )
    sc.tl.umap(adata, neighbors_key="neighbors_resolvi")
    adata.obsm["X_umap_resolvi"] = adata.obsm["X_umap"].copy()
    sc.tl.leiden(
        adata,
        resolution=resolution,
        neighbors_key="neighbors_resolvi",
        key_added="leiden_resolvi",
    )

    logger.info("Comparing UMAPs and clusters before and after ResolVI...")
    # Plot standard UMAP
    adata.obsm["X_umap"] = adata.obsm["X_umap_standard"]
    sc.pl.umap(
        adata,
        color="leiden_standard",
        show=False,
        title="Clusters (Standard - No ResolVI)",
        save=f"_standard_clusters_{resolution}.png",
    )

    # Plot ResolVI UMAP
    adata.obsm["X_umap"] = adata.obsm["X_umap_resolvi"]
    sc.pl.umap(
        adata,
        color="leiden_resolvi",
        show=False,
        title="Clusters (ResolVI)",
        save=f"_resolvi_clusters_{resolution}.png",
    )

    n_clusters = adata.obs["leiden_resolvi"].nunique()
    n_clusters_standard = adata.obs["leiden_standard"].nunique()
    logger.info(f"Found {n_clusters} clusters in leiden_resolvi at {resolution} res")
    logger.info(
        f"Found {n_clusters_standard} clusters in leiden_standard at {resolution} res"
    )

    # UMAP colored by total counts
    sc.pl.umap(
        adata,
        color="total_counts",
        show=False,
        title="Total Counts",
        save="resolVI_total_counts.png",
    )

    sc.pl.embedding(
        adata,
        basis="spatial",
        color="leiden_resolvi",
        show=False,
        title="Spatial Clusters",
        save="_resolVI_embedding.png",
    )

    # Spatial plot colored by cluster
    for roi in adata.obs["ROI"].unique():
        sc.pl.spatial(
            adata[adata.obs["ROI"] == roi],
            color="leiden_resolvi",
            title=f"Spatial Clusters: {roi}",
            show=False,
            save=f"_resolVI_spatial_clusters_{roi}.png",
        )

    return adata


# Define global constants for integration
def run_denoise_resolvi(config: DenoiseResolVIModuleConfig, io_config: IOConfig):
    """Denoising and resolution enhancement using ResolVI.

    Args:
        config (DenoiseResolVIModuleConfig): Denoise with ResolVI module configuration.
        io_config (IOConfig): IO configuration object.

    Returns:
        None
    """
    # Variables

    # Name of the column to store label transfer results in adata.obs
    module_dir = io_config.output_dir / config.module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure directory for this module (overrides global setting)
    sc.settings.figdir = module_dir

    # Set figure settings to ensure consistency across all modules
    configure_scanpy_figures(str(io_config.output_dir))

    # Variables
    RESOLUTION = 1

    logger.info("Staring module to denoise STx data...")

    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(
        io_config.output_dir / "quality_control" / "adata_cell_area.h5ad"
    )

    logger.info("Checking adata layers and obsm for ResolVI compatibility...")
    check_layers(adata)

    logger.info("Run ResolVI model to denoise...")
    model, adata = run_resolvi(
        adata, n_latent=config.n_latent, max_epochs=config.max_epochs
    )

    logger.info(f"Saving model to resolvi_model_{config.n_latent}latent...")
    model.save(module_dir / f"resolvi_model_{config.n_latent}_latent", overwrite=True)
    logger.info(
        f"Model saved to {module_dir / f'resolvi_model_{config.n_latent}_latent'}"
    )

    logger.info("Saving intermediate results after ResolVI...")
    adata.write_h5ad(module_dir / "adata_post_resolvi.h5ad")
    logger.info(f"Intermediate results saved to {module_dir}")

    logger.info("Filter, normalize (cell area), log-transform, and scale data...")
    qc_normalize_data(adata)

    logger.info("Visualize post ResolVI...")
    adata = post_resolvi_analysis(adata, resolution=RESOLUTION, save_path=module_dir)

    # Save final results
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Final results saved to {module_dir}")

    logger.info(f"ResolVI module '{config.module_name}' complete.")
