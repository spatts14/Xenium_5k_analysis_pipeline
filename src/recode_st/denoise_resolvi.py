"""Denoising and resolution enhancement using ResolVI."""

import warnings
from logging import getLogger

import scanpy as sc
import scvi

from recode_st.config import DenoiseResolVIModuleConfig, IOConfig
from recode_st.helper_function import configure_scanpy_figures

# Suppress specific warnings to reduce noise in logs
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)


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

    Parameters
    ----------
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
    -------
    model : RESOLVI
        Trained ResolVI model
    adata : AnnData
        AnnData with denoised representations added
    """
    print("=" * 60)
    print("ResolVI DENOISING")
    print("=" * 60)

    # Store raw counts (ResolVI needs raw counts)
    adata.layers["counts"] = adata.X.copy()

    # Setup AnnData for ResolVI
    scvi.external.RESOLVI.setup_anndata(adata, layer="counts")

    # Initialize model
    model = scvi.external.RESOLVI(
        adata,
        n_latent=n_latent,
        n_hidden=n_hidden,
        n_layers=n_layers,
    )

    print("\nModel architecture:")
    print(f"  Latent dimensions: {n_latent}")
    print(f"  Hidden units: {n_hidden}")
    print(f"  Layers: {n_layers}")
    print(
        f"Training on: {'GPU' if use_gpu and scvi.settings.dl_pin_memory_gpu_training else 'CPU'}"
    )

    # Train model
    print(f"\nTraining ResolVI (max {max_epochs} epochs)...")

    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        early_stopping=early_stopping,
        early_stopping_patience=20,
        check_val_every_n_epoch=5,
        train_size=0.9,
        plan_kwargs={"lr": 1e-3},
    )

    print(f"Training completed at epoch {model.history['elbo_train'].shape[0]}")

    # Extract denoised representations
    print("\nExtracting denoised representations...")

    # 1. Latent representation (for clustering, UMAP, etc.)
    adata.obsm["X_resolvi"] = model.get_latent_representation()

    # 2. Denoised normalized expression
    adata.layers["resolvi_denoised"] = model.get_normalized_expression(
        library_size=1e4  # Normalize to 10,000 counts
    )

    print(f"  Added adata.obsm['X_resolvi']: latent representation ({n_latent} dims)")
    print("  Added adata.layers['resolvi_denoised']: denoised expression")

    return model, adata


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

    logger.info("Starting integration of scRNAseq and spatial transcriptomics data...")

    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "quality_control" / "adata.h5ad")
