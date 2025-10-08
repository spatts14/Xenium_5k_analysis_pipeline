"""Quality control module - Memory optimized version."""

import warnings
from logging import getLogger

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.sparse import csr_matrix, issparse
from zarr.errors import PathNotFoundError

from recode_st.config import IOConfig, QualityControlModuleConfig
from recode_st.helper_function import seed_everything
from recode_st.logging_config import configure_logging

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def run_qc(config: QualityControlModuleConfig, io_config: IOConfig):
    """Run quality control on Xenium data."""
    # Set variables
    module_dir = io_config.output_dir / config.module_name
    min_cells = config.min_cells
    min_counts = config.min_counts

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    try:
        logger.info("Loading Xenium data...")
        combined_path = io_config.adata_dir / "combined_adata.h5ad"

        # Read the file normally - we need the data in memory for QC calculations
        adata = sc.read_h5ad(combined_path)

    except PathNotFoundError as err:
        logger.error(f"File not found (or not a valid AnnData file): {combined_path}")
        raise err

    logger.info("Done")

    # $ Calculate and plot metrics

    # Calculate quality control metrics
    sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)

    # Calculate percent negative DNA probe and percent negative decoding count
    cprobes = (
        adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
    )
    cwords = (
        adata.obs["control_codeword_counts"].sum()
        / adata.obs["total_counts"].sum()
        * 100
    )
    logger.info(f"Negative DNA probe count % : {cprobes}")
    logger.info(f"Negative decoding count % : {cwords}")

    # Calculate averages
    avg_total_counts = np.mean(adata.obs["total_counts"])
    logger.info(f"Average number of transcripts per cell: {avg_total_counts:.2f}")

    avg_total_unique_counts = np.mean(adata.obs["n_genes_by_counts"])
    logger.info(f"Average unique transcripts per cell: {avg_total_unique_counts:.2f}")

    area_max = np.max(adata.obs["cell_area"])
    area_min = np.min(adata.obs["cell_area"])
    logger.info(f"Max cell area: {area_max:.2f}")
    logger.info(f"Min cell area: {area_min:.2f}")

    # Plot the summary metrics
    plot_metrics(module_dir, adata)

    # $ QC data #

    # Filter cells and genes
    logger.info("Filtering cells and genes...")
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Normalize data
    logger.info("Normalize data...")

    # Store raw counts efficiently in a sparse format
    # Sparse matrices only store non-zero values, saving lots of memory
    # Convert to CSR (Compressed Sparse Row) format if not already sparse

    if issparse(adata.X):
        # Already sparse - just copy it
        adata.layers["counts"] = adata.X.copy()
    else:
        # Convert dense matrix to sparse format before copying
        # This saves memory because count data has many zeros
        adata.layers["counts"] = csr_matrix(adata.X)

    # Normalize: converts counts to proportions (each cell sums to same total)
    sc.pp.normalize_total(adata, inplace=True)

    # Log transform: log(x + 1) transformation for better distribution
    sc.pp.log1p(adata)

    # Save data
    logger.info("Saving filtered and normalized data...")
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logger.info("Quality control completed successfully.")


def plot_metrics(module_dir, adata):
    """Generates and saves histograms summarizing key cell metrics.

    This function creates a 1x4 grid of histograms visualizing:
        1. Total transcripts per cell
        2. Unique transcripts per cell
        3. Area of segmented cells
        4. Nucleus-to-cell area ratio

    The resulting figure is saved as 'cell_summary_histograms.png'
    in the specified module directory.

    Args:
        module_dir (Path or str): Directory path where the output plot will be saved.
        adata (anndata.AnnData): Annotated data matrix with cell metrics stored in
            `adata.obs`. Must contain the columns:
            - 'total_counts'
            - 'n_genes_by_counts'
            - 'cell_area'
            - 'nucleus_area'

    Returns:
        None: The function saves the plot to disk and logs the output location.
    """
    # Create 4 subplots in a row
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))

    # Plot 1: Total transcripts per cell
    axs[0].set_title("Total transcripts per cell")
    sns.histplot(
        adata.obs["total_counts"],
        kde=False,
        ax=axs[0],
    )

    # Plot 2: Number of unique genes detected per cell
    axs[1].set_title("Unique transcripts per cell")
    sns.histplot(
        adata.obs["n_genes_by_counts"],
        kde=False,
        ax=axs[1],
    )

    # Plot 3: Cell area distribution
    axs[2].set_title("Area of segmented cells")
    sns.histplot(
        adata.obs["cell_area"],
        kde=False,
        ax=axs[2],
    )

    # Plot 4: Ratio of nucleus to cell area
    axs[3].set_title("Nucleus ratio")
    sns.histplot(
        adata.obs["nucleus_area"] / adata.obs["cell_area"],
        kde=False,
        ax=axs[3],
    )

    # Adjust spacing between plots and save
    plt.tight_layout()
    plt.savefig(
        module_dir / "cell_summary_histograms.png",
        dpi=300,
    )

    # Close the plot to free up memory
    plt.close()

    logger.info(f"Saved plots to {module_dir / 'cell_summary_histograms.png'}")


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.1_qc")

    # Set seed for reproducibility
    seed_everything(21122023)

    run_qc(
        QualityControlModuleConfig(
            module_name="1_quality_control",
            min_counts=10,
            min_cells=5,
        ),
        IOConfig(),
    )
