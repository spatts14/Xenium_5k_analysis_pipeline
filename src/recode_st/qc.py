"""Quality control module - Memory optimized version."""

import warnings
from logging import getLogger

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import scTransform
import seaborn as sns
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
    min_cell_area = config.min_cell_area
    max_cell_area = config.max_cell_area
    norm_approach = config.norm_approach

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

    # Minimum transcripts per cell
    min_transcripts = adata.obs["total_counts"].min()
    num_cells_with_min = (adata.obs["total_counts"] == min_transcripts).sum()

    logger.info(f"Minimum number of transcripts per cell: {min_transcripts}")
    logger.info(f"Number of cells with that minimum: {num_cells_with_min}")

    # Minimum genes per cell
    min_genes = adata.obs["n_genes_by_counts"].min()
    num_cells_with_min_genes = (adata.obs["n_genes_by_counts"] == min_genes).sum()

    logger.info(f"Minimum number of genes per cell: {min_genes}")
    logger.info(f"Number of cells with that minimum: {num_cells_with_min_genes}")

    # Maximum transcripts per cell
    max_transcripts = adata.obs["total_counts"].max()
    num_cells_with_max = (adata.obs["total_counts"] == max_transcripts).sum()
    logger.info(f"Maximum number of transcripts per cell: {max_transcripts}")
    logger.info(f"Number of cells with that maximum: {num_cells_with_max}")

    # Maximum genes per cell
    max_genes = adata.obs["n_genes_by_counts"].max()
    num_cells_with_max_genes = (adata.obs["n_genes_by_counts"] == max_genes).sum()
    logger.info(f"Maximum number of genes per cell: {max_genes}")
    logger.info(f"Number of cells with that maximum: {num_cells_with_max_genes}")

    # Scatter plot of number of genes vs total counts
    sns.set_theme(style="white")
    plt.figure(figsize=(6, 5))
    plt.scatter(
        adata.obs["n_genes_by_counts"],
        adata.obs["total_counts"],
        # hue="ROI",
        s=10,  # size of points
        alpha=0.5,  # transparency
    )
    plt.xlabel("Number of genes detected per cell")
    plt.ylabel("Total transcripts per cell")
    plt.title("QC: Genes vs Total Counts per Cell")
    plt.grid(False)
    plt.savefig(module_dir / "qc_genes_vs_total_counts_pre_filter.png", dpi=300)
    plt.close()

    # Plot the summary metrics
    plot_metrics(module_dir, adata)

    # $ QC data #

    # Filter cells and genes
    logger.info("Filtering cells and genes...")
    # Number of cells and genes before filtering
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Apply filters
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Number of cells and genes after filtering
    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars

    # How many were removed
    cells_removed = n_cells_before - n_cells_after
    genes_removed = n_genes_before - n_genes_after

    print(f"Cells removed: {cells_removed} ({cells_removed / n_cells_before:.1%})")
    print(f"Genes removed: {genes_removed} ({genes_removed / n_genes_before:.1%})")

    # Filter cells by cell area
    logger.info(f"Filtering cells with area outside {min_cell_area}-{max_cell_area}")
    adata = adata[
        (adata.obs["cell_area"] >= min_cell_area)
        & (adata.obs["cell_area"] <= max_cell_area),
        :,
    ].copy()

    # Number of cells and genes after filtering
    n_cells_after_area = adata.n_obs
    n_genes_after_area = adata.n_vars

    # How many were removed
    cells_removed = n_cells_after - n_cells_after_area
    genes_removed = n_genes_after - n_genes_after_area

    print(f"Cells removed: {cells_removed} ({cells_removed / n_cells_before:.1%})")
    print(f"Genes removed: {genes_removed} ({genes_removed / n_genes_before:.1%})")

    # Scatter plot of number of genes vs total counts after filtering
    sns.set_theme(style="white")
    plt.figure(figsize=(6, 5))
    plt.scatter(
        adata.obs["n_genes_by_counts"],
        adata.obs["total_counts"],
        # hue="ROI",
        s=10,  # size of points
        alpha=0.5,  # transparency
    )
    plt.xlabel("Number of genes detected per cell")
    plt.ylabel("Total transcripts per cell")
    plt.title("QC: Genes vs Total Counts per Cell")
    plt.grid(False)
    plt.savefig(module_dir / "qc_genes_vs_total_counts_post_filter.png", dpi=300)
    plt.close()

    # Normalize data
    logger.info(f"Normalize data using {norm_approach}...")
    adata.layers["counts"] = adata.X.copy()  # make copy of raw data

    if norm_approach == "scanpy_log":
        # scRNAseq approach
        sc.pp.normalize_total(adata, inplace=True)  # normalize data
        sc.pp.log1p(adata)  # Log transform data
    elif norm_approach == "sctransform":  # ? Should I be log transforming after this?
        # scTransform approach
        vst_out = scTransform.vst(
            adata.X, gene_names=adata.var_names, cell_names=adata.obs_names
        )
        adata.X = vst_out["y"]  # Use Pearson residuals as normalized expression
    elif norm_approach == "cell_area":
        # Check if cell area is available
        if "cell_area" in adata.obs.columns:
            cell_area_inv = 1 / adata.obs["cell_area"].values  # shape (n_cells,)

            if sp.issparse(adata.X):
                # Sparse-safe multiplication
                scaling = sp.diags(cell_area_inv)
                adata.X = scaling.dot(adata.X)
            else:
                # Dense case
                adata.X = adata.X * cell_area_inv[:, None]

            # Log transform
            sc.pp.log1p(adata)
        else:
            logger.warning(
                "Cell area not found in adata.obs; skipping normalization by cell area"
            )
    elif norm_approach == "none":
        # No normalization
        logger.info("No normalization applied.")
    else:
        raise ValueError(f"Normalization approach {norm_approach} not recognized")

    # Save data
    logger.info("Saving filtered and normalized data...")
    adata.write_h5ad(module_dir / f"adata_{norm_approach}.h5ad")
    logger.info(f"Data saved to {module_dir / f'adata_{norm_approach}.h5ad'}")
    logger.info("Quality control completed successfully.")


def plot_metrics(module_dir, adata):
    """Generates and saves histograms summarizing key cell metrics.

    This function creates a 1x4 grid of histograms visualizing:
        1. Total transcripts per cell
        2. Unique transcripts per cell
        3. Area of segmented cells
        4. Nucleus-to-cell area ratio

    Args:
        module_dir (Path or str): Directory path where the output plot will be saved.
        adata (anndata.AnnData): Annotated data matrix with cell metrics
        stored in `adata.obs`.
            Must contain the columns:
            - 'total_counts'
            - 'n_genes_by_counts'
            - 'cell_area'
            - 'nucleus_area'

    Returns:
        None
    """
    # Create 4 subplots
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))

    axs[0].set_title("Total transcripts per cell")
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])

    axs[1].set_title("Unique transcripts per cell")
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1])

    axs[2].set_title("Area of segmented cells")
    sns.histplot(adata.obs["cell_area"], kde=False, ax=axs[2])

    axs[3].set_title("Nucleus ratio")
    sns.histplot(
        adata.obs["nucleus_area"] / adata.obs["cell_area"], kde=False, ax=axs[3]
    )

    fig.suptitle("QC meterics pre-normalization", fontsize=16)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Save figure
    output_path = module_dir / "cell_summary_histograms.png"
    plt.savefig(output_path, dpi=300)
    plt.close()
    logger.info(f"Saved plots to {output_path}")


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
