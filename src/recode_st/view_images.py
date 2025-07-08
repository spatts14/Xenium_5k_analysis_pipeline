"""Image viewing module."""

# Import packages
import warnings
from logging import getLogger

import scanpy as sc
import squidpy as sq

from recode_st.config import IOConfig, ViewImagesModuleConfig
from recode_st.helper_function import seed_everything
from recode_st.logging_config import configure_logging

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def run_view_images(config: ViewImagesModuleConfig, io_config: IOConfig):
    """Run the image viewing module."""
    # Set variables
    module_dir = io_config.output_dir / config.module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "3_annotate" / "adata.h5ad")

    # View plots
    logger.info("Visualize clusters on tissue...")
    sq.pl.spatial_scatter(
        adata,
        library_id="spatial",
        shape=None,
        outline=False,
        color=["leiden", "total_counts"],
        wspace=0.4,
        size=1,
        save=module_dir / "leiden_clusters.png",
        dpi=300,
    )
    logger.info(f"Saved leiden clusters plot to {module_dir / 'leiden_clusters.png'}")

    # View specific gene expression
    logger.info("Plotting genes of interest on tissue...")
    sq.pl.spatial_scatter(
        adata,
        library_id="spatial",
        color=config.gene_list,
        shape=None,
        size=2,
        img=False,
        save=module_dir / "gene_expression.png",
    )
    logger.info(f"Saved gene expression plot to {module_dir / 'gene_expression.png'}")

    # Save anndata object
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Data saved to {module_dir / 'adata.h5ad'}")
    logger.info("Imaging module completed successfully.")


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.4_view_images")

    # Set seed
    seed_everything(21122023)

    try:
        run_view_images(
            ViewImagesModuleConfig(
                module_name="4_view_images",
                gene_list=("EPCAM", "CD3D", "CD68", "PTPRC", "ACTA2"),
            ),
            IOConfig(),
        )
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
