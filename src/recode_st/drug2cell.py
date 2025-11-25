"""Drug2Cell module."""

import warnings
from logging import getLogger

import drug2cell as d2c
import scanpy as sc
import squidpy as sq

from recode_st.config import Drug2CellModuleConfig, IOConfig
from recode_st.helper_function import configure_scanpy_figures

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def visualize_drug2cell(config, adata, CLUSTER_TYPE, cmap):
    """Visualize drug2cell scores.

    Args:
        config (_type_): _description_
        adata (_type_): _description_
        CLUSTER_TYPE (_type_): _description_
        cmap (_type_): _description_
    Return:
        None
    """
    logger.info("Visualize drug score in UMAP")
    drug_list = config.drug_list
    sc.pl.umap(adata.uns["drug2cell"], color=drug_list, color_map=cmap)

    logger.info("Visualize on tissue")
    sq.pl.spatial_scatter(
        adata.uns["drug2cell"],
        library_id="spatial",
        shape=None,
        color=drug_list,
        wspace=0.4,
        size=2,
        figsize=(6, 6),
        cmap=cmap,
        save=f"_{config.module_name}_{CLUSTER_TYPE}_scatter.png",
    )

    logger.info("Calculating differential expression...")
    sc.tl.rank_genes_groups(
        adata.uns["drug2cell"], method="wilcoxon", groupby=CLUSTER_TYPE
    )
    sc.pl.rank_genes_groups_dotplot(
        adata.uns["drug2cell"],
        swap_axes=True,
        dendrogram=False,
        n_genes=3,
        cmap=cmap,
        save=f"_{config.module_name}_{CLUSTER_TYPE}_dotplot.png",
    )

    logger.info("Visualize only respiratory drugs")
    plot_args = d2c.util.prepare_plot_args(adata.uns["drug2cell"], categories=["R"])
    sc.pl.dotplot(
        adata.uns["drug2cell"],
        groupby=CLUSTER_TYPE,
        swap_axes=True,
        **plot_args,
        cmap=cmap,
        save=f"_{config.module_name}_{CLUSTER_TYPE}_respiratory_drugs_by.png",
    )


def celltype_level(config, adata, CELL_TYPE, CELL_TYPE_LEVEL, cmap):
    """Calculate drug2cell for specified cell type for set cell type level.

    Args:
        config: Configuration
        adata (_type_): _description_
        CELL_TYPE (str): Cell type to subset on
        CELL_TYPE_LEVEL (str): Cell type level or resolution of interest
        cmap (_type_): Color map for visualization.

    Return:
        None
    """
    logger.info(f"Subsetting on {CELL_TYPE} in {CELL_TYPE_LEVEL}")
    subset = adata[adata.obs[CELL_TYPE_LEVEL == CELL_TYPE]]

    logger.info(f"Visualizing on {CELL_TYPE} in {CELL_TYPE_LEVEL}")
    visualize_drug2cell(config, adata=subset, CLUSTER_TYPE=CELL_TYPE, cmap=cmap)


def run_drug2cell(config: Drug2CellModuleConfig, io_config: IOConfig):
    """Calculate drug score using drug2cell.

    Args:
        config (Drug2CellModuleConfig): _description_
        io_config (IOConfig): _description_

    Returns:
        None
    """
    # Variables
    CLUSTER_TYPE = "new_clusters"
    CELL_TYPE_LEVEL = "resolution"
    CELL_TYPE = "macrophage"

    # Name of the column to store label transfer results in adata.obs
    module_dir = io_config.output_dir / config.module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure directory for this module (overrides global setting)
    sc.settings.figdir = module_dir

    # Get shared colormap from global visualization settings
    # This ensures consistency across all modules
    viz_assets = configure_scanpy_figures(str(io_config.output_dir))
    cmap = viz_assets["cmap"]

    logger.info("Starting Drug2Cell module...")

    logger.info("Loading scRNAseq data from HLCA ...")

    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "3_integrate" / "adata.h5ad")

    logger.info("Calculating drug score every cell in adata using ChEMBL database")
    # Computes the mean of the expression of each gene group in each cell.
    # By default, the function will load a set of ChEMBL drugs and their targets
    d2c.score(adata, use_raw=False)

    logger.info("Save csv of all drugs identified in dataset")
    drugs_present = adata.uns["drug2cell"].var
    drugs_present.to_csv(module_dir / "drug2cell_drugs.csv")

    logger.info("Visualize drug2cell...")
    visualize_drug2cell(config, adata, CLUSTER_TYPE, cmap=cmap)

    logger.info(f"Examine drug score in {CELL_TYPE} in {CELL_TYPE_LEVEL}")
    celltype_level(config, adata, CELL_TYPE, CELL_TYPE_LEVEL, cmap=cmap)

    logger.info("Drug2Cell module finished.")
