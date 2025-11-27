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

# Variables
CELL_TYPE_TOP = "ingest_pred_cell_type"
CELL_TYPE_LEVEL = "transf_ann_level_2_label"
KEY_CELL_TYPE = "myeloid cell"


def drug2cell_calculation(adata, module_dir):
    """Calculate drug2cell scores.

    Args:
        adata (anndata.AnnData): Annotated data object
        module_dir (str): Directory to save module outputs

    Returns:
        adata (anndata.AnnData): adata with drug2cell scores added
    """
    # Computes the mean of the expression of each gene group in each cell.
    # By default, the function will load a set of ChEMBL drugs and their targets
    try:
        d2c.score(adata, use_raw=False)
    except Exception as e:
        logger.error(f"Failed to calculate drug2cell scores: {e}")
        return

    # Validate that drug2cell results were generated
    if "drug2cell" not in adata.uns:
        logger.error(
            "drug2cell scores not found in adata.uns. Drug scoring may have failed."
        )
        return

    logger.info(
        f"Successfully calculated drug scores for "
        f"{adata.uns['drug2cell'].shape[1]} drugs"
    )

    logger.info("Save csv of all drugs identified in dataset")
    drugs_present = adata.uns["drug2cell"].var
    drugs_present.to_csv(module_dir / "drug2cell_drugs.csv")
    return adata


def visualize_drug2cell(config, adata, cell_type_top, cmap):
    """Visualize drug2cell scores.

    Args:
        config (Drug2CellModuleConfig): Configuration object
        adata (anndata.AnnData): Annotated data object with drug2cell scores
        cell_type_top (str): Column name for cell type annotations
        cmap: Color map for visualization

    Returns:
        None
    """
    logger.info("Checking drugs in drug_list against calculated drugs...")
    drug_list = config.drug_list
    available_drugs = adata.uns["drug2cell"].var_names.tolist()
    missing_drugs = [drug for drug in drug_list if drug not in available_drugs]
    if missing_drugs:
        logger.warning(
            f"The following drugs from the config are not found in the "
            f"calculated drug2cell scores: {missing_drugs}"
        )
        drug_list = [drug for drug in drug_list if drug in available_drugs]
        if not drug_list:
            logger.error("No valid drugs found for visualization. Skipping.")
            return

    logger.info("Visualize drug score in UMAP")
    sc.pl.umap(
        adata.uns["drug2cell"],
        color=drug_list,
        color_map=cmap,
        save=f"_{config.module_name}_drug2cell_umap.png",
    )
    logger.info("Visualize drug score in spatial scatter plots")
    for roi in adata.obs["ROI"].unique():
        logger.info(f"Visualizing ROI: {roi}")
        subset = adata[adata.obs["ROI"] == roi]
        sq.pl.spatial_scatter(
            subset.uns["drug2cell"],
            library_id="spatial",
            shape=None,
            color=drug_list,
            wspace=0.4,
            cmap=cmap,
            save=f"_{config.module_name}_{cell_type_top}_scatter.png",
        )

    logger.info("Calculating differential expression...")
    sc.tl.rank_genes_groups(
        adata.uns["drug2cell"], method="wilcoxon", groupby=cell_type_top
    )
    sc.pl.rank_genes_groups_dotplot(
        adata.uns["drug2cell"],
        swap_axes=True,
        dendrogram=False,
        n_genes=3,
        cmap=cmap,
        save=f"_{config.module_name}_{cell_type_top}_dotplot.png",
    )

    logger.info("Visualize only respiratory drugs")
    plot_args = d2c.util.prepare_plot_args(adata.uns["drug2cell"], categories=["R"])
    sc.pl.dotplot(
        adata.uns["drug2cell"],
        groupby=cell_type_top,
        swap_axes=True,
        **plot_args,
        cmap=cmap,
        save=f"_{config.module_name}_{cell_type_top}_respiratory_drugs_by.png",
    )


def celltype_level(
    config,
    adata,
    cell_type: str = KEY_CELL_TYPE,
    cell_type_level: str = CELL_TYPE_LEVEL,
    cell_type_top: str = CELL_TYPE_TOP,
    cmap=None,
):
    """Calculate drug2cell for specified cell type for set cell type level.

    Args:
        config: Configuration object
        adata (anndata.AnnData): Annotated data object
        cell_type (str): Cell type to subset on
        cell_type_level (str): Cell type level or resolution of interest
        cell_type_top (str): Column name for cell type annotations
        cmap: Color map for visualization

    Returns:
        None
    """
    logger.info(f"Subsetting on {cell_type} in {cell_type_level}")
    # Fix critical bug: was using assignment (=) instead of comparison (==)
    if cell_type_level not in adata.obs.columns:
        logger.error(f"Column '{cell_type_level}' not found in adata.obs")
        return

    if cell_type not in adata.obs[cell_type_level].values:
        logger.warning(
            f"Cell type '{cell_type}' not found in column '{cell_type_level}'"
        )
        return

    subset = adata[adata.obs[cell_type_level] == cell_type]

    if len(subset) == 0:
        logger.warning(
            f"No cells found for {cell_type} in {cell_type_level}. "
            "Skipping visualization."
        )
        return

    # Copy drug2cell scores to subset if they exist
    if "drug2cell" in adata.uns:
        subset.uns["drug2cell"] = adata.uns["drug2cell"]
    else:
        logger.error("No drug2cell scores found in adata.uns. Cannot visualize subset.")
        return

    logger.info(f"Visualizing on {cell_type} in {cell_type_level}")
    # Fix parameter name mismatch: CLUSTER_TYPE should be CELL_TYPE_LEVEL_ALL
    visualize_drug2cell(config, adata=subset, cell_type_top=cell_type_top, cmap=cmap)


def run_drug2cell(config: Drug2CellModuleConfig, io_config: IOConfig):
    """Calculate drug score using drug2cell.

    Args:
        config (Drug2CellModuleConfig): Drug2Cell module configuration
        io_config (IOConfig): IO configuration object

    Returns:
        None
    """
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

    logger.info("Loading Xenium data...")
    # Fix path - should load from integrate_ingest module output
    input_file = io_config.output_dir / "3_integrate_ingest" / "adata.h5ad"
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        logger.error("Please run the integrate_ingest module first.")
        return

    adata = sc.read_h5ad(input_file)

    logger.info("Calculating drug score every cell in adata using ChEMBL database")
    drug2cell_calculation(adata, module_dir)

    logger.info("Visualize drug2cell...")
    # Validate required columns exist
    if CELL_TYPE_LEVEL not in adata.obs.columns:
        logger.warning(
            f"Column '{CELL_TYPE_LEVEL}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )
        logger.info("Proceeding with cell_type_level_all only...")
    visualize_drug2cell(config, adata=adata, cell_type_top=CELL_TYPE_TOP, cmap=cmap)

    # Only run cell type specific analysis if the column exists
    if CELL_TYPE_LEVEL in adata.obs.columns:
        logger.info(f"Examine drug score in {KEY_CELL_TYPE} in {CELL_TYPE_LEVEL}")
        celltype_level(
            config, adata, KEY_CELL_TYPE, CELL_TYPE_LEVEL, CELL_TYPE_TOP, cmap
        )
    else:
        logger.info(
            f"Skipping cell type specific analysis - "
            f"column '{CELL_TYPE_LEVEL}' not available"
        )

    logger.info("Drug2Cell module finished.")
