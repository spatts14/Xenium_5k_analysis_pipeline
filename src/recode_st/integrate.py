"""Integrate scRNAseq and STx."""

import warnings  # ? what is the best way to suppress warnings from package inputs?
from logging import getLogger

import pandas as pd

# import torch
import scanpy as sc
import seaborn as sns
from scvi.model import SCANVI, SCVI

from recode_st.config import IntegrateModuleConfig, IOConfig

warnings.filterwarnings("ignore")

logger = getLogger(__name__)

# Define global variables - SHOULD THESE BE HERE OR INSIDE THE RUN FUNCTION?
INGEST_LABEL_COL = "ingest_pred_cell_type"
SCANVI_LABEL_COL = "scANVI_pred_cell_type"
REF_CELL_LABEL_COL = "cell type"


def prepare_integrated_datasets(gene_id_dict_path, adata_ref, adata):
    """Ensures ST and reference dataset contain are in correct format for integration.

    Both datasets should have:
    1. The same list of genes
    2. Labeled by the ensembl ID rather than gene symbol
    3. Genes are in the same order.
    This needs to be done prior to integration!

    Args:
        gene_id_dict_path (str | Path): Path to gene ID dictionary .csv file.
        First column: named "gene_symbol" containing the gene ID for all genes in
        STx dataset.
        Second column:  named "ensembl_ID" containing the ensembl ID of the
        corresponding gene.
        adata_ref (anndata.AnnData)): Reference scRNAseq AnnData object.
        adata (anndata.AnnData): STx AnnData object.

    Raises:
        FileNotFoundError: _description_
        ValueError: _description_

    Returns:
        anndata.AnnData: formatted adata_ref and adata_ingest AnnData objects.

    """
    logger.info("Confirm Xenium data and reference data have the same genes...")
    # Replace ensembl ID with gene symbols from adata_ref for matching
    try:
        gene_id_dict = pd.read_csv(
            gene_id_dict_path, index_col=0
        )  # dictionary with ensembl and gene symbols
    except FileNotFoundError as e:
        logger.error(f"Gene ID dictionary file not found: {gene_id_dict_path}")
        raise FileNotFoundError(
            f"Gene ID dictionary file not found: {gene_id_dict_path}"
        ) from e
    except pd.errors.ParserError as e:
        logger.error(f"Error parsing gene ID dictionary file: {gene_id_dict_path}")
        raise ValueError(
            f"Error parsing gene ID dictionary file: {gene_id_dict_path}"
        ) from e
    except Exception as e:
        logger.error(f"Error reading gene ID dictionary file: {gene_id_dict_path}: {e}")
        raise

    # Add ensembl_id to STx data
    adata.var["ensembl_id"] = adata.var.index.map(gene_id_dict["ensembl_id"])

    # List of genes shared between datasets based on ensembl IDs
    var_names = adata_ref.var_names.intersection(adata.var["ensembl_id"])
    logger.info(f"Number of common genes: {len(var_names)}")

    # Subset STx data to common genes
    mask = adata.var["ensembl_id"].isin(
        var_names
    )  # Create mask to filter genes based on ensembl IDs
    adata_ingest = adata[:, mask].copy()

    adata_ingest.var_names = adata_ingest.var[
        "ensembl_id"
    ]  # Rename var_names to the Ensembl IDs
    logger.info(f"First 5 gene names in ST: {adata_ingest.var_names[:5].tolist()}")

    # Subset reference datasets to common genes
    adata_ref = adata_ref[:, var_names].copy()

    # Check that the genes are in the same order
    logger.info(
        f"Checking genes list:{set(adata_ref.var_names) == set(adata_ingest.var_names)}"
    )
    # After subsetting both datasets to common genes, add these checks:

    # 1. Check if gene names are identical
    logger.info(
        f"Gene names match: {adata_ref.var_names.equals(adata_ingest.var_names)}"
    )

    # 2. Check if they're in the same order
    logger.info(
        f"Genes in same order: {(adata_ref.var_names == adata_ingest.var_names).all()}"
    )

    # 3. Show first few genes from each
    logger.info(f"First 5 genes in reference: {adata_ref.var_names[:5].tolist()}")
    logger.info(f"First 5 genes in ST: {adata_ingest.var_names[:5].tolist()}")

    # 4. Check for any missing genes in either direction
    missing_in_ref = set(adata_ingest.var_names) - set(adata_ref.var_names)
    missing_in_ingest = set(adata_ref.var_names) - set(adata_ingest.var_names)
    logger.info(f"Genes in ST but not in reference: {len(missing_in_ref)}")
    logger.info(f"Genes in reference but not in ST: {len(missing_in_ingest)}")

    # 5. If order doesn't match, reorder adata_ingest to match adata_ref
    if not (adata_ref.var_names == adata_ingest.var_names).all():
        logger.info("Reordering adata_ingest genes to match reference...")
        adata_ingest = adata_ingest[:, adata_ref.var_names].copy()
        logger.info("Genes reordered successfully.")

    # Confirm that both datasets have the same genes
    logger.info(f"HLCA: {adata_ref.shape}")
    logger.info(f"ST dataset: {adata_ingest.shape}")
    return adata_ref, adata_ingest


def process_reference_data(config, io_config, adata_ref):
    """Preprocesses the reference scRNA-seq dataset if PCA and UMAP are missing.

    This function checks whether the reference AnnData object already contains
    computed PCA and UMAP embeddings. If not, it performs standard single-cell
    preprocessing steps, including normalization, log-transformation, highly
    variable gene selection, PCA, neighbor graph construction, and UMAP embedding.
    The resulting AnnData object is saved to disk and optionally visualized.

    Args:
        config (SimpleNamespace or dict): Configuration object containing
            module-specific parameters, including the module name used to
            construct figure filenames.
        io_config (SimpleNamespace or dict): I/O configuration object with paths
            to input and output data. Must include ``hlca_path`` to save the
            processed reference dataset.
        adata_ref (anndata.AnnData): Reference single-cell RNA-seq dataset (e.g.,
            HLCA) to be processed or verified. The object is modified in place.

    Side Effects:
        - Writes the processed reference AnnData object to ``io_config.hlca_path``.
        - Saves a UMAP visualization as a PNG file.

    Logs:
        - Whether preprocessing is needed.
        - Each major preprocessing step.
        - Path to the saved output file.

    Returns:
        None
    """
    logger.info("Checking if need to process reference scRNAseq data...")
    # Confirm that PCA and UMAP have been computed for reference data
    if "X_pca" not in adata_ref.obsm or "X_umap" not in adata_ref.obsm:
        logger.info("Preprocessing HLCA reference data...")
        sc.pp.normalize_total(adata_ref, target_sum=1e4)
        sc.pp.log1p(adata_ref)
        sc.pp.highly_variable_genes(adata_ref, n_top_genes=5000, flavor="seurat_v3")
        sc.pp.scale(adata_ref, max_value=10)
        sc.tl.pca(adata_ref, n_comps=75, svd_solver="arpack")
        sc.pp.neighbors(adata_ref, n_neighbors=15, n_pcs=40)
        sc.pp.neighbors(adata_ref, n_neighbors=30, n_pcs=75)
        sc.tl.umap(adata_ref)
        sc.pl.umap(
            adata_ref,
            color="cell_type",
            title="HLCA reference data UMAP",
            save=f"_{config.module_name}_hlca_umap.png",
        )
        logger.info("Finished preprocessing HLCA reference data.")
        adata_ref.write_h5ad(io_config.hlca_path)
        logger.info(f"Processed HLCA reference data saved to {io_config.hlca_path}.")
    else:
        logger.info(
            "HLCA reference data already preprocessed. "
            "Contains PCA and UMAP computation."
        )


def perform_ingest_integration(adata_ref, adata, adata_ingest):
    """Integrates STx data with a reference scRNA-seq dataset using ingest.

    This function performs label transfer from a reference single-cell dataset
    (e.g., HLCA) to a STx dataset (e.g., Xenium).
    The ingest approach uses Scanpys `tl.ingest` to
    project query data into the references PCA space and transfer cell-type labels.

    Args:
        adata_ref (anndata.AnnData): Reference scRNA-seq AnnData object containing
            precomputed embeddings (PCA or UMAP) and cell-type annotations.
        adata (anndata.AnnData): Original STx dataset to which
            transferred labels will be written.
        adata_ingest (anndata.AnnData): Copy of the spatial dataset formatted for
            ingestion, aligned by gene set with the reference.

    Returns:
        None
    """
    logger.info("Starting integration...")
    logger.info("Integrating data using ingest...")

    # Run ingest to map Xenium data onto HLCA reference
    sc.tl.ingest(
        adata_ingest,
        adata_ref,
        obs=REF_CELL_LABEL_COL,  # Annotation column to use in adata_ref.obs
        # For core HLCA, "cell_type"
        embedding_method="pca",  # or 'umap'
        labeling_method="knn",  # how to transfer labels
        neighbors_key=None,  # use default neighbors from reference
        inplace=True,
    )

    # Rename predicted cell type column
    adata_ingest.obs[INGEST_LABEL_COL] = adata_ingest.obs[REF_CELL_LABEL_COL]
    del adata_ingest.obs[REF_CELL_LABEL_COL]

    # Copy cell type predictions back to original adata
    adata.obs[INGEST_LABEL_COL] = adata_ingest.obs.loc[
        adata.obs_names, INGEST_LABEL_COL
    ]
    logger.info("Ingest integration complete.")


def perform_scANVI_integration(adata_ref, adata):
    """Integrates STx data with a reference scRNA-seq dataset using scANVI.

    This function performs label transfer from a reference single-cell dataset
    (e.g., HLCA) to a STx dataset (e.g., Xenium).
    The ingest approach uses scANVI to transfer cell-type labels.

    Args:
        adata_ref (anndata.AnnData): Reference scRNA-seq AnnData object containing
            precomputed embeddings (PCA or UMAP) and cell-type annotations.
        adata (anndata.AnnData): Original STx dataset to which
            transferred labels will be written.
        adata_ingest (anndata.AnnData): Copy of the spatial dataset formatted for
            ingestion, aligned by gene set with the reference.

    Returns:
        None
    """
    logger.info("Integrating data using scANVI...")

    logger.info("ADD CODE")
    LABEL_COL = "cell types"  # your ground-truth column
    # 0) Make an unlabeled category
    adata.obs["cell_types_scvi"] = (
        adata.obs[LABEL_COL].astype("string").fillna("Unknown").astype("category")
    )
    # ensure 'Unknown' is a category
    if "Unknown" not in adata.obs["cell_types_scvi"].cat.categories:
        adata.obs["cell_types_scvi"] = adata.obs["cell_types_scvi"].cat.add_categories(
            ["Unknown"]
        )
        adata.obs["cell_types_scvi"] = adata.obs["cell_types_scvi"].fillna("Unknown")

    # 1) Setup WITH labels_key (this is the key piece you were missing)
    SCVI.setup_anndata(
        adata,
        labels_key="cell_types_scvi",  # REQUIRED for SCANVI.from_scvi_model
    )

    # 2) Train SCVI (latent model)
    m = SCVI(adata, n_latent=30)
    m.train()  # tune epochs as needed

    # 3) Build SCANVI from the pretrained SCVI
    scanvi = SCANVI.from_scvi_model(
        m,
        unlabeled_category="Unknown",
    )

    scanvi.train()  # tune epochs as needed

    # 4) Predict and fill only the missing labels
    pred = scanvi.predict()

    adata.obs["cell types_pred_scvi"] = pred

    logger.info("scANVI integration complete.")


def visualize_integration(config, cmap, adata):
    """Visualize integration results using UMAPs.

    Args:
        config (SimpleNamespace or dict): Configuration object containing
            module-specific parameters, including the module name used to
            construct figure filenames.
        cmap (_type_): _description_
        adata (anndata.AnnData): Original STx dataset to which
            transferred labels will be written.
    """
    logger("Generating UMAPs...")
    logger.info(f"Columns in adata {adata.obs.columns}...")

    for method in [INGEST_LABEL_COL, SCANVI_LABEL_COL]:
        sc.pl.umap(
            adata,
            color=method,
            title="INTEGRATION WITH HLCA",
            save=f"_{config.module_name}_{method}_{REF_CELL_LABEL_COL}.png",
        )

    # NEED TO FIX FORMATTING SO ALL PLOTS ARE VISIBLE AND DONT OVERLAP!
    color_list = ["condition", "ROI", INGEST_LABEL_COL, SCANVI_LABEL_COL]

    logger.info("Plotting UMAPs...")
    sc.pl.umap(
        adata,
        color=color_list,
        title="COMPARING INTEGRATION APPROACHES",
        ncol=2,
        save=f"_{config.module_name}_{method}_{REF_CELL_LABEL_COL}.png",
        cmap=cmap,
    )

    logger.info(f"UMAP plot saved to {sc.settings.figdir}")


def compare(adata):
    """Compare ingest and scANVI integrations.

    Args:
        adata (anndata.AnnData): Original STx dataset to which
            transferred labels will be written.
    """
    SCANVI_LABEL_COL
    SCANVI_LABEL_COL


def run_integration(config: IntegrateModuleConfig, io_config: IOConfig):
    """Integrate scRNAseq and STx data using scANVI and ingest.

    Args:
        config (IntegrateModuleConfig): Integration module configuration object.
        io_config (IOConfig): IO configuration object.
    """
    # Variables

    # Name of the column to store label transfer results in adata.obs
    module_dir = io_config.output_dir / config.module_name

    # Paths to input data
    hcla_path = io_config.hlca_path
    gene_id_dict_path = io_config.gene_id_dict_path

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure parameters
    # CLEAN UP CODE BY PUTTING THIS AS A HELPER FUNCTION? BEST WAY TO DO THIS?
    sc.set_figure_params(
        dpi=300,  # resolution of saved figures
        dpi_save=300,  # resolution of saved plots
        frameon=False,  # remove borders
        vector_friendly=True,  # produce vector-friendly PDFs/SVGs
        fontsize=16,  # adjust font size
        facecolor="white",  # background color
        figsize=(15, 4),  # default single-panel figure size
    )
    # Set figure directory where to save scanpy figures
    sc.settings.figdir = module_dir

    # Set color palette
    cmap = sns.color_palette("crest", as_cmap=True)

    logger.info("Starting integration of scRNAseq and spatial transcriptomics data...")

    logger.info("Loading scRNAseq data from HLCA ...")
    adata_ref = sc.read_h5ad(hcla_path)

    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "2_dimension_reduction" / "adata.h5ad")

    logger.info("Formatting data for ingest integration...")
    adata_ref, adata_ingest = prepare_integrated_datasets(
        gene_id_dict_path, adata_ref, adata
    )

    logger.info("Processing reference data...")
    process_reference_data(config, io_config, adata_ref)

    logger.info("Starting integration methods...")
    perform_ingest_integration(adata_ref, adata, adata_ingest)

    perform_scANVI_integration(adata_ref, adata)

    logger.info("Visualize data following label transfer...")
    visualize_integration(config, cmap, adata)

    logger.info("Comparing approaches")
    compare(adata)

    logger.info("Saving integrated data...")
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Integrated data saved to {module_dir}.")

    logger.info("Integration complete.")
