"""Integrate scRNAseq and spatial transcriptomics."""

import warnings  # ? what is the best way to suppress warnings from package inputs?
from logging import getLogger
from pathlib import Path

import pandas as pd

# import torch
import scanpy as sc

from recode_st.config import IntegrateModuleConfig, IOConfig

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def run_integration(config: IntegrateModuleConfig, io_config: IOConfig):
    """Integrate scRNAseq and spatial transcriptomics data using Scanorama.

    Args:
        config (IntegrateModuleConfig): Integration module configuration object.
        io_config (IOConfig): IO configuration object.
    """
    method = config.method

    module_dir = io_config.output_dir / config.module_name
    hlca_path = Path(io_config.hlca_path)
    gene_id_dict_path = Path(io_config.gene_id_dict_path)

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set the directory where to save the ScanPy figures
    sc.settings.figdir = module_dir

    logger.info("Starting integration of scRNAseq and spatial transcriptomics data...")

    logger.info("Loading scRNAseq data from HLCA ...")
    adata_ref = sc.read_h5ad(hlca_path)

    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(hlca_path / "2_dimension_reduction" / "adata.h5ad")

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

    # Add ensembl_id to spatial transcriptomics data
    adata.var["ensembl_id"] = adata.var.index.map(gene_id_dict["ensembl_id"])

    # List of genes shared between datasets based on ensembl IDs
    var_names = adata_ref.var_names.intersection(adata.var["ensembl_id"])
    logger.info(f"Number of common genes: {len(var_names)}")

    # Subset spatial transcriptomics data to common genes
    mask = adata.var["ensembl_id"].isin(
        var_names
    )  # mask to filter genes based on ensembl IDs
    adata_ingest = adata[:, mask].copy()

    # Subset reference datasets to common genes
    adata_ref = adata_ref[:, var_names].copy()

    # Confirm that both datasets have the same genes
    logger.info(f"HLCA: {adata_ref.shape}")
    logger.info(f"ST dataset: {adata_ingest.shape}")

    logger.info("Checking if need to process reference scRNAseq data...")
    # Confirm that PCA and UMAP have been computed for reference data
    if "X_pca" not in adata_ref.obsm or "X_umap" not in adata_ref.obsm:
        logger.info("Preprocessing HLCA reference data...")
        sc.pp.normalize_total(adata_ref, target_sum=1e4)
        sc.pp.log1p(adata_ref)
        sc.pp.highly_variable_genes(adata_ref, n_top_genes=2000)
        sc.tl.pca(adata_ref, n_comps=50)
        sc.pp.neighbors(adata_ref, n_neighbors=15, n_pcs=40)
        sc.tl.umap(adata_ref)
        sc.pl.umap(
            adata_ref,
            color="cell_type",
            title="HLCA reference data UMAP",
            save=f"_{config.module_name}_hlca_umap.png",
        )
        logger.info("Finished preprocessing HLCA reference data.")
        adata_ref.write_h5ad(hlca_path)
        logger.info(f"Processed HLCA reference data saved to {hlca_path}.")
    else:
        logger.info(
            "HLCA reference data already preprocessed. "
            "Contains PCA and UMAP computation."
        )

    logger.info("Starting integration...")
    if method == "ingest":
        logger.info("Integrating data using ingest...")

        # Run ingest to map Xenium data onto HLCA reference
        sc.tl.ingest(
            adata_ingest,
            adata_ref,
            obs=["cell_type"],  # Annotation column to use in adata_ref.obs
            # For core HLCA, "cell_type"
            embedding_method="pca",  # or 'umap'
            labeling_method="knn",  # how to transfer labels
            neighbors_key=None,  # use default neighbors from reference
            inplace=True,
        )

        # Rename predicted cell type column
        adata_ingest.obs["predicted_cell_type"] = adata_ingest.obs["cell_type"]
        del adata_ingest.obs["cell_type"]

        # Copy cell type predictions back to original adata
        adata.obs["predicted_cell_type"] = adata_ingest.obs.loc[
            adata.obs_names, "predicted_cell_type"
        ]

        logger.info("Ingest integration complete.")

    elif method == "scANVI":  # Placeholder for scANVI integration
        logger.info("Integrating data using scANVI...")
        # Placeholder for scANVI integration code
        raise NotImplementedError("scANVI integration method not yet implemented.")
    else:
        raise NotImplementedError(f"Integration method {method} not implemented.")

    logger.info("Visualize data following label transfer...")

    logger.info(f"Columns in adata {adata.obs.columns}...")

    color_list = ["condition", "ROI", "predicted_cell_type"]

    logger.info("Plotting UMAPs...")
    sc.pl.umap(
        adata,
        color=color_list,
        title="Xenium data mapped to HLCA",
        save=f"_{config.module_name}_{method}.png",
    )  # save figure

    sc.pl.umap(
        adata,
        color=["leiden_clusters", "predicted_cell_type"],
        title="Xenium data mapped to HLCA",
        save=f"_{config.module_name}_{method}_leiden.png",
    )  # save figure

    logger.info(f"UMAP plot saved to {sc.settings.figdir}")

    logger.info("Saving integrated data...")
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Integrated data saved to {module_dir}.")

    logger.info("Integration complete.")
