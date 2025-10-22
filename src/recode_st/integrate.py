"""Integrate scRNAseq and spatial transcriptomics."""

import warnings  # ? what is the best way to suppress warnings from package inputs?
from logging import getLogger

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

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    logger.info("Starting integration of scRNAseq and spatial transcriptomics data...")

    logger.info("Loading scRNAseq data from HLCA ...")
    adata_ref = sc.read_h5ad(io_config.hlca_path)

    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "2_dimension_reduction" / "adata.h5ad")

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
        adata_ref.write_h5ad(io_config.hlca_path)
        logger.info(f"Processed HLCA reference data saved to {io_config.hlca_path}.")
    else:
        logger.info(
            "HLCA reference data already preprocessed. "
            "Contains PCA and UMAP computation."
        )

    logger.info("Starting integration...")
    if method == "ingest":
        logger.info("Integrating data using ingest...")

        logger.info("Enuring Xenium data and reference data have the same genes...")
        var_names = adata_ref.var_names.intersection(adata.var_names)
        adata_ref = adata_ref[:, var_names].copy()
        adata = adata[:, var_names].copy()

        # Run ingest to map Xenium data onto HLCA reference
        sc.tl.ingest(
            adata,
            adata_ref,
            obs=["cell_type"],  # Annotation column to use in adata_ref.obs
            # For core HLCA, "cell_type"
            embedding_method="pca",  # or 'umap'
            labeling_method="knn",  # how to transfer labels
            neighbors_key=None,  # use default neighbors from reference
            inplace=True,
        )
        adata.obs["predicted_cell_type"] = adata.obs["cell_type"]
        del adata.obs["cell_type"]
    elif method == "scANVI":  # Placeholder for scANVI integration
        logger.info("Integrating data using scANVI...")
        # Placeholder for scANVI integration code
        raise NotImplementedError("scANVI integration method not yet implemented.")
    else:
        raise NotImplementedError(f"Integration method {method} not implemented.")

    logger.info("Visualize data following label transfer...")
    sc.pl.umap(
        adata,
        color=["leiden", "predicted_cell_type"],
        title="Xenium data mapped to HLCA",
        save=f"_{config.module_name}_{method}.png",  # save figure
    )

    logger.info("Integration complete.")
