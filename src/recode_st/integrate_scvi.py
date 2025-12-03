"""Integrate scRNAseq and STx."""

import warnings
from logging import getLogger
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

# import torch
import scanpy as sc
from scvi.model import SCVI

from recode_st.config import Config, IntegrateSCVIModuleConfig, IOConfig
from recode_st.helper_function import configure_scanpy_figures

# Suppress specific warnings to reduce noise in logs
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)

# Define global constants for integration
INGEST_LABEL_COL = "ingest_pred_cell_type"
SCANVI_LABEL_COL = "scANVI_pred_cell_type"
REF_CELL_LABEL_COL = "cell_type"  # Column in reference data with cell type labels
BATCH_COL = "dataset_origin"
REFERENCE_KEY = "Ref"
SPATIAL_KEY = "STx"
SCVI_LATENT_KEY = "X_scVI"
SCANVI_LATENT_KEY = "X_scANVI"
HLCA_INT_SAVE = Path(
    "/rds/general/user/sep22/ephemeral/recode_hlca_full_processed.h5ad"
)
SUBSET_HLCA_INT_SAVE = Path(
    "/rds/general/user/sep22/ephemeral/recode_subset_hlca_full_processed.h5ad"
)


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
        anndata.AnnData: formatted adata_ref_subset and adata_ingest AnnData objects.

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
    missing = adata.var["ensembl_id"].isna().sum()
    if missing > 0:
        logger.warning(f"{missing} genes in STx have no matching Ensembl ID.")

    # List of genes shared between datasets based on ensembl IDs
    var_names = adata_ref.var_names.intersection(adata.var["ensembl_id"])
    logger.info(f"Number of common genes: {len(var_names)}")
    if len(var_names) < 500:
        logger.warning(
            "Warning: Less than 500 common genes found between datasets. "
            "Integration may not perform well."
        )

    # Subset STx data to common genes
    mask = adata.var["ensembl_id"].isin(
        var_names
    )  # Create mask to filter genes based on ensembl IDs
    adata_ingest = adata[:, mask].copy()

    adata_ingest.var_names = adata_ingest.var[
        "ensembl_id"
    ]  # Rename var_names to the Ensembl IDs

    # Check for duplicate gene names in ST dataset
    if adata_ingest.var_names.has_duplicates:
        raise ValueError("Duplicate Ensembl IDs detected in ST dataset after mapping.")

    logger.info(f"First 5 gene names in ST: {adata_ingest.var_names[:5].tolist()}")
    # Subset reference datasets to common genes
    adata_ref_subset = adata_ref[:, var_names].copy()

    # After subsetting both datasets to common genes, add these checks:
    # 1. Check the sets of genes are the same
    ref_genes = set(adata_ref_subset.var_names)
    st_genes = set(adata_ingest.var_names)
    same_set = ref_genes == st_genes
    logger.info(f"Gene sets identical: {same_set}")

    if not same_set:
        missing_in_ref = st_genes - ref_genes
        missing_in_ingest = ref_genes - st_genes
        logger.warning(f"Genes in ST but not in reference: {len(missing_in_ref)}")
        logger.warning(f"Genes in reference but not in ST: {len(missing_in_ingest)}")
        raise ValueError("Datasets do not have identical gene sets. Cannot reorder!")

    # 2. Check order
    same_order = adata_ref_subset.var_names.equals(adata_ingest.var_names)
    logger.info(f"Genes in same order: {same_order}")

    # 3. If order doesn't match, reorder adata_ingest to match adata_ref
    if not same_order:
        logger.info("Reordering adata_ingest genes to match reference...")
        adata_ingest = adata_ingest[:, adata_ref_subset.var_names].copy()
        logger.info("Genes reordered successfully.")

    # Confirm shapes
    logger.info(f"Reference dataset shape: {adata_ref_subset.shape}")
    logger.info(f"ST dataset shape: {adata_ingest.shape}")

    return adata_ref_subset, adata_ingest


def verify_counts_layer(adata, data_type="data"):
    """Verify that counts layer exists and contains valid count data.

    Args:
        adata (anndata.AnnData): AnnData object to verify
        data_type (str): Description of data type for logging
            (e.g., "reference", "spatial")

    Raises:
        ValueError: If counts layer is missing or invalid
    """
    logger.info(f"Verifying {data_type} data integrity...")

    if "counts" in adata.layers:
        logger.info(f"✓ {data_type.capitalize()} data has counts layer")
        logger.info(f"✓ {data_type.capitalize()} layers: {list(adata.layers.keys())}")
        logger.info(f"✓ Counts layer shape: {adata.layers['counts'].shape}")
        logger.info(f"✓ {data_type.capitalize()} data shape: {adata.shape}")

        # Verify counts layer matches main data dimensions
        if adata.layers["counts"].shape != adata.X.shape:
            counts_shape = adata.layers["counts"].shape
            x_shape = adata.X.shape
            logger.error(
                f"✗ Counts layer shape mismatch: counts={counts_shape}, X={x_shape}"
            )
            raise ValueError("Counts layer dimensions don't match main data matrix")

        # Validate counts data quality
        logger.info("Validating counts data quality...")
        counts_sample = adata.layers["counts"][:100, :100]
        counts_arr = counts_sample.A if hasattr(counts_sample, "A") else counts_sample
        is_integer = np.allclose(counts_arr, np.round(counts_arr))
        median_total = np.median(counts_arr.sum(axis=1))
        has_positive_values = (counts_arr >= 0).all()

        if is_integer and median_total > 100 and has_positive_values:
            logger.info("✓ Counts layer contains valid integer count data")
            logger.info(f"✓ Median counts per cell: {median_total:.0f}")
        else:
            logger.warning(
                f"⚠ Counts layer validation issues: integer-like={is_integer}, "
                f"median_total={median_total:.1f}, non_negative={has_positive_values}"
            )
            # Don't fail here, but warn user
            if not has_positive_values:
                logger.error("✗ Counts layer contains negative values!")
                raise ValueError("Invalid counts layer: contains negative values")
    else:
        logger.error(f"✗ {data_type.capitalize()} data missing 'counts' layer!")
        logger.error("Available layers: " + str(list(adata.layers.keys())))
        if data_type.lower() == "reference":
            logger.error(
                "Please run hlca_full_filt_process.py first to create "
                "processed data with counts layer"
            )
        else:
            logger.error(
                "Please ensure your data preprocessing preserves the counts layer"
            )
        raise ValueError(
            f"{data_type.capitalize()} data missing required 'counts' layer"
        )


def scVI_integration_check(adata, batch_key=BATCH_COL, cell_type=REF_CELL_LABEL_COL):
    """Check dataset readiness for scVI/scANVI integration.

    Args:
        adata (anndata.AnnData): AnnData object to validate.
        batch_key (str): Column name in adata.obs for batch ID.
        cell_type (str): Column name in adata.obs for cell type labels.

    Returns:
        bool: True if dataset is ready, False otherwise.
    """
    issues = []

    # Expression matrix
    if adata.X is None:
        issues.append("No expression matrix (adata.X) present")

    # Raw count layer
    if "counts" not in adata.layers:
        issues.append("Missing 'counts' layer in adata.layers (required for scVI)")

    # Batch
    if batch_key not in adata.obs.columns:
        issues.append(f"Missing batch column '{batch_key}' in adata.obs")

    if cell_type not in adata.obs.columns:
        issues.append(f"Missing cell type column '{cell_type}' in adata.obs")

    # Outcome
    if issues:
        print("NOT READY FOR scANVI:")
        for issue in issues:
            print(f"  - {issue}")
        return False

    logger.info("Dataset is ready for scVI/scANVI integration.")
    return True


def scVI_integration(config, adata_combined, module_dir):
    """Integrates STx data with a reference scRNA-seq dataset using scVI.

    This function harmonizes a reference single-cell dataset
    (e.g., HLCA) with a STx dataset (e.g., Xenium).

    Args:
        config (SimpleNamespace or dict): Configuration object containing
            module-specific parameters, including the module name used to
            construct figure filenames.
        adata_combined (anndata.AnnData): Combined AnnData object of reference
        and query (STx) data.
        module_dir (Path): Output directory for saving results

    Returns:
        tuple: (adata_combined, scanvi_model) - Combined dataset and trained model
    """
    logger.info("Integrating data using scVI...")

    # Training parameters
    MAX_EPOCHS_SCVI = 10

    # Setup scVI
    logger.info("Setting up scVI model...")
    SCVI.setup_anndata(
        adata_combined,
        layer="counts",  # must be counts layer
        batch_key=BATCH_COL,  # variable you want to perform harmonization over
    )

    logger.info("Initializing scVI model...")
    scvi_model = SCVI(adata_combined, n_layers=1, n_latent=20, n_hidden=128)
    logger.info("Training SCVI model...")
    scvi_model.train(
        max_epochs=MAX_EPOCHS_SCVI,
        batch_size=1024,  # Increase from 128
        plan_kwargs={"lr": 1e-3},  # Explicit learning rate
    )

    logger.info("Obtain and visualize latent representation...")
    adata_combined.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()
    sc.pp.pca(adata_combined, n_comps=50)
    sc.pp.neighbors(adata_combined, use_rep=SCVI_LATENT_KEY, n_neighbors=15)
    sc.tl.umap(adata_combined, min_dist=0.3, spread=1.0)
    sc.pl.umap(
        adata_combined,
        color=[BATCH_COL, REF_CELL_LABEL_COL],
        frameon=False,
        ncols=1,
        save=f"_{config.module_name}_scvi_umap.png",
    )

    logger.info("Saving scANVI model...")
    scvi_model.save(module_dir / "_scvi_ref", overwrite=True)

    return adata_combined, scvi_model


def run_integration(
    config: IntegrateSCVIModuleConfig, io_config: IOConfig, base_config: Config
):
    """Integrate scRNAseq and STx data using scANVI and ingest.

    adata_ref_ingest and adata_ingest are used for ingest integration.
    adata_ref and data are used for scANVI integration.

    Args:
        config (IntegrateSCVIModuleConfig): Integration module configuration object.
        io_config (IOConfig): IO configuration object.
        base_config (Config): Base configuration object.

    Returns:
        None
    """
    # Variables

    # Name of the column to store label transfer results in adata.obs
    module_dir = io_config.output_dir / config.module_name

    # Paths to input data
    ref_path = io_config.ref_path
    gene_id_dict_path = io_config.gene_id_dict_path

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure directory for this module (overrides global setting)
    sc.settings.figdir = module_dir

    # Get shared colormap from global visualization settings
    # This ensures consistency across all modules
    configure_scanpy_figures(str(io_config.output_dir))

    logger.info("Starting integration of scRNAseq and spatial transcriptomics data...")

    logger.info("Loading scRNAseq data from HLCA ...")
    adata_ref = sc.read_h5ad(ref_path)

    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "2_dimension_reduction" / "adata.h5ad")

    # Subsample reference data if it's too large to avoid memory issues
    if adata_ref.n_obs > 300000:  # Only subsample very large references
        logger.info(f"Subsampling reference from {adata_ref.n_obs:,} to ~200k cells...")
        import numpy as np

        rng = np.random.default_rng(base_config.seed)

        # Simple random sampling - can be improved with stratified sampling later
        target_cells = 200000
        sample_size = min(target_cells, adata_ref.n_obs)
        indices = rng.choice(adata_ref.n_obs, size=sample_size, replace=False)
        adata_ref = adata_ref[indices].copy()
        logger.info(f"Reference subsampled to {adata_ref.n_obs:,} cells")

    # Log available layers for debugging
    logger.info(f"Reference layers: {list(adata_ref.layers.keys())}")
    logger.info(f"Spatial layers: {list(adata.layers.keys())}")

    # Log counts layer shapes (already verified to exist)
    logger.info(f"Reference counts shape: {adata_ref.layers['counts'].shape}")
    logger.info(f"Spatial counts shape: {adata.layers['counts'].shape}")

    # Ensure REF_CELL_LABEL_COL column exists in spatial data
    if REF_CELL_LABEL_COL not in adata.obs.columns:
        adata.obs[REF_CELL_LABEL_COL] = "STx_UNKNOWN"  # placeholder for spatial cells

    logger.info("Verifying datasets for scVI integration...")
    scVI_integration_check(adata_ref, batch_key=BATCH_COL, cell_type=REF_CELL_LABEL_COL)
    scVI_integration_check(adata, batch_key=BATCH_COL, cell_type=REF_CELL_LABEL_COL)

    # 1. INTEGRATION USING INGEST
    logger.info("Subsetting to only shared genes for scVI integration...")
    adata_ref_subset, adata_ingest = prepare_integrated_datasets(
        gene_id_dict_path, adata_ref, adata
    )

    logger.info("Combine reference and query data ...")

    # Memory optimization: Force garbage collection before large operations
    import gc

    gc.collect()

    # Create combined dataset with memory efficiency
    logger.info(
        f"Creating combined dataset from {adata_ref_subset.n_obs:,} ref + "
        f"{adata_ingest.n_obs:,} spatial cells"
    )
    adata_combined = anndata.concat(
        [adata_ref_subset, adata_ingest],
        join="inner",  # only keeps genes present in both datasets
        label=BATCH_COL,
        keys=[REFERENCE_KEY, SPATIAL_KEY],
        index_unique="_",
    )

    # Clean up intermediate objects to free memory
    del adata_ref_subset
    del adata_ingest
    gc.collect()

    logger.info(
        f"Combined dataset created: {adata_combined.n_obs:,} cells x "
        f"{adata_combined.n_vars:,} genes"
    )

    # Verify counts layer was preserved during concat
    if "counts" not in adata_combined.layers:
        logger.error("Counts layer lost during concat! Manually preserving...")
        # Manually recreate counts layer from original data
        ref_mask = adata_combined.obs[BATCH_COL] == REFERENCE_KEY
        stx_mask = adata_combined.obs[BATCH_COL] == SPATIAL_KEY

        # Initialize counts layer
        adata_combined.layers["counts"] = adata_combined.X.copy()  # temporary

        # Fill with actual counts
        adata_combined.layers["counts"][ref_mask, :] = adata_ref.layers["counts"]
        adata_combined.layers["counts"][stx_mask, :] = adata.layers["counts"]
        logger.info("Counts layer manually restored!")
    else:
        logger.info("✓ Counts layer preserved during concat")

    logger.info(f"Combined data layers: {list(adata_combined.layers.keys())}")

    # Check dimensions
    logger.info(
        f"Combined dataset shape: {adata_combined.shape} "
        f"(cells x genes: {adata_combined.n_obs} x {adata_combined.n_vars})"
    )

    # Save combined data before integration
    combined_save_path = module_dir / "adata_combined_before_scVI.h5ad"
    adata_combined.write_h5ad(combined_save_path)
    logger.info(f"Combined data saved to {combined_save_path}")

    # INTEGRATION USING SCVI
    adata_combined, scvi_model = scVI_integration(config, adata_combined, module_dir)

    logger.info("Integration complete SO FAR.")
