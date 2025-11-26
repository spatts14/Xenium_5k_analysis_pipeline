"""Integrate scRNAseq and STx."""

import warnings
from logging import getLogger
from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# import torch
import scanpy as sc
import seaborn as sns
from scvi.model import SCANVI, SCVI
from sklearn.metrics import confusion_matrix

from recode_st.config import Config, IntegrateModuleConfig, IOConfig
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


def smart_subsample_reference(
    config_base,
    adata_ref: anndata.AnnData,
    target_cells: int = 200000,
    stratify_col: str = REF_CELL_LABEL_COL,
    min_per_type: int = 100,
    sampling_strategy="proportional",
    rare_cell_boost: float = 2.0,
):
    """Subsample reference data while preserving cell type diversity.

    Args:
        config_base (SimpleNamespace or dict): Configuration object used to set seed.
        adata_ref (anndata.AnnData): Reference AnnData object.
        target_cells (int): Target number of cells after subsampling.
        stratify_col (str): Column in adata_ref.obs to stratify sampling by.
        min_per_type (int): Minimum cells to sample per cell type.
        sampling_strategy (str): One of:
            - "proportional": Sample proportionally to original frequencies
            - "balanced": Give equal weight to all cell types
            - "sqrt": Compromise between proportional and balanced (sqrt of proportions)
        rare_cell_boost (float): Multiplier for min_per_type for rare cell types.

    Returns:
        anndata.AnnData: Subsampled reference AnnData object.
    """
    rng = np.random.default_rng(config_base.seed)
    cell_type_counts = adata_ref.obs[stratify_col].value_counts()
    n_types = len(cell_type_counts)
    total_cells = len(adata_ref)

    # Calculate sampling targets based on strategy
    if sampling_strategy == "proportional":
        # Keep original proportions
        sampling_fracs = cell_type_counts / total_cells
        base_targets = (sampling_fracs * target_cells).round().astype(int)

    elif sampling_strategy == "balanced":
        # Equal representation for all cell types
        base_targets = pd.Series(target_cells // n_types, index=cell_type_counts.index)

    elif sampling_strategy == "sqrt":
        # Square root transformation - compromise between proportional and balanced
        sqrt_counts = np.sqrt(cell_type_counts)
        sampling_weights = sqrt_counts / sqrt_counts.sum()
        base_targets = (sampling_weights * target_cells).round().astype(int)

    else:
        raise ValueError(f"Unknown sampling_strategy: {sampling_strategy}")

    # Apply minimum per type, with boost for very rare cell types
    rare_threshold = total_cells * 0.001  # <0.1% of total
    final_targets = {}

    for ct, count in cell_type_counts.items():
        base_target = base_targets[ct]

        # Determine minimum based on rarity
        if count < rare_threshold:
            effective_min = int(min_per_type * rare_cell_boost)
        else:
            effective_min = min_per_type

        # Don't sample more than available, but ensure minimum
        final_targets[ct] = min(count, max(base_target, effective_min))

    # Adjust if we're over target_cells due to minimums
    total_target = sum(final_targets.values())
    if total_target > target_cells:
        # Scale down proportionally, but preserve minimums
        excess = total_target - target_cells
        adjustable_cts = {
            ct: target for ct, target in final_targets.items() if target > min_per_type
        }

        if adjustable_cts:
            # Remove excess proportionally from adjustable cell types
            total_adjustable = sum(adjustable_cts.values())
            for ct in adjustable_cts:
                reduction = int(excess * (adjustable_cts[ct] / total_adjustable))
                final_targets[ct] = max(min_per_type, final_targets[ct] - reduction)

    # Perform sampling
    sampled_indices = []
    for ct, n_sample in final_targets.items():
        idx = adata_ref.obs.loc[adata_ref.obs[stratify_col] == ct].index
        sampled = rng.choice(idx, size=n_sample, replace=False)
        sampled_indices.extend(sampled)

    # Shuffle to avoid any ordering bias
    rng.shuffle(sampled_indices)

    adata_subsampled = adata_ref[sampled_indices].copy()

    # Report sampling statistics
    print(f"Subsampled from {total_cells:,} to {len(adata_subsampled):,} cells")
    print(f"Cell types: {n_types}")
    print(f"Sampling strategy: {sampling_strategy}")
    print("\nCell type distribution:")

    orig_pcts = (cell_type_counts / total_cells * 100).round(2)
    new_counts = adata_subsampled.obs[stratify_col].value_counts()
    new_pcts = (new_counts / len(adata_subsampled) * 100).round(2)

    for ct in cell_type_counts.index[:10]:  # Show top 10
        print(
            f"  {ct}: {cell_type_counts[ct]:,} ({orig_pcts[ct]:.2f}%) → "
            f"{new_counts.get(ct, 0):,} ({new_pcts.get(ct, 0):.2f}%)"
        )

    if n_types > 10:
        print(f"  ... and {n_types - 10} more cell types")

    return adata_subsampled


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


def scVI_integration(config, adata_ref, adata, module_dir):
    """Integrates STx data with a reference scRNA-seq dataset using scVI.

    This function harmonizes a reference single-cell dataset
    (e.g., HLCA) with a STx dataset (e.g., Xenium).

    Args:
        config (SimpleNamespace or dict): Configuration object containing
            module-specific parameters, including the module name used to
            construct figure filenames.
        adata_ref (anndata.AnnData): Reference scRNA-seq AnnData object containing
            precomputed embeddings (PCA or UMAP) and cell-type annotations.
            Should have matching genes with adata_ingest.
        adata (anndata.AnnData): Spatial dataset formatted for integration,
            aligned by gene set with the reference (same genes, same order).
        module_dir (Path): Output directory for saving results

    Returns:
        tuple: (adata_combined, scanvi_model) - Combined dataset and trained model
    """
    logger.info("Integrating data using scVI...")

    # Training parameters
    MAX_EPOCHS_SCVI = 200

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

    # Combine datasets
    logger.info("Combining reference and spatial data...")
    adata_combined = anndata.concat(
        [adata_ref, adata],
        join="inner",  # only keeps genes present in both datasets
        label=BATCH_COL,
        keys=[REFERENCE_KEY, SPATIAL_KEY],
        index_unique="_",
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
        batch_size=512,  # Increase from 128
        early_stopping=True,
        early_stopping_patience=15,
        use_gpu=True,  # Ensure GPU usage
        plan_kwargs={"lr": 1e-3},  # Explicit learning rate
    )

    logger.info("Obtain and visualize latent representation...")
    adata_combined.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()
    sc.pp.pca(adata_combined)
    sc.pp.neighbors(adata_combined, use_rep=SCVI_LATENT_KEY)
    sc.tl.umap(adata_combined)
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


def scANVI_label_transfer(config, adata_combined, scvi_model, module_dir):
    """Integrates STx data with a reference scRNA-seq dataset using scVI.

    This function performs label transfer from a reference single-cell dataset
    (e.g., HLCA) to a STx dataset (e.g., Xenium).
    scANVI uses a semi-supervised approach to transfer cell-type labels.

    Args:
        adata_combined (anndata.AnnData): Reference scRNA-seq AnnData object containing
            precomputed embeddings (PCA or UMAP) and cell-type annotations.
            Should have matching genes with adata_ingest.
        scvi_model (scvi.model.SCVI): Trained scVI model on combined data.
        config (SimpleNamespace or dict): Configuration object containing
        module_dir (Path): Output directory for saving results

    Returns:
        tuple: (adata_combined, scanvi_model) - Combined dataset and trained model
    """
    # Set constant variables
    MAX_EPOCHS_SCANVI = 50

    # Format labels for scANVI
    # Create scANVI labels (reference has labels, spatial is 'Unknown')
    SCANVI_CELLTYPE_KEY = (
        "training_celltype_scanvi"  # new column that scANVI will use training labels.
    )
    UNLABELED_CATEGORY = (
        "Unknown"  # placeholder for cells that do not have known labels aka STx cells
    )

    # Set the label column (SCANVI_CELLTYPE_KEY) to "Unknown" for all cells initially.
    # This ensures that STx cells are marked as unlabeled before we assign ref labels
    adata_combined.obs[SCANVI_CELLTYPE_KEY] = UNLABELED_CATEGORY

    # Check columns in adata_combined
    logger.info(f"Columns in adata_combined: {list(adata_combined.obs.columns)}")

    # Ensure the batch column exists
    if BATCH_COL not in adata_combined.obs.columns:
        logger.error(f"Batch column '{BATCH_COL}' not found")
        return None, None

    # Assign real labels to reference cells
    ref_mask = adata_combined.obs[BATCH_COL] == REFERENCE_KEY

    # Ensure reference label column exists
    if REF_CELL_LABEL_COL not in adata_combined.obs.columns:
        logger.error(f"Cell type column '{REF_CELL_LABEL_COL}' not found")
        return None, None

    # Assign reference labels safely
    adata_combined.obs.loc[ref_mask, SCANVI_CELLTYPE_KEY] = adata_combined.obs.loc[
        ref_mask, REF_CELL_LABEL_COL
    ].astype(str)

    # Convert the column to categorical (recommended for scANVI)
    adata_combined.obs[SCANVI_CELLTYPE_KEY] = adata_combined.obs[
        SCANVI_CELLTYPE_KEY
    ].astype("category")

    logger.info(
        f"Reference cells with labels: {ref_mask.sum()}"
    )  # number of reference cells with labels.
    logger.info(
        f"Percent ref cells labeled: {100 * ref_mask.sum() / adata_combined.n_obs:.2f}%"
    )
    logger.info(
        f"Spatial cells (unlabeled): {(~ref_mask).sum()}"
    )  # number of spatial (unlabeled) cells.
    logger.info(
        f"Percent Spatial cells: {100 * (~ref_mask).sum() / adata_combined.n_obs:.2f}%"
    )
    logger.info(
        f"Unique cell types: {adata_combined.obs[SCANVI_CELLTYPE_KEY].value_counts()}"
    )

    # Initialize scANVI from trained scVI
    logger.info("Initializing scANVI model...")
    SCANVI.setup_anndata(
        adata_combined,
        labels_key=SCANVI_CELLTYPE_KEY,
        unlabeled_category=UNLABELED_CATEGORY,
        batch_key=BATCH_COL,
    )
    scanvi_model = SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category=UNLABELED_CATEGORY,
        labels_key=SCANVI_CELLTYPE_KEY,
    )

    # Train scANVI
    logger.info("Training scANVI model...")
    scanvi_model.train(
        max_epochs=MAX_EPOCHS_SCANVI,
        batch_size=512,  # Increase from 128
        early_stopping=True,
        early_stopping_patience=15,
        use_gpu=True,
    )
    logger.info("scANVI training completed!")

    logger.info("Get latent representation...")
    adata_combined.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation()

    logger.info("Visualizing scANVI latent space...")
    sc.pp.pca(adata_combined)
    sc.pp.neighbors(adata_combined, use_rep=SCANVI_LATENT_KEY)
    sc.tl.umap(adata_combined)
    sc.pl.umap(
        adata_combined,
        color=[BATCH_COL, REF_CELL_LABEL_COL],
        frameon=False,
        ncols=1,
        save=f"_{config.module_name}_scanvi_umap.png",
    )

    logger.info("Saving scANVI model...")
    scanvi_model.save(module_dir / "_scanvi_ref", overwrite=True)

    return adata_combined, scanvi_model


def extract_predictions_and_visualize(adata_combined, scanvi_model, adata, module_dir):
    """Extract predictions from scANVI model and copy them to original adata.

    Args:
    adata_combined : anndata.AnnData
        Combined reference and spatial data
    scanvi_model : scvi.model.SCANVI
        Trained scANVI model
    adata : anndata.AnnData
        Original spatial dataset to which predictions will be copied
    module_dir : Path
        Output directory for saving results

    Returns:
    adata : anndata.AnnData
        Original adata with scANVI predictions added
    """
    logger.info("Extracting scANVI predictions...")

    # Get predictions (hard labels)
    logger.info("Making predictions...")
    adata_combined.obs[SCANVI_LABEL_COL] = scanvi_model.predict(adata_combined)

    logger.info("Visualize predictions...")
    df = (
        adata_combined.obs.groupby([REF_CELL_LABEL_COL, SCANVI_LABEL_COL])
        .size()
        .unstack(fill_value=0)
    )
    norm_df = df / df.sum(axis=0)

    plt.figure(figsize=(8, 8))
    _ = plt.pcolor(norm_df)
    _ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
    _ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xlabel("Predicted")
    plt.ylabel("Observed")
    plt.savefig(
        module_dir / "scanvi_confusion_matrix.png", dpi=300, bbox_inches="tight"
    )

    # Get prediction probabilities
    prediction_probs = scanvi_model.predict(adata_combined, soft=True)
    adata_combined.obsm["scanvi_probs"] = prediction_probs

    # Calculate prediction uncertainty (entropy)
    entropy = -np.sum(prediction_probs * np.log(prediction_probs + 1e-10), axis=1)
    adata_combined.obs["prediction_entropy"] = entropy

    # Calculate maximum probability (confidence)
    max_prob = np.max(prediction_probs, axis=1)
    adata_combined.obs["prediction_max_prob"] = max_prob

    # Extract spatial predictions
    spatial_mask = adata_combined.obs[BATCH_COL] == SPATIAL_KEY
    adata_spatial_predicted = adata_combined[spatial_mask].copy()

    logger.info(f"Spatial cells annotated: {adata_spatial_predicted.n_obs}")
    logger.info("Predicted cell type distribution:")
    logger.info(adata_spatial_predicted.obs[SCANVI_LABEL_COL].value_counts())

    # Map predictions back to original adata using cell names
    # Remove the suffix added by concat (index_unique="_")
    spatial_indices = adata_spatial_predicted.obs_names.str.replace(
        "_Spatial", "", regex=False
    )

    # Ensure indices match
    matching_mask = spatial_indices.isin(adata.obs_names)
    if not matching_mask.all():
        logger.warning(
            f"Some cell names don't match: {matching_mask.sum()}/{len(matching_mask)}"
        )

    # Create mapping
    cell_prediction_map = dict(
        zip(spatial_indices, adata_spatial_predicted.obs[SCANVI_LABEL_COL])
    )
    probability_map = dict(
        zip(spatial_indices, adata_spatial_predicted.obsm["scanvi_probs"])
    )
    entropy_map = dict(
        zip(spatial_indices, adata_spatial_predicted.obs["prediction_entropy"])
    )
    max_prob_map = dict(
        zip(spatial_indices, adata_spatial_predicted.obs["prediction_max_prob"])
    )

    # Copy predictions to original adata
    adata.obs[SCANVI_LABEL_COL] = adata.obs_names.map(cell_prediction_map)
    adata.obs["scANVI_prediction_probs"] = adata.obs_names.map(probability_map)
    adata.obs["scANVI_prediction_entropy"] = adata.obs_names.map(entropy_map)
    adata.obs["scANVI_prediction_confidence"] = adata.obs_names.map(max_prob_map)

    logger.info("scANVI predictions copied to original adata")

    # Create visualizations
    logger.info("Creating scANVI visualizations...")

    # UMAP for overview
    sc.pp.neighbors(adata_combined, use_rep=SCANVI_LATENT_KEY, n_neighbors=15)
    sc.tl.umap(adata_combined)

    # Plot results
    _, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: UMAP by dataset
    sc.pl.umap(adata_combined, color=BATCH_COL, ax=axes[0, 0], show=False)
    axes[0, 0].set_title(BATCH_COL)

    # Plot 2: UMAP by cell type (reference) and predictions
    sc.pl.umap(adata_combined, color=SCANVI_LABEL_COL, ax=axes[0, 1], show=False)
    axes[0, 1].set_title("scANVI Predictions")

    # Plot 3: Prediction uncertainty
    sc.pl.umap(adata_combined, color="prediction_entropy", ax=axes[1, 0], show=False)
    axes[1, 0].set_title("Prediction Uncertainty (Entropy)")

    # Plot 4: Prediction confidence
    sc.pl.umap(adata_combined, color="prediction_max_prob", ax=axes[1, 1], show=False)
    axes[1, 1].set_title("Prediction Confidence (Max Probability)")

    plt.tight_layout()

    # Save plot
    plot_path = module_dir / "scanvi_integration_results.png"
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()

    logger.info(f"scANVI results plot saved: {plot_path}")

    return adata


def prepare_ingest_datasets(gene_id_dict_path, adata_ref, adata):
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

    # Subset reference datasets to common genes - preserve layers!
    logger.info(f"Original reference layers: {list(adata_ref.layers.keys())}")
    adata_ref_subset = adata_ref[:, var_names].copy()

    # Explicitly preserve counts layer if it exists
    if "counts" in adata_ref.layers:
        logger.info("Preserving counts layer during gene subsetting...")
        adata_ref_subset.layers["counts"] = (
            adata_ref[:, var_names].layers["counts"].copy()
        )
        logger.info("Counts layer preserved in subset!")
    else:
        logger.warning("No counts layer found in original reference data!")

    logger.info(f"Subset reference layers: {list(adata_ref_subset.layers.keys())}")

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


def process_subset_reference_data(config, adata_ref_ingest):
    """Preprocesses the reference scRNA-seq dataset if PCA and UMAP are missing.

    Verify that a reference AnnData object has the components required for
    sc.tl.ingest, and compute PCA/UMAP if they are missing.

    Args:
        config (SimpleNamespace or dict): Configuration object containing
            module-specific parameters, including the module name used to
            construct figure filenames.
        adata_ref_ingest (anndata.AnnData): Subset of reference single-cell RNA-seq
            dataset (e.g., HLCA) to be processed or verified.
            The object is modified in place.

    Side Effects:
        - Writes the processed reference AnnData object to `SUBSET_HLCA_INT_SAVE`.
        - Saves a UMAP visualization as a PNG file.

    Logs:
        - Whether preprocessing is needed.
        - Each major preprocessing step.
        - Path to the saved output file.

    Returns:
        adata_ref_ingest (anndata.AnnData): Processed reference AnnData object.
    """
    logger.info("Checking if reference scRNA-seq data needs processing...")

    if not isinstance(adata_ref_ingest, anndata.AnnData):
        raise TypeError("adata_ref_ingest must be an AnnData object.")

    if adata_ref_ingest.n_vars == 0 or adata_ref_ingest.n_obs == 0:
        raise ValueError("adata_ref_ingest is empty or malformed.")

    if adata_ref_ingest.var_names.has_duplicates:
        raise ValueError("adata_ref_ingest.var_names contains duplicates.")

    # Check label column required for visualization / QC
    if REF_CELL_LABEL_COL not in adata_ref_ingest.obs.columns:
        raise ValueError(
            f"Reference data is missing required label column '{REF_CELL_LABEL_COL}'. "
            "This column is needed for downstream mapping and visualization."
        )

    # Validate config has required attributes
    if not hasattr(config, "module_name") or not config.module_name:
        raise ValueError("config must have a non-empty 'module_name' attribute.")

    # Constants for reproducibility
    N_COMPONENTS = 50  # Number of PCA components
    N_NEIGHBORS = 30  # Number of neighbors for graph construction

    # Check and compute PCA if missing
    if "X_pca" in adata_ref_ingest.obsm:
        logger.info("PCA already computed.")
    else:
        logger.info("Calculating PCA...")
        sc.tl.pca(adata_ref_ingest, n_comps=N_COMPONENTS)

    # Check and compute neighbors if missing
    if "neighbors" in adata_ref_ingest.uns:
        logger.info("Neighbors already computed.")
    else:
        logger.info("Calculating neighbors...")
        sc.pp.neighbors(adata_ref_ingest, n_neighbors=N_NEIGHBORS, n_pcs=N_COMPONENTS)

    # Check and compute UMAP if missing
    if "X_umap" in adata_ref_ingest.obsm:
        logger.info("UMAP already computed.")
    else:
        logger.info("Calculating UMAP...")
        sc.tl.umap(adata_ref_ingest)

    # Generate UMAP visualization
    sc.pl.umap(
        adata_ref_ingest,
        color=REF_CELL_LABEL_COL,
        title="HLCA reference subset for ingest UMAP",
        save=f"_{config.module_name}_hlca_subset_ingest_umap.png",
    )
    logger.info(f"UMAP for {REF_CELL_LABEL_COL} saved for adata_ref_ingest")

    # Save processed reference
    adata_ref_ingest.write_h5ad(SUBSET_HLCA_INT_SAVE)
    logger.info("Finished preprocessing")
    logger.info(f"Processed HLCA reference saved to {SUBSET_HLCA_INT_SAVE}")

    return adata_ref_ingest


def ingest_integration(adata_ref, adata, adata_ingest):
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
        adata (anndata.AnnData): Original STx dataset with transferred labels added.
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

    return adata


def visualize_integration(config, cmap, adata):
    """Visualize integration results using UMAPs.

    Args:
        config (SimpleNamespace or dict): Configuration object containing
            module-specific parameters, including the module name used to
            construct figure filenames.
        cmap (_type_): Color map for visualization.
        adata (anndata.AnnData): Original STx dataset to which
            transferred labels will be written.
    """
    logger.info("Generating UMAPs...")
    logger.info(f"Columns in adata: {list(adata.obs.columns)}")

    # Check if UMAP exists, if not compute it
    if "X_umap" not in adata.obsm:
        logger.info("Computing UMAP for visualization...")
        if "X_pca" not in adata.obsm:
            logger.info("Computing PCA first...")
            sc.tl.pca(adata, n_comps=50)
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
        sc.tl.umap(adata)

    for method in [INGEST_LABEL_COL, SCANVI_LABEL_COL]:
        if method not in adata.obs.columns:
            print(f"Missing column: {method}")
            continue

        sc.pl.umap(
            adata,
            color=method,
            title=f"INTEGRATION WITH HLCA - {method}",
            save=f"_{config.module_name}_{method}.png",
            show=False,  # change to True if you want inline display
        )

    # Plot comparison if both methods are available
    color_list = [INGEST_LABEL_COL, SCANVI_LABEL_COL]
    logger.info("Plotting comparison UMAPs...")
    sc.pl.umap(
        adata,
        color=color_list,
        title="COMPARING INTEGRATION APPROACHES",
        ncols=min(2, len(color_list)),
        save=f"_{config.module_name}_comparison.png",
        cmap=cmap,
    )

    logger.info(f"UMAP plots saved to {sc.settings.figdir}")


def compare(adata, module_dir):
    """Compare ingest and scANVI integrations.

    This function computes various metrics to compare the predictions from
    ingest and scANVI methods, including:
    - Agreement/concordance between methods
    - Prediction confidence metrics
    - Cell type distribution comparisons

    Args:
        adata (anndata.AnnData): Original STx dataset with both ingest and
            scANVI predictions.
        module_dir (Path): Output directory for saving comparison results.

    Returns:
        dict: Dictionary containing comparison metrics
    """
    logger.info("Comparing ingest and scANVI integration results...")

    # Check if both prediction columns exist
    if INGEST_LABEL_COL not in adata.obs.columns:
        logger.warning(f"{INGEST_LABEL_COL} not found in adata.obs")
        return None
    if SCANVI_LABEL_COL not in adata.obs.columns:
        logger.warning(f"{SCANVI_LABEL_COL} not found in adata.obs")
        return None

    # Remove cells with missing predictions
    valid_mask = (
        adata.obs[INGEST_LABEL_COL].notna() & adata.obs[SCANVI_LABEL_COL].notna()
    )
    n_valid = valid_mask.sum()
    n_total = adata.n_obs

    logger.info(
        f"Cells with both predictions: {n_valid}/{n_total} "
        f"({100 * n_valid / n_total:.1f}%)"
    )

    if n_valid == 0:
        logger.error("No cells with both predictions found!")
        return None

    adata_valid = adata[valid_mask].copy()

    # 1. Agreement: How many cells have the same label?
    agreement = (
        adata_valid.obs[INGEST_LABEL_COL] == adata_valid.obs[SCANVI_LABEL_COL]
    ).sum()
    agreement_pct = 100 * agreement / n_valid

    logger.info(
        f"Agreement between methods: {agreement}/{n_valid} ({agreement_pct:.2f}%)"
    )

    # 2. Create confusion matrix
    ingest_labels = adata_valid.obs[INGEST_LABEL_COL].values
    scanvi_labels = adata_valid.obs[SCANVI_LABEL_COL].values

    # Get all unique labels
    all_labels = sorted(set(ingest_labels) | set(scanvi_labels))

    # Create confusion matrix
    cm = confusion_matrix(ingest_labels, scanvi_labels, labels=all_labels)
    cm_df = pd.DataFrame(cm, index=all_labels, columns=all_labels)

    # Save confusion matrix
    cm_path = module_dir / "integration_comparison_confusion_matrix.csv"
    cm_df.to_csv(cm_path)
    logger.info(f"Confusion matrix saved to {cm_path}")

    # 3. Prediction confidence metrics (for scANVI)
    confidence_col = "scANVI_prediction_confidence"
    entropy_col = "scANVI_prediction_entropy"

    metrics = {
        "total_cells": n_total,
        "cells_with_both_predictions": n_valid,
        "agreement_count": agreement,
        "agreement_percentage": agreement_pct,
    }

    if confidence_col in adata_valid.obs.columns:
        mean_conf = adata_valid.obs[confidence_col].mean()
        median_conf = adata_valid.obs[confidence_col].median()
        metrics["scANVI_mean_confidence"] = mean_conf
        metrics["scANVI_median_confidence"] = median_conf
        logger.info(f"scANVI mean confidence: {mean_conf:.3f}")
        logger.info(f"scANVI median confidence: {median_conf:.3f}")

    if entropy_col in adata_valid.obs.columns:
        mean_entropy = adata_valid.obs[entropy_col].mean()
        median_entropy = adata_valid.obs[entropy_col].median()
        metrics["scANVI_mean_entropy"] = mean_entropy
        metrics["scANVI_median_entropy"] = median_entropy
        logger.info(f"scANVI mean entropy: {mean_entropy:.3f}")
        logger.info(f"scANVI median entropy: {median_entropy:.3f}")

    # 4. Cell type distribution comparison
    ingest_dist = adata_valid.obs[INGEST_LABEL_COL].value_counts()
    scanvi_dist = adata_valid.obs[SCANVI_LABEL_COL].value_counts()

    comparison_df = pd.DataFrame(
        {
            "ingest_count": ingest_dist,
            "scANVI_count": scanvi_dist,
        }
    ).fillna(0)
    comparison_df["ingest_percentage"] = 100 * comparison_df["ingest_count"] / n_valid
    comparison_df["scANVI_percentage"] = 100 * comparison_df["scANVI_count"] / n_valid
    comparison_df = comparison_df.sort_values("ingest_count", ascending=False)

    # Save distribution comparison
    dist_path = module_dir / "integration_comparison_distributions.csv"
    comparison_df.to_csv(dist_path)
    logger.info(f"Cell type distribution comparison saved to {dist_path}")

    # 5. Create visualization
    _, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: Agreement percentage
    axes[0, 0].bar(
        ["Agreement", "Disagreement"],
        [agreement_pct, 100 - agreement_pct],
        color=["green", "red"],
        alpha=0.7,
    )
    axes[0, 0].set_ylabel("Percentage (%)")
    axes[0, 0].set_title("Agreement between ingest and scANVI")
    axes[0, 0].set_ylim([0, 100])

    # Plot 2: Top 10 cell types - distribution comparison
    top_n = min(10, len(comparison_df))
    top_cell_types = comparison_df.head(top_n).index

    x = np.arange(top_n)
    width = 0.35
    axes[0, 1].bar(
        x - width / 2,
        comparison_df.loc[top_cell_types, "ingest_percentage"],
        width,
        label="ingest",
        alpha=0.7,
    )
    axes[0, 1].bar(
        x + width / 2,
        comparison_df.loc[top_cell_types, "scANVI_percentage"],
        width,
        label="scANVI",
        alpha=0.7,
    )
    axes[0, 1].set_xlabel("Cell Type")
    axes[0, 1].set_ylabel("Percentage (%)")
    axes[0, 1].set_title("Top 10 Cell Types - Distribution Comparison")
    axes[0, 1].set_xticks(x)
    axes[0, 1].set_xticklabels(top_cell_types, rotation=45, ha="right")
    axes[0, 1].legend()
    axes[0, 1].grid(axis="y", alpha=0.3)

    # Plot 3: scANVI confidence distribution (if available)
    if confidence_col in adata_valid.obs.columns:
        axes[1, 0].hist(
            adata_valid.obs[confidence_col],
            bins=50,
            edgecolor="black",
            alpha=0.7,
        )
        axes[1, 0].axvline(
            mean_conf, color="red", linestyle="--", label=f"Mean: {mean_conf:.3f}"
        )
        axes[1, 0].set_xlabel("Prediction Confidence")
        axes[1, 0].set_ylabel("Number of Cells")
        axes[1, 0].set_title("scANVI Prediction Confidence Distribution")
        axes[1, 0].legend()
        axes[1, 0].grid(alpha=0.3)
    else:
        axes[1, 0].text(
            0.5,
            0.5,
            "Confidence data\nnot available",
            ha="center",
            va="center",
            transform=axes[1, 0].transAxes,
        )
        axes[1, 0].set_title("scANVI Confidence Not Available")

    # Plot 4: Confusion matrix heatmap (top cell types only)
    top_n_cm = min(15, len(all_labels))
    top_labels_cm = all_labels[:top_n_cm]
    cm_subset = cm_df.loc[top_labels_cm, top_labels_cm]

    sns.heatmap(
        cm_subset,
        annot=True,
        fmt="d",
        cmap="Blues",
        ax=axes[1, 1],
        cbar_kws={"label": "Number of Cells"},
    )
    axes[1, 1].set_xlabel("scANVI Prediction")
    axes[1, 1].set_ylabel("ingest Prediction")
    axes[1, 1].set_title("Confusion Matrix (Top 15 Cell Types)")
    plt.setp(axes[1, 1].get_xticklabels(), rotation=45, ha="right")
    plt.setp(axes[1, 1].get_yticklabels(), rotation=0)

    plt.tight_layout()

    # Save comparison plot
    plot_path = module_dir / "integration_comparison.png"
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()

    logger.info(f"Comparison plot saved: {plot_path}")

    # Save metrics to file
    metrics_path = module_dir / "integration_comparison_metrics.txt"
    with open(metrics_path, "w") as f:
        f.write("Integration Comparison Metrics\n")
        f.write("=" * 40 + "\n\n")
        for key, value in metrics.items():
            f.write(f"{key}: {value}\n")
    logger.info(f"Metrics saved to {metrics_path}")

    return metrics


def run_integration(
    config: IntegrateModuleConfig, io_config: IOConfig, config_base: Config
):
    """Integrate scRNAseq and STx data using scANVI and ingest.

    adata_ref_ingest and adata_ingest are used for ingest integration.
    adata_ref and data are used for scANVI integration.

    Args:
        config (IntegrateModuleConfig): Integration module configuration object.
        io_config (IOConfig): IO configuration object.
        config_base (Config): Base configuration object.

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
    viz_assets = configure_scanpy_figures(str(io_config.output_dir))
    cmap = viz_assets["cmap"]

    logger.info("Starting integration of scRNAseq and spatial transcriptomics data...")

    logger.info("Loading scRNAseq data from HLCA ...")
    adata_ref = sc.read_h5ad(ref_path)

    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(io_config.output_dir / "2_dimension_reduction" / "adata.h5ad")

    logger.info("Selecting highly variable genes on reference dataset...")
    sc.pp.highly_variable_genes(
        adata_ref,
        n_top_genes=2000,
        layer="counts",
        flavor="seurat_v3",
        subset=True,  # Subset to highly variable genes for integration
    )
    logger.info(f"Number of highly variable genes selected: {adata_ref.n_vars}")

    logger.info("Sub sample cells from reference dataset for faster integration...")
    if adata_ref.n_obs > 300000:  # Only subsample very large references
        logger.info(f"Subsampling reference from {adata_ref.n_obs} to ~200k cells...")
        adata_ref = smart_subsample_reference(
            config_base=config_base,
            adata_ref=adata_ref,
            target_cells=200000,
            stratify_col=REF_CELL_LABEL_COL,
            min_per_type=100,
            sampling_strategy="proportional",
            rare_cell_boost=2.0,
        )

    # 1. INTEGRATION using scVI and scANVI
    logger.info("Verify adata compatibility for scVI/scANVI integration...")
    verify_counts_layer(adata_ref, "reference")
    verify_counts_layer(adata, "spatial")

    logger.info("Performing integration using: scANVI...")
    logger.info(
        "Step 1. Harmonize scRNAseq reference dataset with STx dataset scVI model..."
    )
    adata_combined, trained_scvi_model = scVI_integration(
        config, adata_ref, adata, module_dir
    )

    logger.info(f"scVI model training complete. adata combined: {adata_combined}")
    logger.info(f"scVI model training complete. Model: {trained_scvi_model}")

    logger.info("Step 2. Transfer labels using scANVI model...")
    adata_combined, trained_scanvi_model = scANVI_label_transfer(
        config, adata_combined, trained_scvi_model, module_dir
    )

    # Extract scANVI predictions and copy to original adata
    logger.info("Extracting predicted labels from scANVI...")
    if adata_combined is not None and trained_scanvi_model is not None:
        adata = extract_predictions_and_visualize(
            adata_combined, trained_scanvi_model, adata, module_dir
        )
    else:
        logger.error("scANVI integration failed. Skipping scANVI predictions.")
        logger.warning("Continuing with ingest integration only...")

    # 2. INTEGRATION USING INGEST
    logger.info("Formatting data for ingest integration...")
    adata_ref_ingest, adata_ingest = prepare_ingest_datasets(
        gene_id_dict_path, adata_ref, adata
    )

    logger.info("Processing subset of reference data for ingest...")
    adata_ref_ingest = process_subset_reference_data(config, adata_ref_ingest)

    logger.info("Starting integration methods...")

    # Perform ingest integration
    logger.info("Performing integration using: sc.tl.ingest...")
    ingest_integration(adata_ref_ingest, adata, adata_ingest)

    # Delete adata_ref_ingest and adata_ingest to free memory
    logger.info("Delete adata_ref_ingest and adata_ingest to free memory...")
    del adata_ref_ingest
    del adata_ingest

    logger.info("Visualize data following label transfer...")
    visualize_integration(config, cmap, adata)

    logger.info("Comparing approaches...")
    compare(adata, module_dir)

    logger.info("Saving integrated data...")
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Integrated data saved to {module_dir}.")

    logger.info("Integration complete.")
