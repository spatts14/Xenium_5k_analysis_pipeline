"""Integrate scRNAseq and STx."""

import warnings  # ? what is the best way to suppress warnings from package inputs?
from logging import getLogger

import anndata
import matplotlib.pyplot as plt
import numpy as np
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
        sc.pp.filter_cells(adata_ref, min_genes=200)
        sc.pp.filter_genes(adata_ref, min_cells=10)
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


def perform_scANVI_integration(adata_ref, adata_ingest):
    """Integrates STx data with a reference scRNA-seq dataset using scANVI.

    This function performs label transfer from a reference single-cell dataset
    (e.g., HLCA) to a STx dataset (e.g., Xenium).
    scANVI uses a semi-supervised approach to transfer cell-type labels.

    Args:
        adata_ref (anndata.AnnData): Reference scRNA-seq AnnData object containing
            precomputed embeddings (PCA or UMAP) and cell-type annotations.
            Should have matching genes with adata_ingest.
        adata_ingest (anndata.AnnData): Spatial dataset formatted for integration,
            aligned by gene set with the reference (same genes, same order).

    Returns:
        tuple: (adata_combined, scanvi_model) - Combined dataset and trained model
    """
    logger.info("Integrating data using scANVI...")

    # Training parameters
    max_epochs_scvi = 200
    max_epochs_scanvi = 100

    # Add dataset labels
    BATCH_COL = "dataset_origin"
    adata_ref.obs[BATCH_COL] = "Reference"
    adata_ingest.obs[BATCH_COL] = "Spatial"

    # Combine datasets
    logger.info("Combining reference and spatial data...")
    adata_combined = anndata.concat(
        [adata_ref, adata_ingest],
        join="inner",
        label=BATCH_COL,
        keys=["Reference", "Spatial"],
        index_unique="_",
    )

    # Create scANVI labels (reference has labels, spatial is 'Unknown')
    labels_key = "cell_type_scanvi"
    unlabeled_category = "Unknown"

    adata_combined.obs[labels_key] = unlabeled_category

    # Set reference labels
    ref_mask = adata_combined.obs[BATCH_COL] == "Reference"
    if REF_CELL_LABEL_COL in adata_combined.obs.columns:
        adata_combined.obs.loc[ref_mask, labels_key] = adata_combined.obs.loc[
            ref_mask, REF_CELL_LABEL_COL
        ]
    else:
        logger.error(f"Cell type column '{REF_CELL_LABEL_COL}' not found")
        return None, None

    logger.info(f"Reference cells with labels: {ref_mask.sum()}")
    logger.info(f"Spatial cells (unlabeled): {(~ref_mask).sum()}")
    logger.info(f"Unique cell types: {adata_combined.obs[labels_key].value_counts()}")

    # Setup scANVI
    logger.info("Setting up scANVI model...")
    SCANVI.setup_anndata(
        adata_combined,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category,
        batch_key=BATCH_COL,
    )

    # Train scVI first (recommended)
    logger.info("Training scVI model...")
    scvi_model = SCVI(adata_combined, n_latent=30, n_hidden=128)
    scvi_model.train(max_epochs=max_epochs_scvi, patience=10, batch_size=128)

    # Initialize scANVI from trained scVI
    logger.info("Initializing scANVI model...")
    scanvi_model = SCANVI.from_scvi_model(
        scvi_model, unlabeled_category=unlabeled_category
    )

    # Train scANVI
    logger.info("Training scANVI model...")
    scanvi_model.train(max_epochs=max_epochs_scanvi, patience=5, batch_size=128)

    logger.info("scANVI training completed!")

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

    # Get latent representation
    logger.info("Getting latent representation...")
    adata_combined.obsm["X_scanvi"] = scanvi_model.get_latent_representation()

    # Get predictions (hard labels)
    logger.info("Making predictions...")
    adata_combined.obs["scanvi_predicted"] = scanvi_model.predict(adata_combined)

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
    spatial_mask = adata_combined.obs["dataset_origin"] == "Spatial"
    adata_spatial_predicted = adata_combined[spatial_mask].copy()

    logger.info(f"Spatial cells annotated: {adata_spatial_predicted.n_obs}")
    logger.info("Predicted cell type distribution:")
    logger.info(adata_spatial_predicted.obs["scanvi_predicted"].value_counts())

    # Map predictions back to original adata using cell names
    # Remove the "_Spatial" suffix added by concat
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
    prediction_map = dict(
        zip(spatial_indices, adata_spatial_predicted.obs["scanvi_predicted"])
    )
    entropy_map = dict(
        zip(spatial_indices, adata_spatial_predicted.obs["prediction_entropy"])
    )
    max_prob_map = dict(
        zip(spatial_indices, adata_spatial_predicted.obs["prediction_max_prob"])
    )

    # Copy predictions to original adata
    adata.obs[SCANVI_LABEL_COL] = adata.obs_names.map(prediction_map)
    adata.obs["scANVI_prediction_entropy"] = adata.obs_names.map(entropy_map)
    adata.obs["scANVI_prediction_confidence"] = adata.obs_names.map(max_prob_map)

    logger.info("scANVI predictions copied to original adata")

    # Create visualizations
    logger.info("Creating scANVI visualizations...")

    # UMAP for overview
    sc.pp.neighbors(adata_combined, use_rep="X_scanvi", n_neighbors=15)
    sc.tl.umap(adata_combined)

    # Plot results
    _, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: UMAP by dataset
    sc.pl.umap(adata_combined, color="dataset_origin", ax=axes[0, 0], show=False)
    axes[0, 0].set_title("Dataset Origin")

    # Plot 2: UMAP by cell type (reference) and predictions
    sc.pl.umap(adata_combined, color="scanvi_predicted", ax=axes[0, 1], show=False)
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

    logger.info(f" scANVI results plot saved: {plot_path}")

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

    # Plot individual methods
    for method in [INGEST_LABEL_COL, SCANVI_LABEL_COL]:
        if method in adata.obs.columns:
            sc.pl.umap(
                adata,
                color=method,
                title=f"INTEGRATION WITH HLCA - {method}",
                save=f"_{config.module_name}_{method}.png",
            )

    # Plot comparison if both methods are available
    color_list = [
        col for col in [INGEST_LABEL_COL, SCANVI_LABEL_COL] if col in adata.obs.columns
    ]

    # Add optional columns if they exist
    for optional_col in ["condition", "ROI"]:
        if optional_col in adata.obs.columns:
            color_list.append(optional_col)

    if len(color_list) > 0:
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
    from sklearn.metrics import confusion_matrix

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

    logger.info("Subset on shared HVGs for integration...")
    # SHOULD WE DO THIS?

    logger.info("Starting integration methods...")

    # 1. Perform ingest integration
    perform_ingest_integration(adata_ref, adata, adata_ingest)

    # 2. Perform scANVI integration (use adata_ingest which has matching genes)
    adata_combined, scanvi_model = perform_scANVI_integration(adata_ref, adata_ingest)

    # 3. Extract scANVI predictions and copy to original adata
    if adata_combined is not None and scanvi_model is not None:
        adata = extract_predictions_and_visualize(
            adata_combined, scanvi_model, adata, module_dir
        )
    else:
        logger.error("scANVI integration failed. Skipping scANVI predictions.")

    logger.info("Visualize data following label transfer...")
    visualize_integration(config, cmap, adata)

    logger.info("Comparing approaches...")
    compare(adata, module_dir)

    logger.info("Saving integrated data...")
    adata.write_h5ad(module_dir / "adata.h5ad")
    logger.info(f"Integrated data saved to {module_dir}.")

    logger.info("Integration complete.")
