"""Generate misc plots."""

import os
import random
from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import torch
from anndata import AnnData


def seed_everything(seed: int):
    """Set random seed on every random module for reproducibility.

    Args:
        seed: The seed value to set for random number generation.
    """
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = True
    elif torch.backends.mps.is_available():
        torch.mps.manual_seed(seed)
    else:
        pass


def create_df_gene(
    adata: AnnData,
    genes: str | Sequence[str],
) -> pd.DataFrame:
    """Create a dataframe for plotting gene expression.

    Args:
        adata: AnnData object.
        genes: A gene name or list of gene names.

    Returns:
        DataFrame containing metadata and gene expression columns.
    """
    if isinstance(genes, str):
        genes = [genes]

    # Check genes exist
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise ValueError(f"Genes not found in adata.var_names: {missing}")

    # Get gene indices
    gene_indices = adata.var_names.get_indexer(genes)

    # Extract expression matrix (cells x genes)
    expr = adata.X[:, gene_indices]

    # Convert sparse → dense only once
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    else:
        expr = np.asarray(expr)

    # Build dataframe
    df = adata.obs.copy()
    expr_df = pd.DataFrame(
        expr,
        index=adata.obs_names,
        columns=genes,
    )

    df = pd.concat([df, expr_df], axis=1)

    return df


# Set random seed for reproducibility
seed_everything(19960915)

# Set variables
res = "leiden_res_0.5"

# Set directories
dir = Path(
    "/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/output/2026-02-19_analysis_run/"
)


fig_dir = dir / "manual_analysis/plots"
fig_dir.mkdir(parents=True, exist_ok=True)

# Configure scanpy to save figures in our custom directory
sc.settings.figdir = fig_dir


# Set colors
cmap = sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True)
color_palette_level_1 = sns.color_palette("hls", 12)
custom_palette = sns.color_palette(
    ["#516a93", "#736da6", "#a06cad", "#cd68a5"], as_cmap=False
)

# Load data
print(f"Loading data from {dir / 'annotate/adata.h5ad'}...")
adata = sc.read_h5ad(dir / "annotate/adata.h5ad")

print("Data loaded successfully.")

# Gene name
gene_names = ["16S", "COL1A1"]

# Create dataframe for plotting
df = create_df_gene(adata, gene_names)
print(f"Successfully created dataframe with genes: {gene_names}")

# Generate violin plots
for gene in gene_names:
    print(f"Visualizing {gene}.")
    plt.figure(figsize=(8, 6))
    sns.violinplot(
        data=df,
        x="condition",
        y=gene,
        hue="timepoint",
        split=False,
        palette=custom_palette,
    )
    plt.title(f"{gene} Expression by Condition and Timepoint")
    plt.tight_layout()
    output_file = fig_dir / f"{gene}_violin_plot.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Saved plot to {output_file}")
    plt.close()  # Close figure to prevent memory leaks

print("All violin plots generated successfully.")
