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

# Set directories
dir = Path(
    "/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/output/2026-02-22_analysis_run_HVG2000/"
)

fig_dir = dir / "manual_analysis/plots"
fig_dir.mkdir(parents=True, exist_ok=True)

# Configure scanpy to save figures in our custom directory
sc.settings.figdir = fig_dir

# Set colors
sns.set_theme(style="ticks", font_scale=1.0)
# cmap = sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True)
cmap = sns.color_palette("coolwarm", as_cmap=True)
color_palette_level_1 = sns.color_palette("hls", 12)
custom_palette = sns.color_palette(
    ["#516a93", "#736da6", "#a06cad", "#cd68a5"], as_cmap=False
)

# Load data
print(f"Loading data from {dir / 'annotate/adata.h5ad'}...")
adata = sc.read_h5ad(dir / "annotate/adata.h5ad")

print("Data loaded successfully.")

# Gene name
gene_names = ["16S"]

# Ensure genes are in adata.var_names and warn about missing ones
genes_found = [g for g in gene_names if g in adata.var_names]
genes_missing = [g for g in gene_names if g not in adata.var_names]
if genes_missing:
    print(f"Warning: genes not found and skipped: {genes_missing}")

# Extract expression + metadata
expr = adata[:, genes_found].X
if hasattr(expr, "toarray"):
    expr = expr.toarray()

df = pd.DataFrame(expr, columns=genes_found)
df["cell_type"] = adata.obs["level_0_annotation"].values
df["condition"] = adata.obs["condition"].values

# Compute mean expression per (cell_type, condition)
mean_expr = df.groupby(["cell_type", "condition"])[genes_found].mean()

# One heatmap per gene, side by side
n_genes = len(genes_found)
conditions = df["condition"].unique()
cell_types = df["cell_type"].unique()

# Shared color scale across all genes for comparability
vmin = mean_expr.min().min()
vmax = mean_expr.max().max()

fig, axes = plt.subplots(
    1, n_genes, figsize=(10, max(4, len(cell_types) * 0.5 + 2)), sharey=True
)
if n_genes == 1:
    axes = [axes]

for ax, gene in zip(axes, genes_found):
    # Pivot: rows = cell_type, columns = condition
    pivot = mean_expr[gene].unstack("condition")

    sns.heatmap(
        pivot,
        ax=ax,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.5,
        linecolor="white",
        annot=True,
        fmt=".2f",
        annot_kws={"size": 9},
        cbar=(gene == genes_found[-1]),  # colorbar only on last panel
        cbar_kws={"label": "Mean expression", "shrink": 0.7},
        square=False,
    )
    ax.set_title(gene, fontsize=12, fontweight="bold", pad=8)
    ax.set_xlabel("Condition", fontsize=10)
    ax.set_ylabel("Cell type" if gene == genes_found[0] else "", fontsize=10)
    ax.tick_params(axis="x", rotation=45)
    ax.tick_params(axis="y", rotation=0)

fig.suptitle("Mean gene expression per cell type and condition", fontsize=13, y=1.02)
plt.tight_layout()

out_png = fig_dir / "gene_expression_heatmap_16S.png"

plt.savefig(out_png, dpi=150, bbox_inches="tight")
plt.show()
print(f"Saved: {out_png}")
