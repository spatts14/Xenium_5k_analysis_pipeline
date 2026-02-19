"""Generate misc plots."""

import os
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
import torch


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


def create_df_gene(adata, gene_name):
    """Create a dataframe for plotting gene expression.

    Args:
        adata: The AnnData object containing the data.
        gene_name: The name of the gene to extract expression values for.

    Returns:
        A dataframe with columns 'gene_expression' and metadata.
    """
    # Get gene index
    gene_idx = adata.var_names.get_loc(gene_name)

    # Extract expression vector
    expr = adata.X[:, gene_idx]

    # Handle sparse matrix
    if hasattr(expr, "toarray"):
        expr = expr.toarray().ravel()
    else:
        expr = np.asarray(expr).ravel()

    # Create dataframe
    df = adata.obs.copy()
    df[gene_name] = expr

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
gene_name = "16S"

# Create dataframe for plotting
df = create_df_gene(adata, gene_name)

print(f"Visualizing {gene_name}.")
sns.violinplot(
    data=df,
    x="condition",
    y=gene_name,
    hue="timepoint",
    split=True,
    palette=custom_palette,
)
plt.savefig(fig_dir / f"{gene_name}_violin_plot.png", dpi=300, bbox_inches="tight")
