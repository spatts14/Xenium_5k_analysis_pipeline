"""Generate composition of celltype plots."""

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


def plot_celltype_composition(
    adata,
    celltype_col: str = "level_0_annotation",
    groupby_col: str = "condition",
    figsize: tuple = (10, 6),
    palette: str = "tab20",
    ylabel: str = "Percentage (%)",
    legend_bbox: tuple = (1.02, 1),
):
    """Plot stacked bar chart showing cell type composition per condition.

    Args:
    adata : AnnData
        Annotated data object with cell type and condition info in .obs
    celltype_col : str
        Column name in adata.obs containing cell type annotations
    groupby_col : str
        Column name in adata.obs to group by (e.g., 'condition', 'sample')
    figsize : tuple
        Figure size (width, height)
    palette : str
        Seaborn/matplotlib color palette name
    title : str
        Plot title (default: auto-generated)
    ylabel : str
        Y-axis label
    legend_bbox : tuple
        Legend position (bbox_to_anchor)

    Returns:
    fig, ax : matplotlib figure and axes objects
    """
    # Validate columns exist
    if celltype_col not in adata.obs.columns:
        raise ValueError(f"Column '{celltype_col}' not found in adata.obs")
    if groupby_col not in adata.obs.columns:
        raise ValueError(f"Column '{groupby_col}' not found in adata.obs")

    # Calculate counts and percentages
    counts = (
        adata.obs.groupby([groupby_col, celltype_col], observed=True)
        .size()
        .reset_index(name="count")
    )
    totals = counts.groupby(groupby_col)["count"].transform("sum")
    counts["percentage"] = (counts["count"] / totals) * 100

    # Pivot for stacking
    pivot_df = counts.pivot(
        index=groupby_col, columns=celltype_col, values="percentage"
    ).fillna(0)

    # Get unique cell types and colors
    cell_types = pivot_df.columns.tolist()

    # Check if color palette exists in adata.uns, otherwise create new one
    palette_key = f"{celltype_col}_colors"
    if palette_key in adata.uns:
        print(f"Using existing color palette from adata.uns['{palette_key}']")
        existing_colors = adata.uns[palette_key]
        # Create mapping from cell types to colors
        cell_type_categories = adata.obs[celltype_col].cat.categories.tolist()
        color_dict = dict(zip(cell_type_categories, existing_colors))
        # Filter to only the cell types present in the plot
        color_dict = {ct: color_dict[ct] for ct in cell_types if ct in color_dict}
    else:
        print(f"No existing palette found, creating new one with '{palette}'")
        colors = sns.color_palette(palette, n_colors=len(cell_types))
        color_dict = dict(zip(cell_types, colors))

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot stacked bars using seaborn barplot for each layer
    x_positions = np.arange(len(pivot_df.index))
    bar_width = 0.7
    bottom = np.zeros(len(pivot_df.index))

    for cell_type in cell_types:
        values = pivot_df[cell_type].values
        sns.barplot(
            x=x_positions,
            y=values,
            color=color_dict[cell_type],
            label=cell_type,
            bottom=bottom,
            ax=ax,
            width=bar_width,
            edgecolor="white",
            linewidth=0.5,
        )
        bottom += values

    # Set x-tick labels
    ax.set_xticks(x_positions)
    ax.set_xticklabels(pivot_df.index, rotation=45, ha="right")

    # Customize plot
    ax.set_title(
        f"Cell Type Composition by {groupby_col.replace('_', ' ').title()}",
        fontsize=14,
    )
    ax.set_xlabel(groupby_col.replace("_", " ").title(), fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_ylim(0, 100)

    # Customize legend
    ax.legend(
        title=celltype_col.replace("_", " ").title(),
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )

    plt.tight_layout()
    plt.savefig(fig_dir / f"celltype_composition_{groupby_col}_{celltype_col}.pdf")

    return fig, ax


# Set random seed for reproducibility
seed_everything(19960915)

# Set variables
color = "level_0_annotation"

# Set directories
dir = Path(
    "/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/output/2026-02-22_analysis_run_HVG2000/"
)

fig_dir = dir / "manual_analysis/plots"
fig_dir.mkdir(parents=True, exist_ok=True)

# Configure scanpy to save figures in our custom directory
sc.settings.figdir = fig_dir

# Set colors
cmap = sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True)
color_palette_level_1 = sns.color_palette("hls", 12)

# Load data
print(f"Loading data from {dir / 'annotate/adata.h5ad'}...")
adata = sc.read_h5ad(dir / "annotate/adata.h5ad")

print("Data loaded successfully.")

# Plot
fig, ax = plot_celltype_composition(
    adata, celltype_col="level_0_annotation", groupby_col="condition"
)

fig, ax = plot_celltype_composition(
    adata, celltype_col="level_0_annotation", groupby_col="ROI"
)

fig, ax = plot_celltype_composition(
    adata, celltype_col="level_0_annotation", groupby_col="timepoint"
)
