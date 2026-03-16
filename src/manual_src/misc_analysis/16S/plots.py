"""Plots for 16S letter."""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq


# Functions
def plot_gene_expression_per_roi(
    adata, gene, roi_col="ROI", condition_col="condition", agg_func="sum"
):
    """Plot ROI-level expression for a given gene.

    Args:
    adata : AnnData
        AnnData object containing expression matrix.
    gene : str
        Gene name to plot.
    roi_col : str
        Column in adata.obs defining ROI.
    condition_col : str
        Column in adata.obs defining condition.
    agg_func : str
        Aggregation across spots within ROI ('sum', 'mean', etc).

    Returns:
    plot_df : pd.DataFrame
        Aggregated dataframe used for plotting.
    """
    if gene not in adata.var_names:
        raise ValueError(f"{gene} not found in adata.var_names")

    # Extract expression
    expr = adata[:, gene].X
    if hasattr(expr, "toarray"):  # sparse matrix
        expr = expr.toarray().flatten()
    else:
        expr = expr.flatten()

    # Temporary dataframe
    temp_df = pd.DataFrame(
        {
            f"{gene}_expression": expr,
            "ROI": adata.obs[roi_col].values,
            "condition": adata.obs[condition_col].values,
        }
    )

    # Aggregate per ROI
    plot_df = temp_df.groupby(["ROI", "condition"], as_index=False, observed=True).agg(
        {f"{gene}_expression": agg_func}
    )

    plot_df.rename(columns={f"{gene}_expression": f"{gene}_total_count"}, inplace=True)

    print(f"ROI-level {gene} counts:")
    print(plot_df.head())

    # Plot
    plt.figure(figsize=(6, 6))

    condition_color_dict = {"PM08": "#BC3C29", "IPF": "#E18727"}

    sns.boxplot(
        data=plot_df,
        x="condition",
        y=f"{gene}_total_count",
        palette=condition_color_dict,
    )

    sns.stripplot(
        data=plot_df,
        x="condition",
        y=f"{gene}_total_count",
        color="black",
        size=4,
        jitter=True,
        alpha=0.6,
    )

    plt.xlabel("Condition")
    plt.ylabel(f"{gene} Expression")
    plt.title(f"{gene} Expression per ROI by Condition")
    plt.tight_layout()
    plt.savefig(fig_dir / f"{gene}_expression_per_ROI.pdf")

    return plot_df


# Set directories
dir = Path(
    "/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/output/2026-02-22_analysis_run_HVG2000"
)
fig_dir = dir / "16S/figures"
os.makedirs(fig_dir, exist_ok=True)

# Set figures parameters
cmap = sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True)

# Set figure directory for scanpy
sc.settings.figdir = fig_dir

# Load data
adata = sc.read_h5ad(dir / "subset_adata/adata_subset_COPD_IPF_PM08.h5ad")

# Subset on IPF and PM08
adata = adata[adata.obs["condition"].isin(["IPF", "PM08"])].copy()

# Remove unknown cell types
adata = adata[~adata.obs["level_0_annotation"].str.contains("Unknown", case=False)]
print(adata.obs["level_0_annotation"].value_counts())

# Look at 16S expression by cell type
sc.pl.dotplot(
    adata,
    ["16S", "KRT5", "MRC1"],
    groupby="level_0_annotation",
    cmap=cmap,
    save="_level_0.pdf",
)

sc.pl.dotplot(
    adata, ["16S"], groupby="level_0_annotation", cmap=cmap, save="_level_0_16S.pdf"
)

# Look at 16S expression by condition
sc.pl.dotplot(
    adata,
    ["16S", "KRT5", "MRC1"],
    groupby="condition",
    cmap=cmap,
    save="_condition_multi.png",
)

# Look at 16S expression by condition
sc.pl.dotplot(adata, ["16S"], groupby="condition", cmap=cmap, save="_condition_16S.pdf")

sc.pl.violin(adata, ["16S"], groupby="condition", save="_violin_16S.png")

sc.pl.matrixplot(adata, ["16S"], groupby="condition", cmap=cmap, save="_matrix_16S.pdf")


# Plot per cell type
cell_type_list = adata.obs["level_0_annotation"].unique()

for cell in cell_type_list:
    resolution = "level_0_annotation"
    adata_subset = adata[adata.obs[resolution] == cell].copy()

    sc.pl.dotplot(
        adata_subset,
        ["16S"],
        groupby="condition",
        title=f"{cell}",
        cmap=cmap,
        save=f"_{cell}_condition.pdf",
    )

# Plot sum per ROI
plot_df = plot_gene_expression_per_roi(adata, "16S")

# Plot spatial plots
color_list = ["level_0_annotation", "level_1_annotation", "16S"]
ROI_list = adata.obs["ROI"].unique().tolist()
for roi in ROI_list:
    # Subset on ROI
    adata_roi = adata[adata.obs["ROI"] == roi]

    # Make a folder for the ROI
    roi_dir = fig_dir / roi
    os.makedirs(roi_dir, exist_ok=True)

    # Set figure directory for sc
    sc.settings.figdir = roi_dir

    # plot for each parameter
    for color in color_list:
        sq.pl.spatial_scatter(
            adata_roi,
            library_id="spatial",
            shape=None,
            color=color,
            size=5.0,
            frameon=False,
            figsize=(5, 5),
            wspace=0.4,
            vmax=5,
            title=f"{roi}: {color}",
            save=f"_{roi}_{color}_spatial.pdf",
            linewidths=0,  # removes the dot edge/border
        )

# set figure directory back to main
sc.settings.figdir = fig_dir
