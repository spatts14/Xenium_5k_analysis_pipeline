"""Plots for 16S letter."""

import os
from pathlib import Path

import scanpy as sc
import seaborn as sns
import squidpy as sq

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

# Remove unknown cell types
adata = adata[adata.obs["level_0_annotation"] != "Unknown"].copy()
print(adata.obs["level_0_annotation"].value_counts())

# Plots
sc.pl.umap(
    adata,
    color=["condition", "timepoint", "level_0_annotation"],
    save="_conditions.png",
)

# Look at 16S expression by cell type
sc.pl.dotplot(
    adata,
    ["16S", "KRT5", "MRC1"],
    groupby="level_0_annotation",
    cmap=cmap,
    save="_level_0.pdf",
)

# Look at 16S expression by cell type
sc.pl.dotplot(
    adata, ["16S"], groupby="level_0_annotation", cmap=cmap, save="_level_0_16S.pdf"
)

# Look at 16S expression by cell type
sc.pl.dotplot(
    adata,
    ["16S", "KRT5", "MRC1"],
    groupby="condition",
    cmap=cmap,
    save="_condition.png",
)

# Look at 16S expression by cell type
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
            size=1.0,
            frameon=False,
            figsize=(10, 10),
            wspace=0.4,
            vmax=5,
            title=f"{roi}: {color}",
            save=f"_{roi}_{color}_spatial.pdf",
            linewidths=0,  # removes the dot edge/border
        )

# set figure directory back to main
sc.settings.figdir = fig_dir
