"""Level 1 cell type annotation for subset of cells in the Xenium data."""

import os
import random
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq
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


def S_score_G2M_score(adata, subset):
    """Calculate S and G2M scores for cell cycle analysis and plot on UMAP.

    Args:
        adata: An AnnData object containing the single-cell data.
        subset: A string indicating the subset of data being analyzed (e.g., "immune").

    Returns:
        adata: The input AnnData object with added S and G2M scores in adata.obs.
    """
    s_genes = [
        "MCM5",
        "PCNA",
        "TYMS",
        "FEN1",
        "MCM2",
        "MCM4",
        "RRM1",
        "UNG",
        "GINS2",
        "MCM6",
        "CDCA7",
        "DTL",
        "PRIM1",
        "UHRF1",
        "HELLS",
        "RFC2",
        "RPA2",
        "NASP",
        "RAD51AP1",
        "GMNN",
        "WDR76",
        "SLBP",
        "CCNE2",
        "UBR7",
        "POLD3",
        "MSH2",
        "ATAD2",
        "RAD51",
        "RRM2",
        "CDC45",
        "CDC6",
        "EXO1",
        "TIPIN",
        "DSCC1",
        "BLM",
        "CASP8AP2",
        "USP1",
        "CLSPN",
        "POLA1",
        "CHAF1B",
        "BRIP1",
        "E2F8",
    ]

    g2m_genes = [
        "HMGB2",
        "CDK1",
        "NUSAP1",
        "UBE2C",
        "BIRC5",
        "TPX2",
        "TOP2A",
        "NDC80",
        "CKS2",
        "NUF2",
        "CKS1B",
        "MKI67",
        "TMPO",
        "CENPF",
        "TACC3",
        "FAM64A",
        "SMC4",
        "CCNB2",
        "CKAP2L",
        "CKAP2",
        "AURKB",
        "BUB1",
        "KIF11",
        "ANP32E",
        "TUBB4B",
        "GTSE1",
        "KIF20B",
        "HJURP",
        "CDC20",
        "TTK",
        "CDC25C",
        "KIF2C",
        "RANGAP1",
        "NCAPD2",
        "DLGAP5",
        "CDCA3",
        "HMMR",
        "AURKA",
        "PSRC1",
        "ANLN",
        "LBR",
        "CKAP5",
        "CENPE",
        "CTCF",
        "NEK2",
        "G2E3",
        "GAS2L3",
        "CBX5",
        "CENPA",
    ]

    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=s_genes,
        g2m_genes=g2m_genes,
    )

    sc.pl.umap(
        adata,
        color=["S_score", "G2M_score"],
        frameon=False,
        cmap=cmap_blue,
        save=f"_{subset}_s_score_G2M.pdf",
    )


# Set random seed for reproducibility
seed_everything(19960915)

# Set variables
h5ad_file = "adata_subset_Alveolar_epithelium.h5ad"
subset = "alveolar_epithelium"
mannual_annotation = "mannual_annotation"
res = 0.3
annotation_level_0 = subset + "_level_0"
annotation_level_1 = subset + "_level_1"
subset_key = subset + "_" + str(res)

# Set directories
dir = Path(
    "/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/airscape_analysis/"
)

# Save figures
fig_dir = dir / "manual_analysis/figs" / subset
fig_dir.mkdir(parents=True, exist_ok=True)  # Create directory if it doesn't exist
sc.settings.figdir = fig_dir  # Assign it to Scanpy figure directory

# Save spatial figures
fig_dir_spatial = fig_dir / "spatial"
fig_dir_spatial.mkdir(parents=True, exist_ok=True)

# Save files
file_dir = dir / "manual_analysis/files" / subset
file_dir.mkdir(parents=True, exist_ok=True)


# Set colors
cmap = sns.color_palette("Spectral", as_cmap=True)
cmap_blue = sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True)
color_palette_level_1 = sns.color_palette("hls", 12)

# Load data
adata = sc.read_h5ad(dir / f"celltype_subset/{h5ad_file}")

# UMAP plot colored by manual annotation
sc.pl.umap(
    adata,
    color=mannual_annotation,
    legend_loc="right margin",
    legend_fontsize=14,
    frameon=False,
    ncols=2,
    wspace=0.4,
    save=f"_{annotation_level_0}.pdf",
)


# Calculate S and G2M scores and plot on UMAP
S_score_G2M_score(adata, subset)


# Cluster with Leiden algorithm and plot UMAP colored by clusters
sc.tl.leiden(
    adata, resolution=res, key_added=subset_key, flavor="igraph", n_iterations=2
)

# Set colors for Leiden clusters - match number of colors to number of clusters
n_clusters = len(adata.obs[subset_key].cat.categories)
color_palette_leiden = sns.color_palette("hls", n_clusters)
adata.uns[f"{subset_key}_colors"] = [color for color in color_palette_leiden.as_hex()]

# Plot UMAP colored by Leiden clusters
sc.pl.umap(
    adata,
    color=[subset_key],
    frameon=False,
    save=f"_{subset_key}.pdf",
)


# Calculate marker genes for Leiden clusters and plot dotplot
rank_subset_key = "rank_genes_leiden_" + subset_key
sc.tl.rank_genes_groups(
    adata,
    groupby=subset_key,
    method="wilcoxon",
    key_added=rank_subset_key,
)
sc.tl.dendrogram(adata, groupby=subset_key)
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby=subset_key,
    standard_scale="var",
    n_genes=4,
    key=rank_subset_key,
    cmap=cmap_blue,
    save=f"_{subset_key}.pdf",
)

# Get the full ranked genes
markers = sc.get.rank_genes_groups_df(adata, key=rank_subset_key, group=None)

# For top N genes per group (like n_genes in the dotplot)
n_genes = 10
top_dotplot_genes = (
    markers.groupby("group")
    .apply(lambda x: x.sort_values("scores", ascending=False).head(n_genes))
    .reset_index(drop=True)
)

# Save the full ranked genes and top dotplot genes for each cluster
for cluster in adata.obs[subset_key].unique():
    # LFC-sorted markers for the cluster
    markers_subset = markers[markers["group"] == cluster].sort_values(
        "logfoldchanges", ascending=False
    )
    markers_subset.to_csv(
        file_dir / f"{subset}_leiden_{res}_markers_cluster_{cluster}.csv", index=False
    )

    # Top scored genes for the cluster (like in the dotplot)
    top_dotplot_genes_subset = top_dotplot_genes[top_dotplot_genes["group"] == cluster]
    top_dotplot_genes_subset.to_csv(
        file_dir / f"{subset}_leiden_{res}_top_dotplot_genes_cluster_{cluster}.csv",
        index=False,
    )


# Annotate clusters based on marker genes and plot UMAP
annotation_dict = {
    "0": "Smooth muscle cells",
    "1": "Plasma cells",
    "2": "Lymphatic endothelial cells",
    "3": "Fibroblasts",
    "4": "T cells",
    "5": "Macrophages",
    "6": "Macrophages",
    "7": "Fibroblasts",
    "8": "Mast cells",
    "9": "Fibroblasts",
    "10": "Fibroblasts",
    "11": "Fibroblasts",
    "12": "Blood endothelial cells",
    "13": "Pericytes",
}

# Map Leiden clusters to annotation_level_1 and save in adata.obs
adata.obs[annotation_level_1] = (
    adata.obs[subset_key].map(annotation_dict).astype("category")
)

# Set colors for annotation_level_1 - match number of colors to number of categories
n_categories = len(adata.obs[annotation_level_1].cat.categories)
color_palette_annotation = sns.color_palette("hls", n_categories)
adata.uns[f"{annotation_level_1}_colors"] = [
    color for color in color_palette_annotation.as_hex()
]

# Plot UMAP colored by annotation_level_1
sc.pl.umap(
    adata,
    color=annotation_level_1,
    wspace=0.4,
    legend_fontsize=16,  # increase for larger/bolder appearance
    legend_fontoutline=2,  # white outline thickness
    show=False,
    frameon=False,
    save=f"_{annotation_level_1}.pdf",
)


# Spatial plot colored by annotation_level_1
ROI_list = adata.obs["ROI"].unique()
for roi in ROI_list:
    adata_roi = adata[adata.obs["ROI"] == roi].copy()

    cats = adata_roi.obs[annotation_level_1].cat.categories
    color_map = dict(
        zip(
            adata.obs[annotation_level_1].cat.categories,
            color_palette_annotation.as_hex(),
        )
    )

    adata_roi.uns[f"{annotation_level_1}_colors"] = [color_map[c] for c in cats]

    sq.pl.spatial_scatter(
        adata_roi,
        library_id="spatial",
        shape=None,
        color=annotation_level_1,
        size=2,
        marker=".",
        frameon=False,
        figsize=(10, 10),
        wspace=0.4,
        save=fig_dir_spatial / f"{roi}_{subset_key}.pdf",
    )


# Calculate number of cells per annotated cell type and save as csv
celltype_counts = pd.DataFrame(adata.obs[annotation_level_1].value_counts())
celltype_counts.to_csv(file_dir / f"{subset}_celltype_counts.csv")


# Number of cells per annotated cell type per condition and ROI
celltype_counts_condition = pd.crosstab(
    adata.obs[annotation_level_1], [adata.obs["condition"], adata.obs["ROI"]]
)
celltype_counts_condition.to_csv(
    file_dir / f"{subset}_celltype_counts_per_condition_ROI.csv"
)

print("Done!")
