"""Level 1 cell type annotation for subset of cells in the Xenium data."""

import os
import random
from pathlib import Path

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


# Set random seed for reproducibility
seed_everything(19960915)

# Set variables
res = "leiden_res_0.5"

# Set directories
dir = Path(
    "/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/airscape_analysis/output/2026-02-18_17-40-00_analysis_run/"
)

fig_dir = dir / "manual_analysis/plots"
fig_dir.mkdir(parents=True, exist_ok=True)


# Set colors
cmap = sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True)
color_palette_level_1 = sns.color_palette("hls", 12)

# Load data
adata = sc.read_h5ad(
    dir / "output/2026-02-18_17-40-00_analysis_run/dimension_reduction/adata.h5ad"
)

sc.tl.umap(adata)

sc.pl.umap(
    adata,
    color=res,
    cmap=cmap,
    wspace=0.4,
    show=False,
    frameon=False,
    save=fig_dir / f"{res}.pdf",
)
