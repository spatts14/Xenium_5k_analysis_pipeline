import logging
import sys
from pathlib import Path

import scanpy as sc

log_file = "hlca_analysis.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),  # prints to console (stdout)
        logging.FileHandler(log_file),  # also saves to file
    ],
)

logging.info("Starting HLCA plotting script...")

base_dir = Path("/rds/general/user/sep22/home/Projects/_Public_datasets/HLCA/data/")
fig_dir = Path(
    "/rds/general/user/sep22/home/Projects/_Public_datasets/HLCA/out/figures/"
)

adata_path = base_dir / "hlca_full_processed.h5ad"
if not adata_path.exists():
    raise FileNotFoundError(f"{adata_path} not found.")

logging.info(f"Reading AnnData file from {adata_path}")
adata = sc.read_h5ad(base_dir / "hlca_full_processed.h5ad")

# Remove nose cells
remove = ["nose", "inferior turbinate", "trachea"]
adata = adata[~adata.obs["tissue_level_2"].isin(remove), :].copy()
logging.info(
    f"Filtered out nose cells. Remaining tissue types: {adata.obs['tissue_level_2'].unique()}"
)

# Set figure directory
sc.settings.figdir = fig_dir / "HLCA_full_figures_tissue_filtered"
sc.settings.figdir.mkdir(parents=True, exist_ok=True)
logging.info(f"Figures will be saved to {sc.settings.figdir}")

logging.info("Generating UMAP: cell_type")
sc.pl.umap(
    adata,
    color="cell_type",
    title="HLCA full reference data: cell type",
    save="_cell_type.png",
)

logging.info("Generating UMAP: tissue_coarse_unharmonized")
sc.pl.umap(
    adata,
    color="tissue_coarse_unharmonized",
    title="HLCA full reference data: tissue",
    save="_tissue_coarse_unharmonized.png",
)

logging.info("Generating UMAP: tissue_level_2")
sc.pl.umap(
    adata,
    color="tissue_level_2",
    title="HLCA full reference data: tissue",
    save="_tissue_level_2.png",
)

logging.info("Generating UMAP: disease")
sc.pl.umap(
    adata,
    color="disease",
    title="HLCA full reference data: disease status",
    save="_disease.png",
)

logging.info("All plots generated successfully.")
