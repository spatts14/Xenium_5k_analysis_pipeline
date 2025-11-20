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

output_dir = Path("/rds/general/user/sep22/home/Projects/_Public_datasets/HLCA/data/")
fig_dir = Path(
    "/rds/general/user/sep22/home/Projects/_Public_datasets/HLCA/out/figures/"
)

# Set figure directory
sc.settings.figdir = fig_dir / "HLCA_full_figures"
sc.settings.figdir.mkdir(parents=True, exist_ok=True)
logging.info(f"Figures will be saved to {sc.settings.figdir}")

# Import data
adata_path = output_dir / "hlca_full_filtered.h5ad"
if not adata_path.exists():
    raise FileNotFoundError(f"{adata_path} not found.")

logging.info("Reading input AnnData file...")
adata = sc.read_h5ad(adata_path)

logging.info("Processing HLCA...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata)

logging.info("Plot UMAP by cell type...")
sc.pl.umap(
    adata,
    color="cell_type",
    title="HLCA full reference data UMAP",
    save="_cell_type.png",
)


adata.write_h5ad(output_dir / "hlca_full_processed.h5ad")
logging.info(
    f"Processed HLCA full and saved {output_dir / 'hlca_full_processed.h5ad'}."
)
