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

adata_path = output_dir / "HLCA_full.h5ad"
if not adata_path.exists():
    raise FileNotFoundError(f"{adata_path} not found.")


logging.info("Reading input AnnData file...")
adata = sc.read_h5ad(adata_path)


for key in ["disease", "tissue_coarse_unharmonized", "cell_type"]:
    if key not in adata.obs:
        raise KeyError(f"{key} not found in adata.obs columns.")

# Subset disease types
adata = adata[
    adata.obs["disease"].isin(
        [
            "normal",
            "pulmonary fibrosis",
            "interstitial lung disease",
            "chronic obstructive pulmonary disease",
        ]
    )
].copy()

# Subset tissue types
adata = adata[
    adata.obs["tissue_coarse_unharmonized"].isin(
        [
            "parenchyma",
            "airway",
            "Intermediate Bronchi",
            "Distal Bronchi",
            "segmental_bronchi",
        ]
    )
].copy()

adata.write_h5ad(output_dir / "hlca_full_filtered.h5ad")

# Process
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata)

sc.pl.umap(adata, color="cell_type", title="HLCA full reference data UMAP")
adata.write_h5ad(output_dir / "hlca_full_processed.h5ad")
