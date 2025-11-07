from pathlib import Path
import scanpy as sc
import scipy.sparse as sp
import gc
import logging
import sys

# Set log files
log_file = "hlca_analysis.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),   # prints to console (stdout)
        logging.FileHandler(log_file)        # also saves to file
    ]
)

logging.info("Starting HLCA plotting script...")

# Set directories and file paths
output_dir = Path("/rds/general/user/sep22/home/Projects/_Public_datasets/HLCA/data/")
adata_path = output_dir / "HLCA_full.h5ad"
filtered_path = output_dir / "hlca_full_filtered.h5ad"
processed_path = output_dir / "hlca_full_processed.h5ad"

if not adata_path.exists():
    raise FileNotFoundError(f"{adata_path} not found.")

logging.info(f"Reading AnnData file from {adata_path}")
adata = sc.read_h5ad(adata_path, backed='r')

# Filter before loading into memory
mask = (
    adata.obs["disease"].isin([
        "normal",
        "pulmonary fibrosis",
        "interstitial lung disease",
        "chronic obstructive pulmonary disease"
    ])
)

adata_subset = adata[mask, :].to_memory()

logging.info(f"Diseases: {adata_subset.obs["disease"].value_counts()}")
logging.info(f"Tissue types: {adata_subset.obs["tissue_coarse_unharmonized"].value_counts()}")

adata = None  # free memory
gc.collect()

logging.info("Saving filtered data...")
adata_subset.write_h5ad(filtered_path)

# logging.info("Processing filtered data...")
# sc.pp.normalize_total(adata_subset, target_sum=1e4)
# sc.pp.log1p(adata_subset)
# sc.pp.highly_variable_genes(adata_subset, n_top_genes=2000)
# sc.tl.pca(adata_subset, n_comps=50)
# sc.pp.neighbors(adata_subset, n_neighbors=15, n_pcs=40)
# sc.tl.umap(adata_subset)

# logging.info("Saving filtered and processed data...")
# adata_subset.write_h5ad(processed_path)

logging.info("HLCA plotting script completed successfully.")
