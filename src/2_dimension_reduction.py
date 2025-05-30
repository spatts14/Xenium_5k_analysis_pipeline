# Import packages
import warnings  # ? what is the best way to suppress warnings from package inputs?

import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc
import torch

warnings.filterwarnings("ignore")

# from helper_function.py import seed_everything

# Set seed
# seed_everything(21122023)

# Set variables
module_name = "2_DR"

# Set directories
input_path = "/Users/sarapatti/Desktop/PhD_projects/Llyod_lab/ReCoDe-spatial-transcriptomics"
output_path = "/Users/sarapatti/Desktop/PhD_projects/Llyod_lab/ReCoDe-spatial-transcriptomics/analysis"
logging_path = "/Users/sarapatti/Desktop/PhD_projects/Llyod_lab/ReCoDe-spatial-transcriptomics/analysis/logging"


# Confirm directories exist
if not Path(input_path).exists():
    raise FileNotFoundError(f"Input path {input_path} does not exist.")
if not Path(output_path).exists():
    raise FileNotFoundError(f"Output path {output_path} does not exist.")


# Create output directories if they do not exist
os.makedirs(Path(output_path) / module_name, exist_ok=True)

# Set up logging
os.makedirs(
    logging_path, exist_ok=True
)  # should set up all these directories at the start of the pipeline?
logging.basicConfig(
    filename=Path(logging_path) / f"{module_name}.txt",  # output file
    filemode="w",  # overwrites the file each time
    format="%(asctime)s - %(levelname)s - %(message)s",  # log format
    level=logging.INFO,  # minimum level to log
)

# Import data
logging.info("Loading Xenium data...")
adata = sc.read_h5ad(Path(output_path) / "1_qc/adata.h5ad")

# Preform dimension reduction analysis
logging.info("compute PCA...")
sc.pp.pca(adata)  # compute principal components

logging.info("compute neighbors...")
sc.pp.neighbors(adata)  # compute a neighborhood graph

logging.info("create UMAPs and cluster cells..")
sc.tl.umap(adata)  # calculate umap
sc.tl.leiden(adata)  # determine cell clusters


# plot UMAP
sc.pl.umap(
    adata,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "leiden",
    ],
    wspace=0.4,
)
plt.tight_layout()
plt.savefig(Path(output_path) / f"{module_name}/UMAP.png", dpi=300)
plt.close()
logging.info(f"Saved plots to {module_name}/UMAP.png")

# Save anndata object
adata.write_h5ad(Path(output_path) / f"{module_name}/adata.h5ad")
logging.info(f"Data saved to {module_name}/adata.h5ad")
