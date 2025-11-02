"""Format for HLCA full data files."""

import scanpy as sc

# Load HLCA full dataset
adata = sc.read_h5ad(
    "/rds/general/user/sep22/home/Projects/_Public_datasets/HLCA/data/HLCA_full.h5ad"
)

print(adata.obs.columns)

# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
# sc.pp.highly_variable_genes(adata, n_top_genes=2000)
# sc.tl.pca(adata, n_comps=50)
# sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
# sc.tl.umap(adata)
# sc.pl.umap(adata, color="cell_type", title="HLCA full reference data UMAP")
