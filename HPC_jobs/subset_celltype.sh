#!/bin/bash
#PBS -l walltime=12:0:0
#PBS -lselect=1:ncpus=1:mem=128gb
#PBS -N subset_celltype
#PBS -j oe

# Load production tools
module load tools/prod

# Load python and bundle
module load Biopython/1.84-foss-2024a

# Change to directory
cd /rds/general/user/sep22/home/Projects/Xenium_5k_analysis_pipeline

# Activate virtual environment
source xenium_5k_venv/bin/activate

# Set directory paths
export ANNOTATE_DIR="/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/output/2026-02-22_analysis_run_HVG2000"
export CELLTYPE_SUBSET_DIR="/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/output/2026-02-22_analysis_run_HVG2000/celltype_subset"

# Run with error logging
echo "Starting at $(date)"

# # make new adata for each major cell type
# python src/manual_src/celltype_subset.py

# # high resolution clustering and plotting for each major cell type
# export H5AD_FILE="adata_subset_Immune_cells.h5ad"
# export SUBSET="Immune_cells"
# python src/manual_src/celltype_level_1.py

# high resolution clustering and plotting for each major cell type
export H5AD_FILE="adata_subset_Airway_epithelial_cells.h5ad"
export SUBSET="Airway_epithelial_cells"
python src/manual_src/celltype_level_1.py

# # high resolution clustering and plotting for each major cell type
# export H5AD_FILE="adata_subset_Endothelial_cells.h5ad"
# export SUBSET="Endothelial_cells"
# python src/manual_src/celltype_level_1.py

# # high resolution clustering and plotting for each major cell type
# export H5AD_FILE="adata_subset_Stromal_cells.h5ad"
# export SUBSET="Stromal_cells"
# python src/manual_src/celltype_level_1.py

# # high resolution clustering and plotting for each major cell type
# export H5AD_FILE="adata_subset_Alveolar_epithelial_cells.h5ad"
# export SUBSET="Alveolar_epithelial_cells"
# python src/manual_src/celltype_level_1.py

echo "Completed at $(date)"
