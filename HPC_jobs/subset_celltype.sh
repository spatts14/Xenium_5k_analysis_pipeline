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

# Run with error logging
echo "Starting at $(date)"

# make new adata for each major cell type
python src/manual_src/celltype_subset.py


# high resolution clustering and plotting for each major cell type
export H5AD_FILE="adata_subset_Immune_cells.h5ad"
export SUBSET="immune"
python src/manual_src/celltype_level_1.py

# high resolution clustering and plotting for each major cell type
export H5AD_FILE="adata_subset_Airway_epithelial_cells.h5ad"
export SUBSET="airway_epithelial"
python src/manual_src/celltype_level_1.py

# high resolution clustering and plotting for each major cell type
export H5AD_FILE="adata_subset_Endothelial_cells.h5ad"
export SUBSET="endothelial"
python src/manual_src/celltype_level_1.py

# high resolution clustering and plotting for each major cell type
export H5AD_FILE="adata_subset_Stromal_cells.h5ad"
export SUBSET="stromal"
python src/manual_src/celltype_level_1.py

# high resolution clustering and plotting for each major cell type
export H5AD_FILE="adata_subset_Alveolar_epithelial_cells.h5ad"
export SUBSET="alveolar_epithelial"
python src/manual_src/celltype_level_1.py

echo "Completed at $(date)"
