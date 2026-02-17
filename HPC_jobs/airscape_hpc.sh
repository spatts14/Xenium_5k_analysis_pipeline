#!/bin/bash
#PBS -l walltime=12:0:0
#PBS -lselect=1:ncpus=8:mem=512gb

# Load production tools
module load tools/prod

# Load python and bundle
module load Biopython/1.84-foss-2024a

# Change to directory
cd /rds/general/user/sep22/home/Projects/Xenium_5k_analysis_pipeline

# Activate virtual environment
source xenium_5k_venv/bin/activate

# Run
echo "Starting at annotation... $(date)"
python -m recode_st config_files/airscape_annotate.toml

# echo "Creating adata for each cell type... $(date)"
# python src/manual_src/celltype_subset.py

# echo "Subset each cell type... $(date)"
# python src/manual_src/celltype_level_1.py
