#!/bin/bash
#PBS -l walltime=12:0:0
#PBS -lselect=1:ncpus=8:mem=512gb
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
# python src/manual_src/celltype_subset.py
python src/manual_src/celltype_level_1.py
echo "Completed at $(date)"
