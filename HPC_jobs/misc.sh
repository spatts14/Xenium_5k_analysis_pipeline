#!/bin/bash
#PBS -l walltime=12:0:0
#PBS -lselect=1:ncpus=1:mem=256gb
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

python src/manual_src/misc_plots/16s.py

echo "Completed at $(date)"
