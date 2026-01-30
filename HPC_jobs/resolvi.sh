#!/bin/bash
#PBS -q gpu72
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=8:mem=128gb:ngpus=1
#PBS -j oe

# Load production tools
module load tools/prod

# Load python and bundles
module load CUDA/12.1
module load Biopython/1.84-foss-2024a

# Change to directory
cd /rds/general/user/sep22/home/Projects/Xenium_5k_analysis_pipeline

# Activate virtual environment
source xenium_5k_venv/bin/activate

# Run
python -m recode_st config_files/resolvi.toml

