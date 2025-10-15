#!/bin/bash
#PBS -l walltime=4:0:0
#PBS -lselect=1:ncpus=8:mem=64gb

# Load production tools
module load tools/prod
 
# Load python and bundle
module load Biopython/1.84-foss-2024a
 
# Change to directory
cd "$PBS_O_WORKDIR"

# Activate virtual environment
source xenium_5k_venv/bin/activate

# Run
python -m recode_st config_HPC.toml
