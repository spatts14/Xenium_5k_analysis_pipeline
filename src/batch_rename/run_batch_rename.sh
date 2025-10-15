#!/bin/bash

# Batch rename script for Xenium morphology_focus files
# This script will process all output directories in your Xenium run

# Set your data directory path
DATA_DIR="/Volumes/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/data/xenium_raw/20251001__141239__SP25164_SARA_PATTI_RUN_1"

# Set the mapping file path
MAPPING_FILE="xenium_mapping.md"

echo "Xenium Batch Rename Script"
echo "=========================="
echo "Data directory: $DATA_DIR"
echo "Mapping file: $MAPPING_FILE"
echo ""

# Check if mapping file exists
if [ ! -f "$MAPPING_FILE" ]; then
    echo "Error: Mapping file '$MAPPING_FILE' not found!"
    echo "Please make sure the mapping file is in the current directory."
    exit 1
fi

# First, run a dry-run to see what will be changed
echo "Running dry-run to preview changes..."
echo "------------------------------------"
python3 xenium_batch_rename.py -d "$DATA_DIR" -m "$MAPPING_FILE" --dry-run

echo ""
echo "Dry-run completed. Review the output above."
echo ""
read -p "Do you want to proceed with the actual renaming? (y/N): " -n 1 -r
echo ""

if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Proceeding with renaming..."
    echo "---------------------------"
    python3 xenium_batch_rename.py -d "$DATA_DIR" -m "$MAPPING_FILE" --rename --backup
    echo ""
    echo "Batch rename completed!"
else
    echo "Operation cancelled."
fi
