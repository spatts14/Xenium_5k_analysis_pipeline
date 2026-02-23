"""Module to subset the AnnData object to include only specified cell types."""

import os
from pathlib import Path

import scanpy as sc


def subset_cell_types(adata, celltype_col: str, output_dir: Path):
    """Subset the AnnData object to include only specified cell types.

    Args:
        adata (AnnData): The input AnnData object containing data.
        celltype_col (str): Name of the column in adata.obs that contains annotations.
        output_dir (Path): The output directory where subsetted objects will be saved.

    Returns:
        AnnData: A new AnnData object containing only the specified cell types.
    """
    # Check if the specified cell types are present in adata.obs[celltype_col]
    celltype_list = adata.obs[celltype_col].unique().tolist()
    print(f"Available cell types in adata.obs['{celltype_col}']: {celltype_list}")

    # Subset the AnnData object
    for celltype in celltype_list:
        print(f"Subsetting for cell type: {celltype}")
        if celltype not in adata.obs[celltype_col].unique():
            raise ValueError(
                f"Cell type '{celltype}' not found in adata.obs['{celltype_col}']"
            )
        subset_adata = adata[adata.obs[celltype_col] == celltype].copy()

        # Make filename safe
        celltype_safe = celltype.replace(" ", "_").replace("/", "_")

        # Save the subsetted AnnData object
        subset_adata.write_h5ad(output_dir / f"adata_subset_{celltype_safe}.h5ad")
        print(f"Subsetted AnnData object saved to adata_subset_{celltype_safe}.h5ad")


manual_annotation = "level_0_annotation"

# Read directories from environment variables
annotate_dir = os.getenv("ANNOTATE_DIR")
if not annotate_dir:
    raise ValueError("ANNOTATE_DIR environment variable must be set")
dir = Path(annotate_dir)

celltype_subset_dir = os.getenv("CELLTYPE_SUBSET_DIR")
if not celltype_subset_dir:
    raise ValueError("CELLTYPE_SUBSET_DIR environment variable must be set")
output_dir = Path(celltype_subset_dir)
output_dir.mkdir(parents=True, exist_ok=True)


# Read in adata file
adata = sc.read_h5ad(dir / "annotate/adata.h5ad")

# Subset and save subsetted adata
subset_cell_types(adata, celltype_col=manual_annotation, output_dir=output_dir)

# List all the saved files in the module directory
print("\nSaved subsetted AnnData files:")
for file in output_dir.glob("adata_subset_*.h5ad"):
    print(file.name)

print("Sub-setting major cell types module completed successfully.")
