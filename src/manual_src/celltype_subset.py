"""Module to subset the AnnData object to include only specified cell types."""

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

# Set directories
dir = Path(
    "/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/output/2026-02-21_analysis_run_HVG500/"
)

output_dir = (
    Path(
        "/rds/general/user/sep22/projects/phenotypingsputumasthmaticsaurorawellcomea1/ephemeral"
    )
    / "celltype_subset"
)
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
