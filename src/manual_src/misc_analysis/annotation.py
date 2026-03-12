from pathlib import Path

import pandas as pd
import scanpy as sc


def collect_all_annotations(base_dir):
    """Traverse directory structure to find all *_annotations.csv files.

    Args:
        base_dir: Base directory containing subdirectories with annotation files

    Returns:
        pd.DataFrame: Combined annotations from all CSV files
    """
    base_path = Path(base_dir)

    if not base_path.exists():
        raise FileNotFoundError(f"Base directory not found: {base_dir}")

    all_annotations = []

    # Find all subdirectories in base_dir
    for subdir in base_path.iterdir():
        if subdir.is_dir():
            # Look for 'files' directory within each subdirectory
            files_dir = subdir / "files"

            if files_dir.exists() and files_dir.is_dir():
                # Find all *_annotations.csv files
                annotation_files = list(files_dir.glob("*_annotations.csv"))

                for csv_file in annotation_files:
                    print(f"Found: {csv_file}")

                    try:
                        # Read the CSV file
                        df = pd.read_csv(csv_file)

                        # Add source information
                        df["source_subdir"] = subdir.name
                        df["source_file"] = csv_file.name

                        all_annotations.append(df)
                        print(f"Loaded {len(df)} annotations from {subdir.name}")

                    except Exception as e:
                        print(f"Error reading {csv_file}: {e}")

    if not all_annotations:
        raise ValueError(f"No *_annotations.csv files found in {base_dir}")

    # Combine all annotations
    combined = pd.concat(all_annotations, ignore_index=True)

    # Rename first column
    combined.rename(columns={combined.columns[0]: "annotation"}, inplace=True)

    print(f"\nTotal annotations loaded: {len(combined)}")
    print(f"Columns: {combined.columns.tolist()}")

    return combined


def rename_cells_from_annotations(
    adata, annotations_df, cell_id_col="cell_id", annotation_col="cell_annotation"
):
    """Rename cells in AnnData object using annotations from CSV.

    Args:
        adata: AnnData object
        annotations_df: DataFrame with cell IDs and annotations
        cell_id_col: Column name in annotations_df containing cell IDs (default: 'cell_id')
        annotation_col: Column name in annotations_df containing annotations (default: 'cell_annotation')

    Returns:
        AnnData object with updated cell annotations
    """
    # Check if columns exist
    if cell_id_col not in annotations_df.columns:
        raise ValueError(
            f"Column '{cell_id_col}' not found in annotations. Available: {annotations_df.columns.tolist()}"
        )

    if annotation_col not in annotations_df.columns:
        raise ValueError(
            f"Column '{annotation_col}' not found in annotations. Available: {annotations_df.columns.tolist()}"
        )

    # Check if cell_id exists in adata.obs
    if "cell_id" not in adata.obs.columns:
        print(
            "Warning: 'cell_id' not found in adata.obs, using adata.obs_names instead"
        )

    # Create mapping dictionary
    cell_id_to_annotation = dict(
        zip(annotations_df[cell_id_col], annotations_df[annotation_col])
    )

    print(f"\nMapping {len(cell_id_to_annotation)} cell IDs to annotations")

    # Map annotations to adata
    adata.obs["cell_annotation"] = adata.obs["cell_id"].map(cell_id_to_annotation)

    # Check for unmapped cells
    unmapped = adata.obs["cell_annotation"].isna().sum()
    if unmapped > 0:
        print(f"Warning: {unmapped} cells could not be mapped to annotations")
        print(f"Total cells in adata: {len(adata)}")
        print(f"Successfully mapped: {len(adata) - unmapped}")
    else:
        print(f"Successfully mapped all {len(adata)} cells!")

    # Print summary
    print("\nAnnotation summary:")
    print(adata.obs["cell_annotation"].value_counts())

    return adata


# Define base directory
base_dir = Path(
    "/Volumes/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/output/2026-02-22_analysis_run_HVG2000/celltype_subset/"
)

# Combine all annotation files and save
annotations_df = collect_all_annotations(base_dir)
output_path = Path(base_dir) / "combined_annotations.csv"
annotations_df.to_csv(output_path, index=False)
print(f"Saved combined annotations to: {output_path}")

# Load adata and rename cells based on annotations
adata = rename_cells_from_annotations(
    adata, annotations_df, cell_id_col=cell_id_col, annotation_col=annotation_col
)

# Verify results
print("\n" + "=" * 80)
print("VERIFICATION")
print("=" * 80)
print(f"Total cells: {len(adata)}")
print(f"Cells with annotations: {adata.obs['cell_annotation'].notna().sum()}")
print(f"\nUnique annotations: {adata.obs['cell_annotation'].nunique()}")
