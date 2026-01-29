"""Concatenate multiple AnnData objects into a single AnnData object."""

import logging
from pathlib import Path

import anndata as ad

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s - %(levelname)s] %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def concatenate_h5ad_files(adata_dir: str, output_path: str = None):
    """
    Load all .h5ad files from a directory and concatenate them.

    Parameters
    ----------
    adata_dir : str
        Path to directory containing individual .h5ad files
    output_path : str, optional
        Path where to save the concatenated file. If None, saves to adata_dir/all_samples.h5ad

    Returns
    -------
    ad.AnnData
        Concatenated AnnData object
    """
    adata_dir = Path(adata_dir)

    if not adata_dir.exists():
        raise FileNotFoundError(f"Directory not found: {adata_dir}")

    # Find all .h5ad files (excluding any existing all_samples.h5ad)
    h5ad_files = [f for f in adata_dir.glob("*.h5ad") if f.name != "all_samples.h5ad"]

    if not h5ad_files:
        raise ValueError(f"No .h5ad files found in {adata_dir}")

    logger.info(f"Found {len(h5ad_files)} .h5ad files to concatenate")

    # Load all AnnData objects
    all_adatas = []
    failed_files = []

    for h5ad_file in sorted(h5ad_files):
        try:
            logger.info(f"Loading {h5ad_file.name}...")
            adata = ad.read_h5ad(h5ad_file)

            # Verify .X is not None
            if adata.X is None:
                logger.error(f"{h5ad_file.name} has None for .X matrix, skipping")
                failed_files.append(h5ad_file.name)
                continue

            logger.info(f"  Shape: {adata.shape}, .X shape: {adata.X.shape}")
            all_adatas.append(adata)

        except Exception as e:
            logger.error(f"Failed to load {h5ad_file.name}: {e}")
            failed_files.append(h5ad_file.name)
            continue

    if not all_adatas:
        raise ValueError("No valid AnnData objects were loaded")

    logger.info(f"\nSuccessfully loaded {len(all_adatas)} AnnData objects")
    if failed_files:
        logger.warning(f"Failed to load {len(failed_files)} files: {failed_files}")

    # Concatenate
    logger.info("\nConcatenating AnnData objects...")
    logger.info("Using join='outer' to keep all genes...")

    try:
        combined = ad.concat(all_adatas, join="outer", fill_value=0, merge="unique")

        # Verify .X is not None after concatenation
        if combined.X is None:
            logger.warning("Combined AnnData has None for .X with join='outer'")
            logger.info("Trying join='inner' instead...")
            combined = ad.concat(all_adatas, join="inner", merge="unique")

            if combined.X is None:
                raise ValueError(
                    "Concatenation failed: .X is None even with join='inner'"
                )

        logger.info(f"\nCombined AnnData shape: {combined.shape}")
        logger.info(f"Combined .X shape: {combined.X.shape}")
        logger.info(f"Number of cells: {combined.n_obs}")
        logger.info(f"Number of genes: {combined.n_vars}")

        # Check for unique obs_names
        if not combined.obs_names.is_unique:
            logger.warning("obs_names are not unique after concatenation!")
            logger.info("Making obs_names unique...")
            combined.obs_names_make_unique()

        # Save concatenated AnnData
        if output_path is None:
            output_path = adata_dir / "all_samples.h5ad"
        else:
            output_path = Path(output_path)

        logger.info(f"\nSaving concatenated AnnData to {output_path}...")
        combined.write(output_path)
        logger.info("Done!")

        return combined

    except Exception as e:
        logger.error(f"Concatenation failed: {e}")
        raise


if __name__ == "__main__":
    # Directory containing individual .h5ad files
    adata_dir = "/Volumes/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/009_ST_Xenium/data/out/adata/"

    # Run concatenation
    try:
        combined_adata = concatenate_h5ad_files(adata_dir)
        print("\n" + "=" * 50)
        print("Concatenation completed successfully!")
        print("=" * 50)
    except Exception as e:
        print(f"\nError: {e}")
        raise
