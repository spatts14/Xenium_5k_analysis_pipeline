"""Muspan module."""

import warnings
from logging import getLogger

import muspan as ms
import scanpy as sc

from recode_st.helper_function import seed_everything
from recode_st.logging_config import configure_logging
from recode_st.paths import area_path, output_path, xenium_path

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def run_muspan():
    """Run Muspan on Xenium data."""
    try:
        import muspan  # noqa: F401
    except ModuleNotFoundError as err:
        logger.error(
            "Could not load necessary MuSpAn package. You can obtain this with:\n"
            "    pip install 'recode_st[muspan] @ git+"
            "https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics.git"
        )
        raise err
    # Set variables
    module_name = "6_muspan"  # name of the module
    module_dir = output_path / module_name
    seed = 21122023  # seed for reproducibility
    domain_name = "Xenium_lung_FFPE"  # name of the domain
    transcripts_of_interest = [
        "EPCAM",
        "CD3D",
        "CD68",
        "VWF",
        "PTPRC",
        "ACTA2",
    ]

    # Set seed
    seed_everything(seed)

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Import data
    logger.info("Loading Xenium data...")
    adata = sc.read_h5ad(output_path / "5_spatial_stats" / "adata.h5ad")

    # Import Xenium data using muspan
    tofi = transcripts_of_interest

    # Create muspan object
    logger.info("Creating MuSpAn domain object...")
    domain = ms.io.xenium_to_domain(
        path_to_xenium_data=str(xenium_path),  # Convert Path to string
        cells_from_selection_csv=str(area_path),  # Convert Path to string
        domain_name=domain_name,
        load_transcripts=True,
        selected_transcripts=tofi,
        load_nuclei=True,
        load_cells_as_shapes=True,
        exclude_no_nuclei_cells=True,
    )
    logger.info(f"Domain object created: {domain}")

    # Filer adata only to include cells within the area of interest
    logger.info("Filtering adata to include only cells in the area of interest...")
    domain_cell_ids = domain.labels["Cell ID"]  # Get cell IDs from the domain
    filt_adata = adata[
        adata.obs_names.isin(domain_cell_ids)
    ].copy()  # filter adata object
    logger.info(f"Filtered adata from {adata.n_obs} to {filt_adata.n_obs} cells")

    # Add cell cluster IDs
    logger.info("All cell_type IDs to domain with cluster labels")
    cell_type = filt_adata.obs["cell_type"]  # filter on cell type column
    domain.add_labels(label_name="cell_type", labels=cell_type)
    logger.info(f"Label keys in domain: {domain.labels.keys()}")


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.6_muspan")  # re-name the logger to match the module

    try:
        run_muspan()
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
