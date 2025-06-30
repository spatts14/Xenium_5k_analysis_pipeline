"""Muspan module."""

import warnings
from logging import getLogger

import matplotlib.pyplot as plt
import scanpy as sc

from recode_st.config import MuspanModuleConfig
from recode_st.helper_function import seed_everything
from recode_st.logging_config import configure_logging
from recode_st.paths import area_path, output_path, xenium_path

warnings.filterwarnings("ignore")

logger = getLogger(__name__)

# TODO: replace "cell_type" with variable


def run_muspan(config: MuspanModuleConfig):
    """Run Muspan on Xenium data."""
    try:
        import muspan as ms
    except ModuleNotFoundError as err:
        logger.error(
            "Could not load necessary MuSpAn package. You can obtain this with:\n"
            "    pip install 'recode_st[muspan] @ git+"
            "https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics.git"
        )
        raise err
    # Set variables
    module_dir = output_path / config.module_name
    domain_name = "Xenium_lung_FFPE"  # name of the domain
    transcripts_of_interest = [
        "EPCAM",
        "CD3D",
        "CD68",
        "VWF",
        "PTPRC",
        "ACTA2",
    ]

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
    # Get cell IDs in their original order (preserving duplicates if any)
    domain_cell_ids_ordered = [
        str(cell_id) for cell_id in domain.labels["Cell ID"]["labels"]
    ]

    # Get unique cell IDs for filtering adata
    domain_cell_ids_unique = set(domain_cell_ids_ordered)

    logger.info(f"Number of unique cells in the domain: {len(domain_cell_ids_unique)}")
    logger.info(
        f"Total cell entries in domain (including duplicates): "
        f"{len(domain_cell_ids_ordered)}"
    )

    # Filter adata to include only cells in the area of interest
    filt_adata = adata[
        adata.obs["cell_id"].isin(domain_cell_ids_unique)
    ]  # filter adata object

    logger.info(f"Filtered adata from {adata.n_obs} to {filt_adata.n_obs} cells")

    # Add cell cluster IDs
    logger.info("Adding cell_type IDs to domain with cluster labels")

    # Create a mapping from cell_id to cell_type
    cell_id_to_type = dict(zip(filt_adata.obs["cell_id"], filt_adata.obs["cell_type"]))

    # Get cell types in the same order as domain cell IDs
    # (preserving order and duplicates)
    cell_types_ordered = [
        cell_id_to_type.get(cell_id, "Unknown") for cell_id in domain_cell_ids_ordered
    ]

    # Add cell_type label to the domain
    domain.add_labels(label_name="cell_type", labels=cell_types_ordered)
    logger.info(f"Label keys in domain: {domain.labels.keys()}")
    logger.info(f"Length of cell_types_ordered: {len(cell_types_ordered)}")
    logger.info(f"Length of domain cell IDs: {len(domain_cell_ids_ordered)}")
    logger.info(
        f"Number of 'Unknown' cell types: {cell_types_ordered.count('Unknown')}"
    )

    logger.info("Queries to isolate the different objects within the MuSpAn domain")

    qCells = ms.query.query(
        domain, ("Collection",), "is", "Cell boundaries"
    )  # Query to isolate cell boundaries

    qTrans = ms.query.query(
        domain, ("Collection",), "is", "Transcripts"
    )  # Query to isolate transcripts
    qNuc = ms.query.query(
        domain, ("Collection",), "is", "Nucleus boundaries"
    )  # Query to isolate nucleus boundaries

    logger.info(f"Visualize the MuSpAn domain: {domain_name}")
    fig, ax = plt.subplots(figsize=(20, 15), nrows=2, ncols=2)
    # Visualise all objects in the MuSpAn domain
    ms.visualise.visualise(domain, ax=ax[0, 0], marker_size=0.05)
    ax[0, 0].set_title("All objects")

    ms.visualise.visualise(
        domain,
        color_by=("label", "cell_type"),
        ax=ax[0, 1],
        objects_to_plot=qCells,
    )
    ax[0, 1].set_title("Cell type")

    # Visualise transcripts, colored by 'Transcript'
    ms.visualise.visualise(
        domain,
        color_by=("label", "Transcript ID"),
        ax=ax[1, 0],
        objects_to_plot=qTrans,
        marker_size=5,
    )
    ax[1, 0].set_title("Transcripts")

    # Visualise nuclei, colored by 'Nucleus Area'
    ms.visualise.visualise(
        domain,
        color_by=("label", "Nucleus Area"),
        ax=ax[1, 1],
        objects_to_plot=qNuc,
    )
    ax[1, 1].set_title("Nuclei")

    plt.tight_layout()
    plt.savefig(module_dir / "muspan_domain_visualization.png")
    logger.info("Muspan Domain Visualization plotted and saved")

    logger.info("Convert cell boundaries to cell centres (centroids)")
    domain.convert_objects(
        population=("Collection", "Cell boundaries"),
        object_type="point",
        conversion_method="centroids",
        collection_name="Cell centroids",
        inherit_collections=False,
    )

    # Plot the cell centres with color based on 'cell_type'
    plt.figure(figsize=(10, 6))
    ms.visualise.visualise(
        domain,
        objects_to_plot=("collection", "Cell centroids"),
        color_by="cell_type",
        ax=plt.gca(),
    )
    plt.tight_layout()
    plt.savefig(module_dir / "muspan_cell_centroids.png")
    logger.info(
        "Cell centroids visualized with cell_type color coding plotted and saved"
    )

    # Overlay the cell centres with cell type
    plt.figure(figsize=(10, 6))
    ms.visualise.visualise(
        domain, objects_to_plot=qCells, color_by="cell_type", ax=plt.gca()
    )
    ms.visualise.visualise(
        domain,
        objects_to_plot=("collection", "Cell centroids"),
        color_by=("constant", "black"),
        ax=plt.gca(),
        marker_size=1,
        add_cbar=False,
    )
    plt.tight_layout()
    plt.savefig(module_dir / "muspan_cell_centroids_n_boundaries.png")
    logger.info(
        "Cell centroids & cell boundaries with cell type color plotted and saved"
    )

    # Save domain
    ms.io.save_domain(
        domain,
        name_of_file="muspan_object",
        path_to_save=str(module_dir),
    )
    logger.info("Domain saved")


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.6_muspan")  # re-name the logger to match the module

    # Set seed
    seed_everything(21122023)

    try:
        run_muspan(MuspanModuleConfig(module_name="6_muspan"))
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
