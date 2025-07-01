"""Muspan module - spatial graph construction."""

import warnings
from logging import getLogger

import matplotlib.pyplot as plt
import seaborn as sns

from recode_st.helper_function import seed_everything
from recode_st.logging_config import configure_logging
from recode_st.paths import output_path

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def run_muspan():
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
    module_name = "6_muspan"  # name of the module
    module_dir = output_path / module_name
    muspan_object = "muspan_object.muspan"
    seed = 21122023  # seed for reproducibility
    color_map = sns.color_palette("Blues", as_cmap=True)
    min_edge_distance = 0
    max_edge_distance = 45
    distance_list = [10, 20, 50]
    min_edge_distance_shape = 0
    max_edge_distance_shape = 1
    k_list = [2, 5, 10, 15]

    # Set seed
    seed_everything(seed)

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    print(module_dir)

    # Import data
    logger.info("Loading MuSpAn object...")
    domain = ms.io.load_domain(path_to_domain=str(module_dir / muspan_object))

    # Create delanuay triangulation spatial graph
    logger.info("Creating Delaunay triangulation spatial graph...")
    ms.networks.generate_network(
        domain,
        network_name="Delaunay CC filtered",
        network_type="Delaunay",
        objects_as_nodes=("collection", "Cell centroids"),
        min_edge_distance=min_edge_distance,
        max_edge_distance=max_edge_distance,
    )

    logger.info("Delaunay CC unfiltered:", domain.networks["Delaunay CC unfiltered"])
    logger.info("Delaunay CC filtered:", domain.networks["Delaunay CC filtered"])

    logger.info("Plotting filtered and unfiltered Delaunay networks...")
    # Plot the original Delaunay network
    fig, ax = plt.subplots(1, 2, figsize=(20, 12))
    ax[0].set_title("Delaunay CC unfiltered")
    ms.visualise.visualise_network(
        domain,
        network_name="Delaunay CC unfiltered",
        ax=ax[0],
        edge_weight_name=None,
        visualise_kwargs=dict(
            objects_to_plot=("collection", "Cell centroids"),
            marker_size=5,
            add_cbar=False,
            color_by=("constant", "#4D7EAB"),
        ),
    )

    # Plot the filtered Delaunay network
    ax[1].set_title("Delaunay CC filtered")
    ms.visualise.visualise_network(
        domain,
        network_name="Delaunay CC filtered",
        ax=ax[1],
        edge_weight_name=None,
        visualise_kwargs=dict(
            objects_to_plot=("collection", "Cell centroids"),
            marker_size=5,
            add_cbar=False,
            color_by=("constant", "#4D7EAB"),
        ),
    )

    plt.tight_layout()
    plt.savefig(module_dir / "muspan_delaunay.png")
    logger.info("Delaunay networks plotted and saved")

    # Create proximity triangulation spatial graph with point-like objects
    logger.info("Creating Proximity based networks (point-like objects)...")

    for distance in distance_list:
        ms.networks.generate_network(
            domain,
            network_name=f"prox network centroids {distance}",
            network_type="Proximity",
            objects_as_nodes=("collection", "Cell centroids"),
            max_edge_distance=distance,
            min_edge_distance=0,
        )

    logger.info("Plotting Proximity networks (point-like objects)...")
    # Plot the original Proximity network
    fig, ax = plt.subplots(1, 3, figsize=(20, 6))
    ax[0].set_title(f"Proximity network: {distance_list[0]} max distance")
    ms.visualise.visualise_network(
        domain,
        network_name=f"prox network centroids {distance_list[0]}",
        ax=ax[0],
        edge_cmap=color_map,
        edge_vmin=0,
        edge_vmax=distance_list[0],
        add_cbar=False,
        visualise_kwargs=dict(
            objects_to_plot=("collection", "Cell centroids"),
            marker_size=0.5,
            color_by=("constant", "black"),
        ),
    )

    # Plot the proximity network with 10μm max distance
    ax[1].set_title(f"Proximity network: {distance_list[1]} max distance")
    ms.visualise.visualise_network(
        domain,
        network_name=f"prox network centroids {distance_list[1]}",
        ax=ax[1],
        edge_cmap=color_map,
        edge_vmin=0,
        edge_vmax=distance_list[1],
        add_cbar=False,
        visualise_kwargs=dict(
            objects_to_plot=("collection", "Cell centroids"),
            marker_size=0.5,
            color_by=("constant", "black"),
        ),
    )

    # Plot the proximity network with 10μm max distance
    ax[2].set_title(f"Proximity network: {distance_list[2]} max distance")
    ms.visualise.visualise_network(
        domain,
        network_name=f"prox network centroids {distance_list[2]}",
        ax=ax[2],
        edge_cmap=color_map,
        edge_vmin=0,
        edge_vmax=distance_list[2],
        add_cbar=True,
        visualise_kwargs=dict(
            objects_to_plot=("collection", "Cell centroids"),
            marker_size=0.5,
            color_by=("constant", "black"),
        ),
    )

    plt.tight_layout()
    plt.savefig(module_dir / "muspan_proximity_point.png")
    logger.info("Proximity networks (point-like objects) plotted and saved")

    # Create proximity triangulation spatial graph with shape-like objects
    logger.info("Creating Proximity based networks (shape-like objects)...")

    ms.networks.generate_network(
        domain,
        network_name="Contact network",
        network_type="Proximity",
        objects_as_nodes=("collection", "Cell boundaries"),
        min_edge_distance=min_edge_distance_shape,
        max_edge_distance=max_edge_distance_shape,
    )

    logger.info("Plotting Proximity networks (shape-like objects)...")
    # Plot the cell boundaries underneath the network
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ms.visualise.visualise(
        domain,
        ax=ax,
        objects_to_plot=("collection", "Cell boundaries"),
        marker_size=10,
        add_cbar=False,
        shape_kwargs={"color": "#4D7EAB", "alpha": 0.4, "edgecolor": "black"},
    )

    # Plot the contact network on top of the cell boundaries
    ms.visualise.visualise_network(
        domain,
        network_name="Contact network",
        ax=ax,
        edge_weight_name=None,
        visualise_kwargs=dict(
            objects_to_plot=("collection", "Cell centroids"),
            marker_size=0.5,
            color_by=("constant", "black"),
        ),
    )

    plt.tight_layout()
    plt.savefig(module_dir / "muspan_proximity_shape.png")
    logger.info("Proximity networks (shape-like objects) plotted and saved")

    # Create proximity triangulation spatial graph with shape-like objects
    logger.info("Creating KNN based networks ...")

    for k in k_list:
        ms.networks.generate_network(
            domain,
            network_name=f"{k}-NN network",
            network_type="KNN",
            objects_as_nodes=("collection", "Cell centroids"),
            number_of_nearest_neighbours=k,
        )

    fig, axes = plt.subplots(1, len(k_list), figsize=(6 * len(k_list), 6))

    if len(k_list) == 1:
        axes = [axes]  # Ensure axes is iterable

    for i, k in enumerate(k_list):
        axes[i].set_title(f"{k}-NN network")
        ms.visualise.visualise_network(
            domain,
            network_name=f"{k}-NN network",
            ax=axes[i],
            edge_weight_name="Distance",
            edge_cmap=color_map,
            visualise_kwargs=dict(
                objects_to_plot=("collection", "Cell centroids"),
                marker_size=1,
                color_by=("constant", "black"),
            ),
        )

    plt.tight_layout()
    plt.savefig(module_dir / "muspan_knn.png")
    logger.info("KNN networks plotted and saved")

    # Confirm the domain has the expected labels
    logger.info(f"Networks in domain: {domain.networks.keys()}")

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

    try:
        run_muspan()
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
