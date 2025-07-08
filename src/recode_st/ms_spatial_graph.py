"""Muspan module - spatial graph construction."""

import warnings
from logging import getLogger

import matplotlib.pyplot as plt
import seaborn as sns

from recode_st.config import IOConfig, MuspanSpatialGraphModuleConfig
from recode_st.helper_function import seed_everything
from recode_st.logging_config import configure_logging

warnings.filterwarnings("ignore")

logger = getLogger(__name__)

# Import muspan at module level
try:
    import muspan as ms
except ModuleNotFoundError:
    ms = None


def run_muspan_graph(config: MuspanSpatialGraphModuleConfig, io_config: IOConfig):
    """Run Muspan spatial graph analysis on Xenium data."""
    if ms is None:
        logger.error(
            "Could not load necessary MuSpAn package. You can obtain this with:\n"
            "    pip install 'recode_st[muspan] @ git+"
            "https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics.git"
        )
        raise ModuleNotFoundError("MuSpAn package not found")

    # Set variables from config
    module_dir = io_config.output_dir / config.module_name
    muspan_object = config.muspan_object
    min_edge_distance = config.min_edge_distance
    max_edge_distance = config.max_edge_distance
    distance_list = config.distance_list
    min_edge_distance_shape = config.min_edge_distance_shape
    max_edge_distance_shape = config.max_edge_distance_shape
    k_list = config.k_list
    color_map = sns.color_palette("Blues", as_cmap=True)

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Import data
    logger.info("Loading MuSpAn object...")
    domain = ms.io.load_domain(path_to_domain=str(module_dir / muspan_object))

    # Create delanuay triangulation spatial graph
    logger.info("Creating Delaunay triangulation spatial graph...")
    create_delaunay_unfiltered_network(domain)
    create_delaunay_filt_network(domain, min_edge_distance, max_edge_distance)
    logger.info("Plotting filtered and unfiltered Delaunay networks...")
    plot_delaunay_networks(domain, module_dir)

    # Create proximity triangulation spatial graph with point-like objects
    logger.info("Creating Proximity based networks (point-like objects)...")
    create_proximity_point_networks(domain, distance_list)
    logger.info("Plotting Proximity networks (point-like objects)...")
    plot_proximity_networks(domain, module_dir, color_map, distance_list)

    # Create proximity triangulation spatial graph with shape-like objects
    logger.info("Creating Proximity based networks (shape-like objects)...")
    create_proximity_shape(domain, min_edge_distance_shape, max_edge_distance_shape)
    logger.info("Plotting Proximity networks (shape-like objects)...")
    plot_proximity_shape(domain, module_dir)

    # Create KNN based networks
    logger.info("Creating KNN based networks ...")
    create_knn_networks(domain, k_list)
    logger.info("Plotting KNN based networks...")
    plot_knn_networks(domain, module_dir, color_map, k_list)

    # Confirm the domain has the expected labels
    logger.info(f"Networks in domain: {domain.networks.keys()}")

    # Save domain
    ms.io.save_domain(
        domain,
        name_of_file="muspan_object",
        path_to_save=str(module_dir),
    )
    logger.info("Domain saved")


def plot_delaunay_networks(domain, module_dir):
    """Plot the unfiltered and filtered Delaunay cell-cell (CC) networks for domain.

    This function creates a side-by-side visualization of two Delaunay networks:
    one unfiltered and one filtered, using the provided domain object. The networks
    are visualized with cell centroids as markers and edges colored in a constant color.

    Args:
        domain (muspan object): The spatial domain object containing network
            information and cell centroids. Must be compatible with
            `ms.visualise.visualise_network`.
        module_dir (Path or str): Directory path where the output image will be saved.

    Returns:
        None: Displays the generated matplotlib figure with the two network plots.
    """
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


def create_delaunay_filt_network(domain, min_edge_distance, max_edge_distance):
    """Creates a filtered Delaunay network for a given spatial domain.

    This function generates a Delaunay triangulation-based network using the domain,
    where nodes represent cell centroids. Edges are filtered based on the provided
    minimum and maximum edge distance thresholds.

    Args:
        domain (muspan object):
            The spatial domain or object collection on which to construct the network.
        min_edge_distance (float): The minimum distance for edges in the network.
        max_edge_distance (float): The maximum distance for edges in the network.

    Returns:
        None. The network is generated and stored within the provided domain context.

    Raises:
        Any exceptions raised by `ms.networks.generate_network` will propagate.
    """
    ms.networks.generate_network(
        domain,
        network_name="Delaunay CC filtered",
        network_type="Delaunay",
        objects_as_nodes=("collection", "Cell centroids"),
        min_edge_distance=min_edge_distance,
        max_edge_distance=max_edge_distance,
    )


def create_delaunay_unfiltered_network(domain):
    """Creates an unfiltered Delaunay network for the given domain.

    This function generates a Delaunay network using the cell centroids
    from the specified domain. The resulting network is named "Delaunay CC unfiltered"
    and is not filtered by any criteria.

    Args:
        domain: (muspan object)
            The spatial domain or dataset containing cell centroid information.

    Returns:
        None. The function creates the network as a side effect.
    """
    ms.networks.generate_network(
        domain,
        network_name="Delaunay CC unfiltered",
        network_type="Delaunay",
        objects_as_nodes=("collection", "Cell centroids"),
    )


def create_proximity_shape(domain, min_edge_distance_shape, max_edge_distance_shape):
    """Creates a proximity-based contact network for a given domain.

    Args:
        domain: (muspan object)
            The spatial domain or context in which the network will be generated.
        min_edge_distance_shape: float
            The minimum distance between cell boundaries to consider
            an edge in the network.
        max_edge_distance_shape: float
            The maximum distance between cell boundaries to consider
            an edge in the network.

    Returns:
        None

    Side Effects:
        Modifies the given domain by adding a "Contact network" of type "Proximity"
        based on the specified distance thresholds.
    """
    ms.networks.generate_network(
        domain,
        network_name="Contact network",
        network_type="Proximity",
        objects_as_nodes=("collection", "Cell boundaries"),
        min_edge_distance=min_edge_distance_shape,
        max_edge_distance=max_edge_distance_shape,
    )


def create_proximity_point_networks(domain, distance_list):
    """Creates proximity-based point networks for specified distances.

    For each distance in `distance_list`, this function generates a proximity network
    using the centroids of cells as nodes. The network is created within the
    specified `domain`and includes edges between nodes whose centroids are within
    the given maximum distance.

    Args:
        domain: (muspan object)
            The domain or region within which to generate the networks.
        distance_list (list of float): List of maximum edge distances for which to
        generate proximity networks.

    Returns:
        Calculates and saves proximity networks for each distance in `distance_list`.
    """
    for distance in distance_list:
        ms.networks.generate_network(
            domain,
            network_name=f"prox network centroids {distance}",
            network_type="Proximity",
            objects_as_nodes=("collection", "Cell centroids"),
            max_edge_distance=distance,
            min_edge_distance=0,
        )


def plot_proximity_shape(domain, module_dir):
    """Plots the proximity network overlaid on cell boundaries for the domain.

    This function visualizes the cell boundaries and the contact
    network (proximity network) for the provided spatial domain object.
    The resulting plot is saved as a PNG file in the specified module directory.

    Args:
        domain: (muspan object)
            The spatial domain object containing cell boundary and network information.
        module_dir (Path or str): Directory path where the output image will be saved.

    Side Effects:
        - Saves PNG image named 'muspan_proximity_shape.png' in the specified directory.
        - Logs an info message upon successful plot creation.
    """
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


def plot_knn_networks(domain, module_dir, color_map, k_list):
    """Plots and saves k-nearest neighbor (k-NN) networks for a given spatial domain.

    For each value of k in `k_list`, this function visualizes the corresponding
    k-NN network using the provided color map for edge weights and saves the
    resulting figure as 'muspan_knn.png' in the specified module directory.

    Args:
        domain: (muspan object)
            The spatial domain object containing network and cell centroid information.
        module_dir (Path or str): Directory path where the output image will be saved.
        color_map: Colormap to use for visualizing edge weights in the network.
        k_list (list of int): List of k values for which to plot k-NN networks.

    Returns:
        None

    Side Effects:
        - Saves the generated plot as 'muspan_knn.png' in `module_dir`.
        - Logs an info message upon successful plot saving.
    """
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


def create_knn_networks(domain, k_list):
    """Creates multiple k-nearest neighbor (KNN) networks the domain using a list of k.

    Args:
        domain: (muspan object)
            The spatial domain or dataset on which to generate the KNN networks.
        k_list (list of int): A list of integers specifying the number of nearest
        neighbors (k) for each network to be created.

    Returns:
        None

    Side Effects:
        Modifies the given domain by adding KNN networks with the specified k values.
        Each network is named as "{k}-NN network".
    """
    for k in k_list:
        ms.networks.generate_network(
            domain,
            network_name=f"{k}-NN network",
            network_type="KNN",
            objects_as_nodes=("collection", "Cell centroids"),
            number_of_nearest_neighbours=k,
        )


def plot_proximity_networks(domain, module_dir, color_map, distance_list):
    """Plot proximity networks for point-like objects at specified distances.

    This function generates and saves visualizations of proximity networks the domain,
    using a list of maximum distances. Each subplot corresponds to a different
    proximity threshold, displaying the network of connections between
    point-like objects (e.g., cell centroids) within the specified distance.

    Args:
        domain: (muspan object)
            The spatial domain containing point-like objects to be analyzed.
        module_dir (Path or str): Directory path where the output plot image is saved.
        color_map: Colormap used for visualizing edge distances in the network.
        distance_list (list of float): List of maximum distances to
        define proximity networks.

    Returns:
        None. The function saves the generated plot as 'muspan_proximity_point.png'
        in the specified directory and logs the completion of the plotting process.
    """
    fig, axes = plt.subplots(1, len(distance_list), figsize=(20, 6))

    if len(distance_list) == 1:
        axes = [axes]  # Ensure axes is iterable

    for i, distance in enumerate(distance_list):
        axes[i].set_title(f"Proximity network: {distance} max distance")
        ms.visualise.visualise_network(
            domain,
            network_name=f"prox network centroids {distance}",
            ax=axes[i],
            edge_cmap=color_map,
            edge_vmin=0,
            edge_vmax=distance,
            add_cbar=False,
            visualise_kwargs=dict(
                objects_to_plot=("collection", "Cell centroids"),
                marker_size=0.5,
                color_by=("constant", "black"),
            ),
        )
    plt.tight_layout()
    plt.savefig(module_dir / "muspan_proximity_point.png")
    logger.info("Proximity networks (point-like objects) plotted and saved")


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.ms_spatial_graph")

    # Set seed
    seed_everything(21122023)

    try:
        run_muspan_graph(
            MuspanSpatialGraphModuleConfig(
                module_name="6_muspan",
                muspan_object="muspan_object.muspan",
                min_edge_distance=0,
                max_edge_distance=45,
                distance_list=(10, 20, 50),
                min_edge_distance_shape=0,
                max_edge_distance_shape=1,
                k_list=(2, 5, 10, 15),
            ),
            IOConfig(),
        )
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
    except ModuleNotFoundError as err:
        logger.error(f"Module not found: {err}")
