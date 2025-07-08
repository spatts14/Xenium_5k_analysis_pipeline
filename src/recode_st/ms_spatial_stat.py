"""Muspan module - spatial statistics and graph analysis."""

import warnings
from logging import getLogger

import matplotlib.pyplot as plt
import numpy as np

from recode_st.config import IOConfig, MuspanSpatialStatModuleConfig
from recode_st.helper_function import seed_everything
from recode_st.logging_config import configure_logging

warnings.filterwarnings("ignore")

logger = getLogger(__name__)


def run_muspan_stats(config: MuspanSpatialStatModuleConfig, io_config: IOConfig):
    """Run Muspan spatial statistics analysis on Xenium data."""
    try:
        import muspan as ms
    except ModuleNotFoundError as err:
        logger.error(
            "Could not load necessary MuSpAn package. You can obtain this with:\n"
            "    pip install 'recode_st[muspan] @ git+"
            "https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics.git"
        )
        raise err
    # Set variables from config
    module_dir = io_config.output_dir / config.module_name
    muspan_object = config.muspan_object
    cluster_labels = config.cluster_labels

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Import data
    logger.info("Loading MuSpAn object...")
    domain = ms.io.load_domain(path_to_domain=str(module_dir / muspan_object))

    # Calculate pairwise cross-PCF
    logger.info("Calculating pairwise cross-PCF for all cell types...")

    # Define the cell types to be analyzed
    all_cluster_labels = domain.labels["cell_type"]["labels"].tolist()
    unique_cluster = np.unique(all_cluster_labels).astype(str).tolist()

    cell_type_1 = unique_cluster[0]  # TODO: how can I move this to the top?
    cell_type_2 = unique_cluster[1]  # TODO: how can I move this to the top?

    logger.info(f"Found {len(unique_cluster)} unique number of cell types")

    logger.info("Calculating cross-PCF for pairwise cell types...")

    # Create a n x n plot for visualizing the cross-PCF for combination of cell types
    fig, axes = plt.subplots(len(unique_cluster), len(unique_cluster), figsize=(40, 40))

    # Loop through each combination of cell types
    for i in range(len(unique_cluster)):
        for j in range(len(unique_cluster)):
            pop_A = ms.query.query(
                domain, ("label", cluster_labels), "is", unique_cluster[i]
            )
            pop_B = ms.query.query(
                domain, ("label", cluster_labels), "is", unique_cluster[j]
            )
            r, PCF = ms.spatial_statistics.cross_pair_correlation_function(
                domain,
                pop_A,
                pop_B,
                max_R=200,
                annulus_step=5,
                annulus_width=25,
            )

            # Select the current subplot
            ax = axes[i, j]

            # Plot the cross-PCF
            ax.plot(r, PCF)

            # Add a horizontal line at y=1 to indicate the CSR baseline
            ax.axhline(1, color="k", linestyle=":")

            # Set the y-axis limit
            ax.set_ylim([0, 7])

            # Label the y-axis with the cross-PCF notation
            ax.set_ylabel(f"$g_{{{unique_cluster[i]}{unique_cluster[j]}}}(r)$")

            # Label the x-axis with the distance r
            ax.set_xlabel("$r$")

    # Adjust the layout to prevent overlap
    plt.tight_layout()
    plt.savefig(module_dir / "cross_pair_correlation_function_all.png")
    logger.info("Cross-PCF pairwise for all cell types plot saved.")

    # Cross pair correlation function (cross-PCF)
    logger.info(f"Calculating cross-PCF for {cell_type_1} and {cell_type_2}...")
    pop_A = ms.query.query(domain, ("label", cluster_labels), "is", str(cell_type_1))
    pop_B = ms.query.query(domain, ("label", cluster_labels), "is", str(cell_type_2))

    # Calculate the cross-PCF for points of Celltype D with themselves
    # max_R: maximum radius to consider
    # annulus_step: step size for the annulus
    # annulus_width: width of the annulus
    # visualise_output: whether to visualise the output
    r, PCF = ms.spatial_statistics.cross_pair_correlation_function(
        domain=domain,
        population_A=pop_A,
        population_B=pop_B,
        max_R=200,
        annulus_step=5,
        annulus_width=25,
        visualise_output=True,
    )

    # Save the cross-PCF plot
    plt.savefig(
        module_dir / f"cross_pair_correlation_function_{cell_type_1}_{cell_type_2}.png"
    )
    logger.info(f"Cross-PCF for {cell_type_1} and {cell_type_2} plotted and saved'")

    # Query points with labels cell_type_1 or cell_type_2
    query_1_2 = ms.query.query(
        domain,
        ("label", cluster_labels),
        "in",
        [str(cell_type_1), str(cell_type_2)],
    )

    ms.visualise.visualise(
        domain, color_by=("label", cluster_labels), objects_to_plot=query_1_2
    )
    # Plot the cross-PCF
    plt.savefig(module_dir / f"visualize_{cell_type_1}_{cell_type_2}.png")
    logger.info(f"{cell_type_1} and {cell_type_2} plotted and saved'")


if __name__ == "__main__":
    # Set up logger
    configure_logging()
    logger = getLogger("recode_st.ms_spatial_stat")

    # Set seed
    seed_everything(21122023)

    try:
        run_muspan_stats(
            MuspanSpatialStatModuleConfig(
                module_name="6_muspan",
                muspan_object="muspan_object.muspan",
                cluster_labels="cell_type",
            ),
            IOConfig(),
        )
    except FileNotFoundError as err:
        logger.error(f"File not found: {err}")
