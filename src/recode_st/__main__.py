"""The main entry point for the recode_st package."""

import sys
from logging import getLogger
from pathlib import Path

from recode_st.config import Config, load_config

logger = getLogger(__package__)


def main(config: Config):
    """Main function to run the recode_st package."""
    from recode_st.annotate import run_annotate
    from recode_st.dimension_reduction import run_dimension_reduction
    from recode_st.format_data import run_format
    from recode_st.helper_function import seed_everything
    from recode_st.logging_config import configure_logging
    from recode_st.ms_spatial_graph import run_muspan_graph
    from recode_st.ms_spatial_stat import run_muspan_stats
    from recode_st.muspan import run_muspan
    from recode_st.qc import run_qc
    from recode_st.spatial_statistics import run_spatial_statistics
    from recode_st.view_images import run_view_images

    configure_logging(config.io.logging_dir, config.log_level)

    logger.info("Seeding everything...")
    seed_everything(config.seed)

    logger.info("Starting recode_st pipeline...")

    if config.modules.format_data:
        logger.info("Running Module 0 - Format")
        run_format(config.io)

    if config.modules.quality_control:
        logger.info("Running Module 1 - Quality Control")
        run_qc(config.modules.quality_control, config.io)

    if config.modules.dimension_reduction:
        logger.info("Running Module 2 - Dimension Reduction")
        run_dimension_reduction(config.modules.dimension_reduction, config.io)

    if config.modules.annotate:
        logger.info("Running Module 3 - Annotate")
        run_annotate(config.modules.annotate, config.io)

    if config.modules.view_images:
        logger.info("Running Module 4 - View Images")
        run_view_images(config.modules.view_images, config.io)

    if config.modules.spatial_statistics:
        logger.info("Running Module 5 - Spatial Statistics")
        run_spatial_statistics(config.modules.spatial_statistics, config.io)

    if config.modules.muspan:
        logger.info("Running Module 6 - MuSpAn")
        run_muspan(config.modules.muspan, config.io)

    if config.modules.muspan_spatial_graph:
        logger.info("Running Module 7 - MuSpAn Spatial Graph")
        run_muspan_graph(config.modules.muspan_spatial_graph, config.io)

    if config.modules.muspan_spatial_stat:
        logger.info("Running Module 8 - MuSpAn Spatial Stats")
        run_muspan_stats(config.modules.muspan_spatial_stat, config.io)


if __name__ == "__main__":
    try:
        config_file = Path(sys.argv[1])
    except IndexError:
        raise Exception("No config file specified.")

    config = load_config(config_file)

    main(config)
