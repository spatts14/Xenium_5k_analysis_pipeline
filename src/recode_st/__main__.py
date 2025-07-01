"""The main entry point for the recode_st package."""

from logging import getLogger

from recode_st.annotate import run_annotate
from recode_st.dimension_reduction import run_dimension_reduction
from recode_st.format import run_format
from recode_st.logging_config import configure_logging
from recode_st.ms_spatial_graph import run_muspan_graph
from recode_st.ms_spatial_stat import run_muspan_stats
from recode_st.muspan import run_muspan
from recode_st.qc import run_qc
from recode_st.spatial_statistics import run_spatial_statistics
from recode_st.view_images import run_view_images

logger = getLogger("recode_st")


def main():
    """Main function to run the recode_st package."""
    configure_logging()

    logger.info("Starting recode_st pipeline...")

    logger.info("Running Module 0 - Format")
    run_format()

    logger.info("Running Module 1 - Quality Control")
    run_qc()

    logger.info("Running Module 2 - Dimension Reduction")
    run_dimension_reduction()

    logger.info("Running Module 3 - Annotate")
    run_annotate()

    logger.info("Running Module 4 - View Images")
    run_view_images()

    logger.info("Running Module 5 - Spatial Statistics")
    run_spatial_statistics()

    logger.info("Running Module 6 - MuSpAn")
    run_muspan()

    logger.info("Running Module 7 - MuSpAn Spatial Graph")
    run_muspan_graph()

    logger.info("Running Module 8 - MuSpAn Spatial Stats")
    run_muspan_stats()


if __name__ == "__main__":
    main()
