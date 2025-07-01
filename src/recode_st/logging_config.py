"""Module for configuring logging in the recode_st package."""

import logging
from pathlib import Path


def configure_logging(logging_dir: Path, log_level=logging.INFO):
    """Configure logging for the package.

    Args:
        logging_dir: The directory where log files will be stored.
        log_level: The logging level to set. Defaults to logging.INFO.
    """
    logger = logging.getLogger("recode_st")
    logger.setLevel(log_level)

    ch = logging.StreamHandler()
    fh = logging.FileHandler(logging_dir / "recode_st.log", mode="w")

    formatter = logging.Formatter(
        "[%(asctime)s - %(name)s - %(levelname)s] %(message)s",
        datefmt="%Y/%m/%d %H:%M:%S",
    )

    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    logger.addHandler(ch)
    logger.addHandler(fh)
