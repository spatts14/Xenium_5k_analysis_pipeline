"""Test for the logging_config module."""

import logging
from pathlib import Path


def test_configure_logging(tmp_path, caplog):
    """Test that logging is configured correctly."""
    from recode_st.logging_config import configure_logging

    # Create a temporary directory for logging
    logging_dir = tmp_path / "logs"
    logging_dir.mkdir(parents=True, exist_ok=True)

    # Configure logging
    configure_logging(logging_dir, "DEBUG")

    # Check if the logger is set up correctly
    logger = logging.getLogger("recode_st")
    logger.debug("Test debug message.")
    log_file = Path(logger.handlers[-1].baseFilename)

    assert logger.name == "recode_st"
    assert logger.level == logging.DEBUG
    assert len(logger.handlers) == 2
    assert log_file == next(logging_dir.glob("*"))
    assert caplog.record_tuples[-1] == (
        "recode_st",
        logging.DEBUG,
        "Test debug message.",
    )
    assert log_file.read_text().endswith(" - recode_st - DEBUG] Test debug message.\n")
