"""Tests for the config module."""

from pathlib import Path

import pytest
from pydantic import ValidationError


def test_config():
    """Test that the configuration can be created and has default values."""
    from recode_st.config import Config

    with pytest.raises(ValidationError):
        Config()

    config = Config(seed=42)

    assert config.log_level == "INFO"
    assert config.seed == 42


def test_io_config_resolve_paths(tmp_path):
    """Test that the IOConfig can be created and has default values."""
    from recode_st.config import IOConfig

    with pytest.raises(ValidationError):
        IOConfig(base_dir="non_existent_directory")

    # Check that the default paths are set correctly
    io_config = IOConfig()

    assert io_config.base_dir == Path()
    assert io_config.data_dir == Path() / "data"
    assert io_config.output_dir == Path() / "analysis"
    assert io_config.xenium_dir == Path() / "data" / "xenium"
    assert io_config.zarr_dir == Path() / "data" / "xenium.zarr"
    assert io_config.area_path == Path() / "data" / "selected_cells_stats.csv"
    assert io_config.logging_dir == Path() / "analysis" / "logs"

    # Check that paths outside of the base directory are resolved correctly
    base_dir = tmp_path / "base"
    base_dir.mkdir()
    io_config = IOConfig(base_dir=base_dir, data_dir=tmp_path / "data")

    assert io_config.base_dir == base_dir
    assert io_config.data_dir == tmp_path / "data"
    assert io_config.output_dir == base_dir / "analysis"
    assert io_config.xenium_dir == tmp_path / "data" / "xenium"
    assert io_config.zarr_dir == tmp_path / "data" / "xenium.zarr"
    assert io_config.area_path == tmp_path / "data" / "selected_cells_stats.csv"
    assert io_config.logging_dir == base_dir / "analysis" / "logs"


def test_load_config(tmp_path):
    """Test loading configuration from a file."""
    from recode_st.config import load_config

    # Check that an error is raised for a non-existent config file
    with pytest.raises(ValueError):
        load_config("non_existent_config.toml")

    # Check that an error is raised for an empty config file
    config_file = tmp_path / "config.toml"
    config_file.write_text("", encoding="utf-8")

    with pytest.raises(ValidationError):
        load_config(config_file)

    # Check that a valid config file can be loaded
    config_file.write_text(
        f"seed=42\nlog_level='DEBUG'\n[io]\nbase_dir='{tmp_path}'", encoding="utf-8"
    )
    config = load_config(config_file)

    assert config.seed == 42
    assert config.log_level == "DEBUG"
    assert config.io.base_dir == tmp_path
