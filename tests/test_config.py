"""Tests for the config module."""

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
