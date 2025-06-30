"""The configuration module for recode_st."""

from pathlib import Path
from typing import Literal

import tomllib
from pydantic import BaseModel


class Config(BaseModel):
    """The possible configuration options for recode_st."""

    log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "INFO"
    """The logging level to use for the package."""


def load_config(config_file: str | Path) -> Config:
    """Load the configuration from a TOML file.

    Args:
        config_file: The path to the configuration TOML file.

    Returns:
        An instance of Config with the loaded settings.
    """
    config_path = Path(config_file)
    if not config_path.exists():
        raise ValueError(f"Config file {config_path} not found.")

    with config_path.open("rb") as f:
        data = tomllib.load(f)

    return Config.model_validate(data)
