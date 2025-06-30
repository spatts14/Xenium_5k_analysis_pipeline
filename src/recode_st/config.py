"""The configuration module for recode_st."""

from pathlib import Path
from typing import Literal

import tomllib
from pydantic import BaseModel


class BaseModuleConfig(BaseModel):
    """Configuration for a specific module in recode_st."""

    module_name: str
    """The name of the module."""


class QualityControlModuleConfig(BaseModuleConfig):
    """Configuration for the Quality Control module."""

    min_counts: int
    """Minimum number of counts required for a cell to pass filtering."""

    min_cells: int
    """Minimum number of cells expressed required for a gene to pass filtering"""


class DimensionReductionModuleConfig(BaseModuleConfig):
    """Configuration for the Quality Control module."""


class AnnotateModuleConfig(BaseModuleConfig):
    """Configuration for the Annotate module."""


class ViewImagesModuleConfig(BaseModuleConfig):
    """Configuration for the View Images module."""

    gene_list: list[str]
    """List of genes to visualize on tissue."""


class SpatialStatisticsModuleConfig(BaseModuleConfig):
    """Configuration for the Spatial Statistics module."""


class MuspanModuleConfig(BaseModuleConfig):
    """Configuration for the Muspan module."""


class Config(BaseModel):
    """The possible configuration options for recode_st."""

    log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "INFO"
    """The logging level to use for the package."""

    seed: int
    """A random seed to use for reproducibility."""

    quality_control: QualityControlModuleConfig
    """Configuration for the Quality Control module."""

    dimension_reduction: DimensionReductionModuleConfig
    """Configuration for the Dimension Reduction module."""

    annotate: AnnotateModuleConfig
    """Configuration for the Annotate module."""

    view_images: ViewImagesModuleConfig
    """Configuration for the View Images module."""

    spatial_statistics: SpatialStatisticsModuleConfig
    """Configuration for the Spatial Statistics module."""

    muspan: MuspanModuleConfig
    """Configuration for the Muspan module."""


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
