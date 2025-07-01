"""The configuration module for recode_st."""

from pathlib import Path
from typing import Literal, Self

import tomllib
from pydantic import BaseModel, DirectoryPath, model_validator


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

    gene_list: tuple[str, ...]
    """List of genes to visualize on tissue. A tuple of any length."""


class SpatialStatisticsModuleConfig(BaseModuleConfig):
    """Configuration for the Spatial Statistics module."""


class MuspanModuleConfig(BaseModuleConfig):
    """Configuration for the Muspan module."""


class FormatDataModuleConfig(BaseModuleConfig):
    """Configuration for the Format Data module."""


class ModulesConfig(BaseModel):
    """Configuration for all modules with optional fields."""

    format_data: FormatDataModuleConfig | None = None
    """Configuration for the Format Data module."""

    quality_control: QualityControlModuleConfig | None = None
    """Configuration for the Quality Control module."""

    dimension_reduction: DimensionReductionModuleConfig | None = None
    """Configuration for the Dimension Reduction module."""

    annotate: AnnotateModuleConfig | None = None
    """Configuration for the Annotate module."""

    view_images: ViewImagesModuleConfig | None = None
    """Configuration for the View Images module."""

    spatial_statistics: SpatialStatisticsModuleConfig | None = None
    """Configuration for the Spatial Statistics module."""

    muspan: MuspanModuleConfig | None = None
    """Configuration for the Muspan module."""


class Config(BaseModel):
    """The possible configuration options for recode_st."""

    log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "INFO"
    """The logging level to use for the package."""

    seed: int
    """A random seed to use for reproducibility."""

    base_dir: DirectoryPath = Path()
    """The base directory for all the input and output files."""

    data_dir: Path = Path("data")
    """The data directory for all the input files."""

    output_dir: Path = Path("analysis")
    """The output directory for all the output files."""

    xenium_dir: Path = Path("xenium")
    """The directory containing the Xenium input data."""

    zarr_dir: Path = Path("xenium.zarr")
    """The directory containing the Zarr-formatted input data."""

    logging_dir: Path = Path("logs")
    """The directory for logging output."""

    modules: ModulesConfig = ModulesConfig()
    """Configuration for all modules."""

    @model_validator(mode="after")
    def resolve_paths(self) -> Self:
        """Resolve the relative paths so they are relative to the base directory."""
        self.data_dir = self.base_dir / self.data_dir
        self.output_dir = self.base_dir / self.output_dir
        self.xenium_dir = self.data_dir / self.xenium_dir
        self.zarr_dir = self.data_dir / self.zarr_dir
        self.logging_dir = self.output_dir / self.logging_dir

        # Ensure directories exists
        if not Path(self.data_dir).exists():
            self.data_dir.mkdir(parents=True, exist_ok=True)
        if not Path(self.output_dir).exists():
            self.output_dir.mkdir(parents=True, exist_ok=True)
        if not Path(self.logging_dir).exists():
            self.logging_dir.mkdir(parents=True, exist_ok=True)

        return self


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
