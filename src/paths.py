"""Module for defining inputs and outputs paths."""

from pathlib import Path

# Set directories
base_dir = Path(  # TODO: Change to a dynamically set base directory
    "/Users/sarapatti/Desktop/PhD_projects/Llyod_lab/ReCoDe-spatial-transcriptomics"
)
data_path = base_dir / "data"
output_path = base_dir / "analysis"
xenium_path = data_path / "xenium"
zarr_path = data_path / "xenium.zarr"
logging_path = output_path / "logging"
