# Configuration Guide

ReCoDe uses TOML configuration files to control which analysis modules to run and their parameters. Only modules defined in the configuration file will be executed.

## Basic Configuration Structure

```toml
# Global settings
log_level = "INFO"  # DEBUG, INFO, WARNING, ERROR, CRITICAL (default: INFO)
seed = 21122023     # Random seed for reproducibility

# IO config settings
[io]
base_dir = "."      # Base directory for input/output files (default: Current Working Directory)

# Module configurations (only include modules you want to run)
[modules.module_name]
module_name = "output_folder_name"
# module-specific parameters...
```

## IO Config

The IO config allows you to specify where to get data from and where to save it. The options
and their defaults are:

```toml
[io]
base_dir = "."
data_dir = "data"        # relative to base_dir
output_dir = "analysis"  # relative to base_dir
xenium_dir = "xenium"    # relative to data_dir
zarr_dir = "xenium.zarr" # relative to data_dir
logging_path = "logs"    # relative to output_dir
```

This resolves to the following default directory structure:

```bash
base_dir
├── analysis
│   ├── output_folder_name/
│   ├── ...
│   └── logs/
├── config.toml
└── data
    ├── xenium/
    └── xenium.zarr/
```

## Available Modules

- `format_data` - Data formatting and preprocessing
- `quality_control` - Cell and gene filtering (requires `min_counts`, `min_cells`)
- `dimension_reduction` - PCA, UMAP, etc.
- `annotate` - Cell type annotation
- `view_images` - Tissue visualization (requires `gene_list`)
- `spatial_statistics` - Spatial analysis
- `muspan` - Advanced spatial analysis (requires MuSpAn license)

## Example 1: Basic Quality Control Pipeline

```toml
log_level = "INFO"
seed = 12345
base_dir = "."

[modules.format_data]
module_name = "0_format"

[modules.quality_control]
module_name = "1_qc"
min_counts = 10
min_cells = 5

[modules.dimension_reduction]
module_name = "2_dr"
```

## Example 2: Full Analysis with Visualization and specifying the IO

```toml
log_level = "INFO"
seed = 12345
base_dir = "."
data_dir = "data"
output_dir = "analysis"
xenium_dir = "xenium"
zarr_dir = "xenium.zarr"
logging_path = "logs"

[modules.format_data]
module_name = "0_format"

[modules.quality_control]
module_name = "1_qc"
min_counts = 15
min_cells = 3

[modules.dimension_reduction]
module_name = "2_dr"

[modules.annotate]
module_name = "3_annotate"

[modules.view_images]
module_name = "4_images"
gene_list = ["EPCAM", "CD3D", "CD68", "PTPRC"]

[modules.spatial_statistics]
module_name = "5_spatial"
```

## Running with Configuration

```bash
# Run with a configuration file
python -m recode_st config.toml
```

The pipeline will only run the modules you specify in your configuration file, allowing you to customize your analysis workflow.
