"""Identify doublets and spatial overlap with ovrly."""

import pickle
import polars as pl
import warnings
from logging import getLogger
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ovrlpy
import pandas as pd
import scanpy as sc
import seaborn as sns

from recode_st.config import DoubletIdentificationModuleConfig, IOConfig
from recode_st.helper_function import configure_scanpy_figures

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)


def extract_roi_name(folder_name: str) -> str:
    """Extract ROI name from folder name."""
    # Example: "output-...__IPF_RBH_19__20251001__141533" → "IPF_RBH_19"
    parts = folder_name.split("__")
    for part in parts:
        if part.startswith(("IPF", "COPD", "PM08")):
            return part
    raise ValueError(f"Could not parse ROI name from: {folder_name}")


def load_transcripts(io_config: IOConfig):
    """Load transcripts file for each ROI."""
    xenium_path = Path(io_config.xenium_dir)

    output_data_dir = Path(io_config.output_data_dir)
    output_data_dir.mkdir(parents=True, exist_ok=True)

    if not xenium_path.exists():
        raise FileNotFoundError(f"Xenium root directory not found: {xenium_path}")

    # Loop over date/run folders directly inside Xenium root
    for date_dir in xenium_path.iterdir():
        if not date_dir.is_dir():
            continue

        run_name = date_dir.name  # e.g. "20251001__141239__SP25164_SARA_PATTI_RUN_1"
        logger.info(f"Processing run: {run_name}")

        # Loop over ROI folders
        for roi_folder in date_dir.iterdir():
            if roi_folder.is_dir() and roi_folder.name.startswith("output-"):
                roi_name = extract_roi_name(roi_folder.name)
                logger.info(f"Processing ROI: {roi_name}")

                try:
                    df = ovrlpy.io.read_Xenium(roi_folder / "transcripts.parquet")
                    logger.info(f"Loaded Xenium data for {roi_name}")
                except Exception as err:
                    logger.error(f"Failed loading Xenium data for {roi_name}: {err}")
                    continue

        return df


def run_doublet_id(io_config: IOConfig, config: DoubletIdentificationModuleConfig):
    """Identify doublets using ovrl.

    Args:
        io_config (_type_): _description_
        config (_type_): _description_

    Returns:
        _type_: _description_
    """
    # Name of the column to store label transfer results in adata.obs
    module_dir = io_config.output_dir / config.module_name

    # Create output directories if they do not exist
    module_dir.mkdir(exist_ok=True)

    # Set figure directory for this module (overrides global setting)
    sc.settings.figdir = module_dir

    # Set figure settings to ensure consistency across all modules
    configure_scanpy_figures(str(io_config.output_dir))
    cmap = sns.color_palette("Spectral", as_cmap=True)

    logger.info("Loading transcripts form each ROI into a df...")
    # Need to calculate doublets for each ROI independently
    # Load "transcripts.parquet" for each ROI
    # Make a list of dfs and each df in the list is an df for each ROI
    # ? SHOULD I I have a function that 1. Loads the df from the folder and then 2. combined all the df?
    df_all = []
    df = load_transcripts(io_config)
    df_all = df_all.append(df)

    for roi_df in df_all:
        logger.info(f"Identifying doublets for {roi_df} using ovrly ...")
        n = 1_000  # Number of transcripts to show
        logger.info(f"Showing every {n}th transcript for {roi_df}...")

        # Plot transcripts
        fig, ax = plt.subplots()
        ax.scatter(roi_df[::n, "x"], roi_df[::n, "y"], s=0.1)
        _ = ax.set(aspect="equal")
        # TODO: Save figure

        # TODO: Subset on HVG? Cell type specific genes? NEED TO THINK

        logger.info("Running the ovrlpy pipeline...")
        roi_ob = ovrlpy.Ovrlp(
            roi_df,
            n_workers=8,  # I think this is the number of CPUs used?
            random_state=42,  # TODO: Determine if I need this since seed everything should deal with this
        )
        roi_ob.analyse()

        logger.info("Visualizing results ...")
        _ = ovrlpy.plot_pseudocells(roi_ob)
        # TODO: Save figure

        # plot SVI per ROI
        fig = ovrlpy.plot_signal_integrity(roi_ob, signal_threshold=3)
        # TODO: Save figure

        logger.info("Detecting doublets...")
        doublets = roi_ob.detect_doublets(min_signal=3, integrity_sigma=2)
        logger.info(f"Detected f{len(doublets)} in {roi_df}")

        logger.info("Visualizing doublets...")
        # Overview of tissue
        fig, ax = plt.subplots()
        _scatter = ax.scatter(
            doublets["x"], doublets["y"], c=doublets["integrity"], s=0.2, cmap=cmap
        )
        _ = ax.set_aspect("equal")
        _ = fig.colorbar(_scatter, ax=ax)
        # TODO: Save figure

        # Visualize specific doublet case
        doublet_case = np.random(
            len(doublets)
        )  # select a random doublet from length 1:len(doublets)
        x, y = doublets["x", "y"].row(doublet_case)
        _ = ovrlpy.plot_region_of_interest(roi_ob, x, y, window_size=60)
        # TODO: Save figure
        plot_names = f"{roi_df}_doublet_{doublet_case}.png"

        logger.info("Visualize cells in Z stacks using pickle...")
        # Make folder to save pickle files
        pickle_dir = module_dir / "pickle_doublets"
        pickle_dir.mkdir(exist_ok=True)

        # Name pickle file
        file = pickle_dir / f"{roi_df}.pickle"

        # If pickle file exists, load it; if not, create it
        # TODO: figure out how to save and load this correctly
        if file.exists(): #TODO: check this works
            logger.info(f"Loading {file}...")
            with open("my_analysis.pickle", "rb") as file:
            roi_ob = pickle.load(file) #TODO: make sure it loads the variable as the name
        else:
            logger.info(f"Creating {file}...")
            with open("my_analysis.pickle", "wb") as file:
                pickle.dump(roi_ob, file)  # TODO: Determine how to save to specific location

        logger.info("Plotting cells in Z stacks...")
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        plt.title(f"3D visualization of {roi_df}")
        for i in range(-2, 3):
            subset = roi_ob.transcripts.filter(
                (pl.col("z") - pl.col("z_center")).is_between(i, i + 1)
            )
            # downsample the number of transcripts
            subset = subset[::100]
            ax.scatter(subset["x"], subset["y"], i, s=1, alpha=0.1)

        ratio = roi_ob.transcripts["x"].max() / roi_ob.transcripts["y"].max()
        ax.set_box_aspect([ratio, 1, 0.75])
