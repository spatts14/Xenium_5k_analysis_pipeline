"""Module containing helper functions."""

import os
import random

import numpy as np
import scanpy as sc
import torch


def seed_everything(seed: int):
    """Set random seed on every random module for reproducibility.

    Args:
        seed: The seed value to set for random number generation.
    """
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = True
    elif torch.backends.mps.is_available():
        torch.mps.manual_seed(seed)
    else:
        pass


def set_figure_params(config, module_dir, method):
    """Set figure parameters for scanpy plots."""
    sc.set_figure_params(
        dpi=300,
        dpi_save=300,
        frameon=False,
        vector_friendly=True,
        fontsize=16,
        facecolor="white",
        figsize=(15, 4),
    )
    sc.settings.figdir = module_dir
