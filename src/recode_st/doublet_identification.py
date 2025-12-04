"""Identify doublets and spatial overlap with ovrly."""

import warnings
from logging import getLogger

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logger = getLogger(__name__)


def run_doublet_id(adata, module_dir):
    """Identify doublets using ovrl.

    Args:
        adata (_type_): _description_
        module_dir (_type_): _description_

    Returns:
        _type_: _description_
    """
    return None
