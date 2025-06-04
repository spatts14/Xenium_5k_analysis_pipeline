"""Tests for the helper_functions module."""

import os


def test_seed_everything():
    """Test the seed_everything function."""
    from recode_st.helper_function import seed_everything

    assert os.environ.get("PYTHONHASHSEED") is None

    seed_everything(42)
    assert os.environ.get("PYTHONHASHSEED") == "42"
