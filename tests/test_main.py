"""Tests for the __main__ module."""


def test_main(mocker, caplog):
    """Test that the main function runs without errors."""
    from recode_st.__main__ import main
    from recode_st.config import (
        Config,
        FormatDataModuleConfig,
        SpatialStatisticsModuleConfig,
    )

    configure_logging_mock = mocker.patch("recode_st.logging_config.configure_logging")
    seed_everything_mock = mocker.patch("recode_st.helper_function.seed_everything")

    # Test with no modules enabled in config
    config = Config(seed=42)
    with caplog.at_level("INFO"):
        main(config)

    configure_logging_mock.assert_called_once_with(
        config.io.logging_dir, config.log_level
    )
    seed_everything_mock.assert_called_once_with(config.seed)
    assert caplog.record_tuples == [
        ("recode_st", 20, "Seeding everything..."),
        ("recode_st", 20, "Starting recode_st pipeline..."),
    ]

    # Test with format_data module enabled
    convert_xenium_to_zarr_mock = mocker.patch(
        "recode_st.format_data.convert_xenium_to_zarr"
    )
    run_qc_mock = mocker.patch("recode_st.qc.run_qc")
    config.modules.format_data = FormatDataModuleConfig(module_name="0_format")
    with caplog.at_level("INFO"):
        main(config)

    convert_xenium_to_zarr_mock.assert_called_once_with(
        config.io.xenium_dir, config.io.zarr_dir
    )
    run_qc_mock.assert_not_called()
    assert caplog.record_tuples[-2:] == [
        ("recode_st", 20, "Running Module 0 - Format"),
        ("recode_st.format_data", 20, "Finished formatting data."),
    ]

    # Test with format_data and spatial_statistics modules enabled
    convert_xenium_to_zarr_mock.reset_mock()
    run_spatial_statistics_mock = mocker.patch(
        "recode_st.spatial_statistics.run_spatial_statistics"
    )
    config.modules.spatial_statistics = SpatialStatisticsModuleConfig(
        module_name="5_spatial_stats"
    )

    with caplog.at_level("INFO"):
        main(config)

    convert_xenium_to_zarr_mock.assert_called_once_with(
        config.io.xenium_dir, config.io.zarr_dir
    )
    run_spatial_statistics_mock.assert_called_once_with(
        config.modules.spatial_statistics, config.io
    )
    run_qc_mock.assert_not_called()
    assert caplog.record_tuples[-3:] == [
        ("recode_st", 20, "Running Module 0 - Format"),
        ("recode_st.format_data", 20, "Finished formatting data."),
        ("recode_st", 20, "Running Module 5 - Spatial Statistics"),
    ]
