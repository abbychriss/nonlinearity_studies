"""
Nonlinearity Studies Package
=============================

A Python package for analyzing nonlinearity in detector charge measurements
using FITS data and statistical analysis.

Main Modules:
- nonlinearity_studies: Core analysis functions for nonlinearity studies
- run_nonlinearity_studies: Command-line interface for running analyses
- stitch_fits: Utility for stitching multi-extension FITS files
"""

__version__ = "0.1.0"
__author__ = "Abby Chriss"

# Import main analysis functions for convenience
from .nonlinearity_studies import (
    convert_to_electrons,
    calculate_noise_gain,
    get_fits,
    get_zero_one_peaks_ext,
    get_all_peaks_ext,
    get_nonlinearity_ext,
    get_nonlinearity_at_ext,
    resolution_at_charge,
    resolution_at_charge_ext,
    classify_resolution,
    format_resolution_table,
    summarize_resolution,
    build_resolution_summary_wide,
    build_extension_summary,
    format_extension_summary,
    summarize_extensions,
    plot_zero_one_peaks,
    plot_all_peaks,
    plot_nonlinearity,
    plot_resolution,
    pedestal_subtract_ext_cached,
    estimate_optimal_fit_range_right_ext,
    estimate_fit_range_right_by_noise_onset_ext,
    estimate_fit_range_right_changepoint_ext,
)

# Import stitch_fits utility
from .stitch_fits import stitch_fits

__all__ = [
    "convert_to_electrons",
    "calculate_noise_gain",
    "get_fits",
    "get_zero_one_peaks_ext",
    "get_all_peaks_ext",
    "get_nonlinearity_ext",
    "get_nonlinearity_at_ext",
    "resolution_at_charge",
    "resolution_at_charge_ext",
    "classify_resolution",
    "format_resolution_table",
    "summarize_resolution",
    "build_resolution_summary_wide",
    "build_extension_summary",
    "format_extension_summary",
    "summarize_extensions",
    "plot_zero_one_peaks",
    "plot_all_peaks",
    "plot_nonlinearity",
    "plot_resolution",
    "stitch_fits",
]
