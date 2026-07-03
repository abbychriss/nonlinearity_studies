"""
Nonlinearity Studies Package
=============================

A Python package for analyzing nonlinearity in detector charge measurements
using FITS data and statistical analysis.

Analysis library, organized by responsibility:
- fit_data         : fit models & primitives (parabola, Gaussians, peak finding)
- fit_range        : automatic estimation of the nonlinearity fit range
- nonlinearity     : nonlinearity computation per extension / at a charge
- resolution       : single-electron resolution fitting and classification
- plotting         : peak, nonlinearity, and resolution figures
- summary          : per-extension and wide summary tables (CSV)

Other modules:
- run_nonlinearity_studies : command-line interface for running analyses
- cli_config               : config loading & argument-normalization helpers
- nonlinearity_studies     : backwards-compatibility shim re-exporting the above

The zero/one-electron peak fitting, pedestal subtraction, and FITS I/O come from
the ``pedestal_subtract`` package; FITS stitching from ``analysis_tools``.
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

# Import stitch_fits utility (now lives in the shared analysis_tools package)
from analysis_tools import stitch_fits

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
    "pedestal_subtract_ext_cached",
    "estimate_optimal_fit_range_right_ext",
    "estimate_fit_range_right_by_noise_onset_ext",
    "estimate_fit_range_right_changepoint_ext",
    "stitch_fits",
]
