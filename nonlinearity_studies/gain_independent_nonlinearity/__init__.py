"""
Gain-independent nonlinearity analysis.

The production analysis converts the charge spectrum to electrons with a gain
measured from the zero and one electron peaks alone, then reports how far each
peak's charge falls from its peak number. This package leaves the spectrum in
ADU and derives the nonlinearity without ever dividing by that gain, so the
measurement does not assume a constant gain a priori.

Modules, one per topic:
- peaks        : sub-bin peak localization and peak numbering (anchoring n=0 on
                 the zero-electron peak)
- deltas       : the differential gain delta_g(n), its cumulative sum I(n),
                 S(n) = I(n)/n, and the line fits through them
- spectrum_fit : polynomial fits to q(n) -- a0 + g_fit*n + alpha*n^2 -- the
                 running fit g_fit(i), and the gain implied by the fit
- direct_gain  : g_dir(n) = (q(n) - q(0))/n and its deviation from g
- gain_scan    : recomputing the standard nonlinearity with each candidate gain
- plots        : one figure per calculation
- summary      : per-extension summary table (CSV)
- run          : command-line interface tying it all together
"""

from .peaks import (
    sub_bin_log_parabola,
    gaussian_refit_failed_peaks,
    find_zero_peak_ext,
    peak_number_mask,
    n_max_per_ext,
)
from .deltas import (
    get_standard_gain_ext,
    get_gain_indep_nonlinearity_ext,
    fit_delta_g_line_ext,
    get_cumulative_delta_ext,
    fit_alpha_from_S_ext,
)
from .spectrum_fit import (
    estimate_fit_n_max_ext,
    fit_charge_spectrum_ext,
    fit_charge_spectrum_cubic_ext,
    get_cubic_truncation_ext,
    get_gain_from_fit_ext,
    fit_charge_spectrum_running_ext,
    get_gain_fit_residual_ext,
)
from .direct_gain import get_direct_gain_ext, get_direct_delta_ext
from .gain_scan import scan_gain_nonlinearity_ext
from .summary import (
    build_gain_indep_summary,
    format_gain_indep_summary,
    write_gain_indep_summary,
)
from . import plots


def __getattr__(name):
    """Expose ``main`` / ``init_argparse`` without importing the CLI at package
    import time -- importing it here would make ``python -m
    nonlinearity_studies.gain_independent_nonlinearity.run`` load run.py twice
    and warn about it."""
    if name in ('main', 'init_argparse'):
        from . import run
        return getattr(run, name)
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')

__all__ = [
    "sub_bin_log_parabola",
    "gaussian_refit_failed_peaks",
    "find_zero_peak_ext",
    "peak_number_mask",
    "n_max_per_ext",
    "get_standard_gain_ext",
    "get_gain_indep_nonlinearity_ext",
    "fit_delta_g_line_ext",
    "get_cumulative_delta_ext",
    "fit_alpha_from_S_ext",
    "estimate_fit_n_max_ext",
    "fit_charge_spectrum_ext",
    "fit_charge_spectrum_cubic_ext",
    "get_cubic_truncation_ext",
    "get_gain_from_fit_ext",
    "fit_charge_spectrum_running_ext",
    "get_gain_fit_residual_ext",
    "get_direct_gain_ext",
    "get_direct_delta_ext",
    "scan_gain_nonlinearity_ext",
    "build_gain_indep_summary",
    "format_gain_indep_summary",
    "write_gain_indep_summary",
    "plots",
    "main",
    "init_argparse",
]
