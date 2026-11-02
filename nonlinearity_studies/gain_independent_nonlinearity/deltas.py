"""
Differential-gain observables: delta_g(n) and everything derived from it.

With q(n) the ADU charge of peak number n (n = 0 is the zero-electron peak,
see ``peaks.find_zero_peak_ext``) and g the standard gain:

    delta_g(n) = q(n) - q(n-1) - g          (the differential of the gain)
    I(n)       = sum_{i=1..n} delta_g(i)     = q(n) - q(0) - n*g
    S(n)       = I(n) / n

I(n) is the nonlinearity in ADU (theoretically alpha * n^2), so S(n) is linear
in n with slope alpha' -- an estimate of the nonlinearity parameter that never
divides the spectrum by a gain.
"""

import numpy as np

from .peaks import peak_number_mask, n_max_per_ext


def get_standard_gain_ext(q_ext):
    """g = q(1) - q(0), the standard gain, per extension (NaN if < 2 peaks)."""
    return [float(q[1] - q[0]) if len(q) >= 2 else np.nan for q in q_ext]


def get_gain_indep_nonlinearity_ext(q_ext, g_ext):
    """delta_g(n) = q(n) - q(n-1) - g, in ADU, for n = 1 .. len(q)-1.

    Args:
        q_ext: list (one per extension) of ADU charges q(n), n = 0, 1, 2, ...
        g_ext: list (one per extension) of standard gains g.

    Returns:
        list of arrays; element k of an extension's array is delta_g(n=k+1).
    """
    return [np.diff(np.asarray(q, dtype=float)) - g for q, g in zip(q_ext, g_ext)]


def fit_delta_g_line_ext(delta_g_ext, n_max=None):
    """Straight-line fit through the delta_g(n) scatter.

    Returns a list (one per extension) of dicts with ``slope``, ``intercept``
    (both ADU per peak / ADU) and ``pcov``. Extensions with too few points get
    NaNs rather than raising.
    """
    fits = []
    for delta_g, nm in zip(delta_g_ext, n_max_per_ext(n_max, len(delta_g_ext))):
        delta_g = np.asarray(delta_g, dtype=float)
        n = np.arange(1, delta_g.size + 1)
        mask = peak_number_mask(n, nm) & np.isfinite(delta_g)
        fits.append(_line_fit(n[mask], delta_g[mask]))
    return fits


def get_cumulative_delta_ext(delta_g_ext):
    """I(n) = cumulative sum of delta_g, and S(n) = I(n)/n.

    Returns (I_ext, S_ext), both lists of arrays aligned with ``delta_g_ext``
    (element k is n = k+1).
    """
    I_ext = []
    S_ext = []
    for delta_g in delta_g_ext:
        I = np.cumsum(np.asarray(delta_g, dtype=float))
        n = np.arange(1, I.size + 1)
        I_ext.append(I)
        S_ext.append(I / n)
    return I_ext, S_ext


def fit_alpha_from_S_ext(S_ext, n_max=None):
    """Fit S(n) = alpha' * n + const; the slope is the nonlinearity parameter.

    Returns a list (one per extension) of dicts with ``slope`` (alpha', in ADU
    per peak^2), ``intercept`` and ``pcov``.

    Note ``n_max`` can be a list (one per extension) or a scalar
    """
    fits = []
    for S, nm in zip(S_ext, n_max_per_ext(n_max, len(S_ext))):
        S = np.asarray(S, dtype=float)
        n = np.arange(1, S.size + 1)
        mask = peak_number_mask(n, nm) & np.isfinite(S)
        fits.append(_line_fit(n[mask], S[mask]))
    return fits


def _line_fit(x, y):
    """Unweighted straight-line least squares; NaNs instead of an exception."""
    if x.size < 2:
        return dict(slope=np.nan, intercept=np.nan, pcov=None)
    coeff, pcov = np.polyfit(x, y, 1, cov=True) if x.size > 2 else (np.polyfit(x, y, 1), None)
    return dict(slope=float(coeff[0]), intercept=float(coeff[1]), pcov=pcov)
