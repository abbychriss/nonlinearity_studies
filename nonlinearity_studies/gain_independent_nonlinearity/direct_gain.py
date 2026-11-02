"""
The gain read straight off the raw spectrum, with no fit at all.

Writing the charge of the n-th peak as q(n) = n * g_dir(n) defines the direct
gain g_dir(n) = (q(n) - q(0)) / n -- the average ADU per electron out to peak n,
using nothing but the peak number and its measured charge. Its deviation from
the standard gain,

    delta_dir(n) = g - g_dir(n),

is the direct counterpart of delta_g. For a CCD whose standard nonlinearity is
negative (measured charge in electrons falls below the peak number), q/g < n
while q/g_dir(n) = n by construction, so g_dir(n) < g and delta_dir(n) > 0 for
all n > 1.

Also computes the cumulative sum of delta_dir(n), exactly as it's done for differential of the gain delta_g,
and 

    I_dir(n)       = sum_{i=1..n} delta_dir(i) 
    S_dir(n)       = I_dir(n) / n
"""

import numpy as np

from .deltas import _line_fit, get_cumulative_delta_ext


def get_direct_gain_ext(q_ext):
    """g_dir(n) = (q(n) - q(0)) / n for n = 1 .. len(q)-1, in ADU per electron.

    Returns a list of arrays; element k is g_dir(n=k+1).
    """
    g_dir_ext = []
    for q in q_ext:
        q = np.asarray(q, dtype=float)
        n = np.arange(1, q.size)
        g_dir_ext.append((q[1:] - q[0]) / n if q.size >= 2 else np.zeros(0))

    return g_dir_ext


def get_direct_delta_ext(g_dir_ext, g_ext):
    """delta_dir(n) = g - g_dir(n), in ADU per electron, aligned with g_dir."""

    delta_dir_ext = [g - np.asarray(g_dir, dtype=float) for g_dir, g in zip(g_dir_ext, g_ext)]

    I_dir_ext, S_dir_ext = get_cumulative_delta_ext(delta_dir_ext)

    return delta_dir_ext, I_dir_ext, S_dir_ext


def fit_direct_delta_line_ext(delta_dir_ext):
    """Straight-line fit through the delta_dir(n) scatter, one per extension.

    Same dicts (``slope``, ``intercept``, ``pcov``) as ``deltas.fit_delta_g_line_ext``.
    """
    fits = []
    for delta_dir in delta_dir_ext:
        delta_dir = np.asarray(delta_dir, dtype=float)
        n = np.arange(1, delta_dir.size + 1)
        mask = np.isfinite(delta_dir)
        fits.append(_line_fit(n[mask], delta_dir[mask]))
    return fits
