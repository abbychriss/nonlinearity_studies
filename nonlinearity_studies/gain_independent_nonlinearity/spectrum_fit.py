"""
Polynomial fits to the ADU charge spectrum q(n).

    Q(n) = a0 + g_fit * n + alpha * n^2      (+ beta * n^3 for the cubic)

a0 is any residual pedestal left over by the pedestal subtraction, g_fit is the
average gain over the fit domain (the linear component), and alpha is the
nonlinearity. Fitting up to peak number i instead of the whole spectrum turns
the linear term into g_fit(i), the average gain over the first i peaks.
"""

import numpy as np

from ..fit_range import _changepoint_break_index
from .peaks import peak_number_mask, n_max_per_ext


def fit_charge_spectrum_ext(q_ext, n_max=None, deg=2):
    """Fit Q(n) to a polynomial of degree ``deg`` over peak numbers 0 .. n_max-1.

    ``n_max`` is a scalar applied to every extension, or one value per
    extension (what ``estimate_fit_n_max_ext`` returns); None fits every peak.

    Returns a list (one per extension) of dicts:
        coeffs   polynomial coefficients in ASCENDING order (a0, g_fit, alpha, ...)
        a0, g_fit, alpha, beta   named coefficients (beta is NaN for deg < 3)
        pcov     covariance matrix (ascending order to match ``coeffs``), or None
        n_used   number of peaks the fit actually used
    """
    fits = []
    for q, nm in zip(q_ext, n_max_per_ext(n_max, len(q_ext))):
        q = np.asarray(q, dtype=float)
        n = np.arange(q.size)
        mask = peak_number_mask(n, nm) & np.isfinite(q)
        fits.append(_poly_fit(n[mask], q[mask], deg))
    return fits


def fit_charge_spectrum_cubic_ext(q_ext, n_max=None):
    """Cubic version of ``fit_charge_spectrum_ext``; the beta * n^3 term
    quantifies what truncating the series at n^2 leaves out."""
    return fit_charge_spectrum_ext(q_ext, n_max=n_max, deg=3)


def estimate_fit_n_max_ext(q_ext, g_ext, win=25, factor=4.0, k=4.0, persist=4,
                           floor_e=0.15, verbose=False):
    """Per-extension fit domain, found by the production changepoint detector.

    Feeds x = n, y = q(n) -- the raw ADU charge spectrum -- to
    ``fit_range._changepoint_break_index``, the same two-stage local-roughness
    detector the production pipeline uses to place fit_range_right. Both of its
    stages difference out any linear and quadratic component of y, so the raw
    spectrum goes in as-is: no need to subtract the fitted line a0 + g_fit*n
    (which would be circular, that fit being the thing n_max defines) and no
    division by a gain, so the estimate stays gain independent.

    ``floor_e`` is the detector's absolute residual floor in electrons -- the
    unit the production defaults are tuned in -- and is scaled by g to ADU per
    extension.

    Returns a list of n_max values, counting peaks *including* the
    zero-electron peak (the ``peak_number_mask`` convention). Extensions with
    too few peaks to run the detector get None, i.e. every peak is fit.
    """
    out = []
    for ext, (q, g) in enumerate(zip(q_ext, g_ext)):
        q = np.asarray(q, dtype=float)
        n = np.arange(q.size)
        good = np.isfinite(q)
        floor = floor_e * g if np.isfinite(g) else floor_e
        try:
            break_idx, _ = _changepoint_break_index(
                n[good], q[good], win=win, factor=factor, k=k,
                persist=persist, floor=floor)
        except ValueError as exc:
            out.append(None)
            if verbose:
                print(f'EXT {ext + 1}: changepoint fit range unavailable ({exc}); '
                      f'fitting every peak')
            continue
        n_max = int(n[good][break_idx]) + 1
        out.append(n_max)
        if verbose:
            print(f'EXT {ext + 1}: changepoint fit range n_max = {n_max} '
                  f'of {q.size} peaks (floor = {floor:.3f} ADU)')
    return out


def get_cubic_truncation_ext(quad_fits_ext, cubic_fits_ext, q_ext, n_max=None):
    """How much the quadratic and cubic fits disagree across the fit domain.

    Returns a list (one per extension) of dicts:
        beta          the cubic coefficient
        max_diff_adu  max |Q_cubic(n) - Q_quad(n)| over the fit domain, in ADU
        diff_at_nmax  Q_cubic - Q_quad at the last fitted peak, in ADU
        alpha_quad, alpha_cubic   the n^2 coefficient from each fit
    """
    out = []
    for quad, cubic, q, nm in zip(quad_fits_ext, cubic_fits_ext, q_ext,
                                  n_max_per_ext(n_max, len(q_ext))):
        n = np.arange(np.asarray(q).size)
        n = n[peak_number_mask(n, nm)]
        if n.size == 0 or not np.isfinite(quad['g_fit']) or not np.isfinite(cubic['g_fit']):
            out.append(dict(beta=np.nan, max_diff_adu=np.nan, diff_at_nmax=np.nan,
                            alpha_quad=quad['alpha'], alpha_cubic=cubic['alpha']))
            continue
        diff = _polyval_asc(cubic['coeffs'], n) - _polyval_asc(quad['coeffs'], n)
        out.append(dict(beta=cubic['beta'],
                        max_diff_adu=float(np.max(np.abs(diff))),
                        diff_at_nmax=float(diff[-1]),
                        alpha_quad=quad['alpha'], alpha_cubic=cubic['alpha']))
    return out


def get_gain_from_fit_ext(fits_ext, g_ext):
    """The gain implied by the parabola fit, and how it compares to the standard one.

    G = Q(1) - Q(0) = g_fit + alpha is the fit's own zero-to-one-electron step,
    i.e. the fitted counterpart of the standard gain g.

    Returns a list (one per extension) of dicts with ``G``, ``g``, ``G_minus_g``
    and ``g_fit_minus_g``, all in ADU.
    """
    out = []
    for fit, g in zip(fits_ext, g_ext):
        G = fit['g_fit'] + fit['alpha']
        out.append(dict(G=float(G), g=float(g), G_minus_g=float(G - g),
                        g_fit_minus_g=float(fit['g_fit'] - g)))
    return out


def fit_charge_spectrum_running_ext(q_ext, n_min=2, n_max=None, deg=2):
    """Refit Q(n) using only the first i peaks, for every i from n_min upwards.

    i counts peaks including the zero-electron peak, so i=2 uses peaks 0 and 1
    -- two points, where the parabola is degenerate and a straight line is fit
    instead, reproducing the standard gain exactly (g_fit(2) = q(1) - q(0) = g).
    That makes the scan cover the whole space, standard gain included.

    Returns a list (one per extension) of dicts of arrays: ``i``, ``a0``,
    ``g_fit``, ``alpha`` (alpha is 0 for the two-point line fits).
    """
    out = []
    for q, nm in zip(q_ext, n_max_per_ext(n_max, len(q_ext))):
        q = np.asarray(q, dtype=float)
        n_total = q.size if nm is None else min(q.size, int(nm))
        i_values = np.arange(max(2, int(n_min)), n_total + 1)
        a0 = np.full(i_values.size, np.nan)
        g_fit = np.full(i_values.size, np.nan)
        alpha = np.full(i_values.size, np.nan)

        n_all = np.arange(q.size)
        for k, i in enumerate(i_values):
            fit = _poly_fit(n_all[:i], q[:i], min(deg, i - 1), cov=False)
            a0[k], g_fit[k], alpha[k] = fit['a0'], fit['g_fit'], fit['alpha']

        # A line fit has no n^2 term; report it as exactly zero rather than NaN
        # so the running alpha(i) curve starts at the standard-gain point.
        alpha = np.where(np.isnan(alpha) & np.isfinite(g_fit), 0.0, alpha)
        out.append(dict(i=i_values, a0=a0, g_fit=g_fit, alpha=alpha))
    return out


def fit_running_gain_window_ext(running_ext, n_min=None, n_max=None):
    """Summarize g_fit(i) over a fixed window of fit ranges n_min <= i <= n_max.

    The running curve itself is unchanged -- this only characterizes the stretch
    of it the window selects: its average, and the least-squares line through it.

    Returns a list (one per extension) of dicts:
        i                  the window's i values (empty if the window misses the data)
        mean               the average of g_fit(i) over the window
        slope, intercept   the line g_fit(i) = slope * i + intercept
    """
    out = []
    for run in running_ext:
        i = np.asarray(run['i'], dtype=float)
        g_fit = np.asarray(run['g_fit'], dtype=float)
        mask = np.isfinite(g_fit)
        if n_min is not None:
            mask &= i >= n_min
        if n_max is not None:
            mask &= i <= n_max
        line = _poly_fit(i[mask], g_fit[mask], 1, cov=False)
        out.append(dict(i=i[mask],
                        mean=float(np.mean(g_fit[mask])) if mask.any() else np.nan,
                        slope=line['g_fit'], intercept=line['a0']))
    return out


def get_gain_fit_residual_ext(running_ext, g_ext):
    """How the fitted average gain drifts away from the standard gain.

    residual(i) = g_fit(i) - g is the fitted-spectrum counterpart of delta_g;
    its cumulative sum over i is another estimate of the nonlinearity in ADU.

    Returns (residual_ext, cumsum_ext), lists of arrays aligned with
    ``running_ext[ext]['i']``.
    """
    residual_ext = [np.asarray(run['g_fit'], dtype=float) - g
                    for run, g in zip(running_ext, g_ext)]
    return residual_ext, [np.cumsum(np.nan_to_num(r)) for r in residual_ext]


def fit_gain_fit_residual_cumsum_ext(running_ext, cumsum_ext, deg=2):
    """Parabola through the integrated residual sum_i (g_fit(i) - g) vs i.

    Returns the same dicts ``_poly_fit`` returns, one per extension, so the
    ascending ``coeffs`` (a0, linear, quadratic) evaluate straight onto i.
    """
    fits = []
    for run, cumsum in zip(running_ext, cumsum_ext):
        i = np.asarray(run['i'], dtype=float)
        cumsum = np.asarray(cumsum, dtype=float)
        mask = np.isfinite(cumsum)
        fits.append(_poly_fit(i[mask], cumsum[mask], deg, cov=False))
    return fits


def _poly_fit(n, q, deg, cov=True):
    """np.polyfit wrapper returning ASCENDING coefficients and named terms."""
    named = dict(coeffs=np.full(deg + 1, np.nan), pcov=None,
                 a0=np.nan, g_fit=np.nan, alpha=np.nan, beta=np.nan, n_used=int(n.size))
    if n.size < deg + 1 or deg < 1:
        return named

    # np.polyfit's Bayesian covariance needs more points than order + 2.
    want_cov = cov and n.size > deg + 2
    if want_cov:
        coeff, pcov = np.polyfit(n, q, deg, cov=True)
        named['pcov'] = pcov[::-1, ::-1]
    else:
        coeff = np.polyfit(n, q, deg)

    coeffs = coeff[::-1]                      # ascending: a0, g_fit, alpha, beta
    named['coeffs'] = coeffs
    for key, idx in (('a0', 0), ('g_fit', 1), ('alpha', 2), ('beta', 3)):
        if idx < coeffs.size:
            named[key] = float(coeffs[idx])
    return named


def _polyval_asc(coeffs, n):
    """Evaluate a polynomial given ASCENDING coefficients."""
    return np.polyval(np.asarray(coeffs, dtype=float)[::-1], n)
