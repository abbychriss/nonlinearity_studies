"""
Gain optimization: which single gain minimizes the standard nonlinearity?

The standard analysis divides the whole spectrum by one gain g measured from
the zero and one electron peaks. Here that step is repeated with every
candidate gain from the running fit g_fit(i) -- i.e. the average gain over the
first i peaks, from i=2 (which reproduces the standard gain exactly) upwards --
and the standard nonlinearity is recomputed each time with the production
routine ``fit_data.fit_nonlinearity``, unchanged. The candidate whose
nonlinearity at ``eval_charge`` is smallest in absolute value is the gain that
leaves the least analysis-induced nonlinearity behind.
"""

import numpy as np

from ..fit_data import fit_nonlinearity, parabola


def scan_gain_nonlinearity_ext(q_ext, running_ext, fit_range_right=1000.0,
                               eval_charge=500.0,
                               fit_bounds_low=-100, fit_bounds_high=100,
                               fix_intercept=True):
    """Calibrate the ADU spectrum with each candidate gain and refit the nonlinearity.

    Args:
        q_ext: list (one per extension) of ADU charges q(n), n = 0, 1, 2, ...
        running_ext: output of ``spectrum_fit.fit_charge_spectrum_running_ext``;
            its ``i`` and ``g_fit`` arrays supply the candidate gains.
        fit_range_right: right edge (in electrons) of the parabola fit, passed
            straight to ``fit_nonlinearity``. A single number applies to every
            extension; a list gives one right edge per extension (what
            ``--gain_indep_scan_fit_range_right auto`` produces).
        eval_charge: charge (in electrons) the nonlinearity is reported at.

    Returns:
        list (one per extension) of dicts with arrays ``i``, ``gain``,
        ``nonlin`` (nonlinearity at ``eval_charge``, in electrons), plus the
        scalars ``best_i``, ``best_gain``, ``best_nonlin`` for the candidate
        with the smallest |nonlin|, and ``g_dir_eval`` / ``nonlin_dir`` -- the
        direct gain g_dir(eval_charge) and the nonlinearity it leaves behind,
        for comparison against the scanned candidates.
    """
    out = []

    fit_range_right_ext = (list(fit_range_right)
                           if isinstance(fit_range_right, (list, tuple))
                           else [fit_range_right] * len(q_ext))

    for fit_right, q, run in zip(fit_range_right_ext, q_ext, running_ext):
        q = np.asarray(q, dtype=float)
        i_values = np.asarray(run['i'])
        gains = np.asarray(run['g_fit'], dtype=float)
        nonlin = np.full(gains.size, np.nan)
        # Peak numbers are the x-axis of the standard nonlinearity, so the
        # spectrum is referenced to the zero-electron peak before dividing --
        # exactly what convert_to_electrons does with the pedestal.
        q_rel = q - q[0] if q.size else q

        def _nonlin_at(gain):
            """Nonlinearity at eval_charge when this extension is calibrated by gain."""
            if not np.isfinite(gain) or gain <= 0:
                return np.nan
            peak_charge_e = q_rel / gain
            try:
                coeff, _, _, _ = fit_nonlinearity(
                    np.arange(peak_charge_e.size), peak_charge_e, 0.0, 1.0,
                    fit_right, do_convert_to_electrons=False,
                    fit_bounds_low=fit_bounds_low, fit_bounds_high=fit_bounds_high,
                    fix_intercept=fix_intercept)
            except (ValueError, RuntimeError):
                return np.nan
            return float(parabola(float(eval_charge), *coeff))

        for k, gain in enumerate(gains):
            nonlin[k] = _nonlin_at(gain)

        # The direct gain out to the peak eval_charge itself sits at, which the
        # scan never visits: it comes from the raw spectrum, not the running fit.
        n_eval = int(round(eval_charge))
        g_dir_eval = float(q_rel[n_eval] / n_eval) if q_rel.size > n_eval > 0 else np.nan

        best = int(np.nanargmin(np.abs(nonlin))) if np.any(np.isfinite(nonlin)) else None
        out.append(dict(
            i=i_values, gain=gains, nonlin=nonlin,
            best_i=int(i_values[best]) if best is not None else None,
            best_gain=float(gains[best]) if best is not None else np.nan,
            best_nonlin=float(nonlin[best]) if best is not None else np.nan,
            g_dir_eval=g_dir_eval, nonlin_dir=_nonlin_at(g_dir_eval),
        ))

    return out
