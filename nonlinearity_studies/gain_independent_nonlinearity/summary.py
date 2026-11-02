"""
Per-extension summary of every scalar the gain-independent analysis produces.
"""

import csv
from pathlib import Path

import numpy as np

from .peaks import n_max_per_ext


_FIELDS = [
    'ext', 'n_peaks', 'fit_n_max', 'peaks_dropped_left_of_zero',
    'g', 'a0', 'g_fit', 'alpha', 'G', 'G_minus_g', 'g_fit_minus_g',
    'delta_g_slope', 'delta_g_intercept', 'alpha_prime',
    'beta_cubic', 'cubic_minus_quad_max_adu',
    'gain_residual_integral', 'best_gain', 'best_i', 'nonlin_at_eval_charge',
]


def build_gain_indep_summary(q_ext, zero_idx_ext, g_ext, quad_fits_ext, gain_from_fit_ext,
                             delta_line_fits_ext, alpha_fits_ext, cumsum_ext,
                             cubic_trunc_ext=None, scan_ext=None, n_max=None):
    """One row per extension, in the order of ``_FIELDS``.

    ``n_max`` is the fit domain the fits were given -- a scalar, or the list the
    changepoint estimator returns -- and is recorded per extension so the row
    says how many peaks each number was fit over.
    """
    rows = []
    n_max_ext = n_max_per_ext(n_max, len(q_ext))
    for ext in range(len(q_ext)):
        quad = quad_fits_ext[ext]
        trunc = cubic_trunc_ext[ext] if cubic_trunc_ext is not None else {}
        scan = scan_ext[ext] if scan_ext is not None else {}
        cumsum = np.asarray(cumsum_ext[ext], dtype=float)
        rows.append({
            'ext': ext + 1,
            'n_peaks': len(q_ext[ext]),
            'fit_n_max': n_max_ext[ext] if n_max_ext[ext] is not None else len(q_ext[ext]),
            'peaks_dropped_left_of_zero': zero_idx_ext[ext],
            'g': g_ext[ext],
            'a0': quad['a0'],
            'g_fit': quad['g_fit'],
            'alpha': quad['alpha'],
            'G': gain_from_fit_ext[ext]['G'],
            'G_minus_g': gain_from_fit_ext[ext]['G_minus_g'],
            'g_fit_minus_g': gain_from_fit_ext[ext]['g_fit_minus_g'],
            'delta_g_slope': delta_line_fits_ext[ext]['slope'],
            'delta_g_intercept': delta_line_fits_ext[ext]['intercept'],
            'alpha_prime': alpha_fits_ext[ext]['slope'],
            'beta_cubic': trunc.get('beta', np.nan),
            'cubic_minus_quad_max_adu': trunc.get('max_diff_adu', np.nan),
            'gain_residual_integral': cumsum[-1] if cumsum.size else np.nan,
            'best_gain': scan.get('best_gain', np.nan),
            'best_i': scan.get('best_i', None),
            'nonlin_at_eval_charge': scan.get('best_nonlin', np.nan),
        })
    return rows


def format_gain_indep_summary(rows, eval_charge=500.0):
    """Readable text version of the summary rows."""
    lines = []
    for row in rows:
        lines.append(f"\nEXT {row['ext']}: {row['n_peaks']} peaks "
                     f"({row['peaks_dropped_left_of_zero']} dropped left of the zero peak)")
        lines.append(f"  fit domain      n_max      = {row['fit_n_max']} peaks")
        lines.append(f"  standard gain   g          = {row['g']:.6f} ADU/e-")
        lines.append(f"  parabola fit    a0         = {row['a0']:.6f} ADU")
        lines.append(f"                  g_fit      = {row['g_fit']:.6f} ADU/e-")
        lines.append(f"                  alpha      = {row['alpha']:.6e} ADU")
        lines.append(f"                  G = g_fit+alpha = {row['G']:.6f} ADU/e-"
                     f"   (G - g = {row['G_minus_g']:+.6f})")
        lines.append(f"  delta_g line    slope      = {row['delta_g_slope']:.6e} ADU/peak"
                     f"   (intercept {row['delta_g_intercept']:+.6f} ADU)")
        lines.append(f"  S(n) line       alpha'     = {row['alpha_prime']:.6e} ADU")
        lines.append(f"  cubic term      beta       = {row['beta_cubic']:.6e} ADU"
                     f"   (max |cubic - quad| = {row['cubic_minus_quad_max_adu']:.4f} ADU)")
        lines.append(f"  sum(g_fit(i)-g)            = {row['gain_residual_integral']:.4f} ADU")
        if row['best_i'] is not None:
            lines.append(f"  gain scan       best g_fit(i={row['best_i']}) = {row['best_gain']:.6f} ADU/e-"
                         f"   -> nonlinearity at {eval_charge:g} e- = {row['nonlin_at_eval_charge']:+.4f} e-")
    return '\n'.join(lines)


def write_gain_indep_summary(rows, save_path):
    """Write the summary rows to CSV."""
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    with save_path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f'Saved gain-independent summary to {save_path}')
