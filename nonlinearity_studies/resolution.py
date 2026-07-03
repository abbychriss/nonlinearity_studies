"""
resolution — split out of nonlinearity_studies.py.
"""
import csv
import numpy as np
from pathlib import Path
from scipy.optimize import curve_fit

from .fit_data import _gauss_single, _make_comb_model
from .summary import _normalize_charges

_RESOLUTION_FIELDS = ['ext', 'q', 'window', 'n_peaks_fit', 'expected_peaks',
                      'sigma_e', 'sigma_e_err', 'reduced_chi2', 'delta_aic', 'verdict']


def _resolution_poisson_weights(counts):
    """Per-bin Poisson sigma with a floor of 1 so empty bins don't blow up."""
    return np.sqrt(np.maximum(counts, 1.0))


def _resolution_reduced_chi2(counts, model, n_params):
    err = _resolution_poisson_weights(counts)
    chi2 = float(np.sum(((counts - model) / err) ** 2))
    dof = max(len(counts) - n_params, 1)
    return chi2, dof, chi2 / dof


def _resolution_poisson_deviance(counts, model):
    """Cash/Poisson deviance 2*sum(mu - n + n*ln(n/mu)); equals -2 ln L up to a
    data-only constant that cancels in delta-AIC, so it is valid for model
    comparison even when bin counts are low."""
    mu = np.maximum(model, 1e-12)
    n = np.asarray(counts, dtype=float)
    term = np.zeros_like(n)
    nz = n > 0
    term[nz] = n[nz] * np.log(n[nz] / mu[nz])
    return float(2.0 * np.sum(mu - n + term))


def resolution_at_charge(counts, centers, peaks, q, window, gain, s0):
    """Fit the constrained comb in [q - window/2, q + window/2] and score it.

    Args:
        counts, centers : the all-peaks histogram for one extension (electrons).
        peaks           : indices into centers of the detected peaks.
        q               : center charge of the window (electrons).
        window          : window width in electrons (~ number of peaks).
        gain            : ADU/e- for this extension.
        s0              : zero-peak sigma in ADU (double_gauss_popts[ext][0]);
                          s0/gain seeds the shared sigma in electrons.

    Returns a result dict (sigma_e, sigma_e_err, reduced_chi2, delta_aic,
    n_components, expected_peaks, plus the windowed data and fit curves for
    plotting). Raises ValueError if the window has too few bins or peaks to fit.
    """
    centers = np.asarray(centers, dtype=float)
    counts = np.asarray(counts, dtype=float)
    half = window / 2.0
    lo, hi = q - half, q + half

    bin_mask = (centers >= lo) & (centers <= hi)
    xw, yw = centers[bin_mask], counts[bin_mask]
    if xw.size < 5:
        raise ValueError(f'q={q}: only {xw.size} histogram bins in window '
                         f'[{lo:.0f}, {hi:.0f}] -- widen the window, raise bin_factor, '
                         f'or extend the all-peaks histogram range')

    peak_charges = centers[np.asarray(peaks, dtype=int)]
    means = peak_charges[(peak_charges >= lo) & (peak_charges <= hi)]
    if means.size < 2:
        raise ValueError(f'q={q}: only {means.size} detected peaks in window '
                         f'[{lo:.0f}, {hi:.0f}] -- nothing to decompose')

    bin_width = float(np.median(np.diff(xw)))
    sigma0 = max(s0 / gain, bin_width)

    # ---- comb fit: params = [sigma, A_0, ..., A_{K-1}] ----
    comb, K = _make_comb_model(means)
    amp0 = np.interp(means, xw, yw)
    amp0 = np.where(amp0 > 0, amp0, 1.0)
    p0 = [sigma0, *amp0]
    lower = [1e-6] + [0.0] * K
    upper = [window] + [np.inf] * K
    err = _resolution_poisson_weights(yw)
    popt, pcov = curve_fit(comb, xw, yw, p0=p0, sigma=err, absolute_sigma=True,
                           bounds=(lower, upper), maxfev=20000)
    sigma_e = float(popt[0])
    sigma_e_err = float(np.sqrt(pcov[0, 0])) if np.all(np.isfinite(pcov)) else np.nan
    comb_model = comb(xw, *popt)
    k_comb = K + 1

    # ---- null: single broad Gaussian ----
    p0_null = [float(np.max(yw)), float(np.average(xw, weights=np.maximum(yw, 1e-9))),
               max(window / 4.0, bin_width)]
    try:
        popt_null, _ = curve_fit(_gauss_single, xw, yw, p0=p0_null, sigma=err,
                                 absolute_sigma=True,
                                 bounds=([0.0, lo, 1e-6], [np.inf, hi, 5 * window]),
                                 maxfev=20000)
        null_model = _gauss_single(xw, *popt_null)
    except (RuntimeError, ValueError):
        popt_null, null_model = None, np.full_like(xw, np.mean(yw))
    k_null = 3

    # ---- scores ----
    chi2_comb, dof_comb, rchi2_comb = _resolution_reduced_chi2(yw, comb_model, k_comb)
    dev_comb = _resolution_poisson_deviance(yw, comb_model)
    dev_null = _resolution_poisson_deviance(yw, null_model)
    aic_comb = dev_comb + 2 * k_comb
    aic_null = dev_null + 2 * k_null
    nbin = len(yw)
    bic_comb = dev_comb + k_comb * np.log(nbin)
    bic_null = dev_null + k_null * np.log(nbin)

    # Integer electron peaks the window *should* contain if every electron
    # produced a resolvable peak. A find_peaks count well below this is itself a
    # resolution failure (the peaks were too smeared to detect).
    expected_peaks = int(np.floor(hi)) - int(np.ceil(lo)) + 1

    return {
        'q': float(q), 'window': float(window), 'lo': lo, 'hi': hi,
        'n_components': K, 'expected_peaks': expected_peaks, 'means': means,
        'sigma_e': sigma_e, 'sigma_e_err': sigma_e_err, 'sigma_seed_e': sigma0,
        'reduced_chi2': rchi2_comb, 'chi2': chi2_comb, 'dof': dof_comb,
        'delta_aic': aic_null - aic_comb,   # > 0 favors the comb (resolved)
        'delta_bic': bic_null - bic_comb,
        'xw': xw, 'yw': yw, 'bin_width': bin_width,
        'comb_popt': popt, 'comb_model': comb_model,
        'null_popt': popt_null, 'null_model': null_model,
    }


def classify_resolution(res, sigma_well=0.25, sigma_limit=0.5, min_peak_frac=0.6):
    """Three-tier resolution verdict from the fitted sigma and the peak count.

    Two independent signals, either of which can condemn an extension:

      * sigma (e-) relative to the 1 e- peak spacing. Unit-spaced Gaussians show
        a central valley only when sigma < ~0.5 e- (separation > 2 sigma) and
        separate cleanly when sigma < ~0.25 e- (separation > 4 sigma).
      * the fraction of expected peaks actually detected/fit. If find_peaks only
        recovered a handful of the ~window+1 peaks the window should hold, the
        peaks were too smeared to detect -- a resolution failure regardless of
        what sigma the (under-determined) fit returned. This catches an
        extension that fit to only 2 peaks.

    Returns one of 'well resolved', 'marginal', 'unresolved'.
    """
    frac = res['n_components'] / max(res['expected_peaks'], 1)
    if frac < min_peak_frac or res['sigma_e'] >= sigma_limit:
        return 'unresolved'
    if res['sigma_e'] < sigma_well:
        return 'well resolved'
    return 'marginal'


def resolution_at_charge_ext(counts_ext, centers_ext, peaks_ext, gains, double_gauss_popts,
                             charges, window=10.0,
                             sigma_well=0.25, sigma_limit=0.5, min_peak_frac=0.6,
                             verbose=True):
    """Per-extension single-electron resolution at one or more charges.

    Returns ``results_ext``: a list (one per extension) of lists (one per charge)
    of result dicts (or None where the window could not be fit), each augmented
    with 'ext' and 'verdict'.
    """
    charges = _normalize_charges(charges)
    results_ext = []
    for ext in range(len(peaks_ext)):
        s0 = double_gauss_popts[ext][0]
        rows = []
        for q in charges:
            try:
                res = resolution_at_charge(counts_ext[ext], centers_ext[ext],
                                           peaks_ext[ext], q, window, gains[ext], s0)
            except ValueError as exc:
                if verbose:
                    print(f'EXT {ext + 1}, q={q:g}: {exc}')
                rows.append(None)
                continue
            res['ext'] = ext + 1
            res['verdict'] = classify_resolution(res, sigma_well, sigma_limit, min_peak_frac)
            rows.append(res)
        results_ext.append(rows)
    return results_ext


def _resolution_record(res):
    return {
        'ext': res['ext'], 'q': res['q'], 'window': res['window'],
        'n_peaks_fit': res['n_components'], 'expected_peaks': res['expected_peaks'],
        'sigma_e': res['sigma_e'], 'sigma_e_err': res['sigma_e_err'],
        'reduced_chi2': res['reduced_chi2'], 'delta_aic': res['delta_aic'],
        'verdict': res['verdict'],
    }


def format_resolution_table(results_ext):
    """Format resolution results (from resolution_at_charge_ext) as an aligned table."""
    headers = ['EXT', 'q [e-]', 'win', 'peaks', 'sigma [e-]', 'sig_err',
               'red_chi2', 'dAIC', 'verdict']
    cells = []
    for rows in results_ext:
        for res in rows:
            if res is None:
                continue
            cells.append([
                f'{res["ext"]}', f'{res["q"]:.0f}', f'{res["window"]:.0f}',
                f'{res["n_components"]}/{res["expected_peaks"]}',
                f'{res["sigma_e"]:.3f}', f'{res["sigma_e_err"]:.3f}',
                f'{res["reduced_chi2"]:.2f}', f'{res["delta_aic"]:.0f}',
                res['verdict'],
            ])
    widths = [max(len(h), *(len(row[i]) for row in cells)) if cells else len(h)
              for i, h in enumerate(headers)]
    fmt = lambda vals: '  '.join(v.rjust(widths[i]) for i, v in enumerate(vals))
    lines = [fmt(headers), '-' * (sum(widths) + 2 * (len(widths) - 1))]
    lines += [fmt(c) for c in cells]
    return '\n'.join(lines)


def summarize_resolution(results_ext, save_path=None):
    """Print a per-extension/charge resolution table, optionally saving it as CSV."""
    print('\n***********************************************************')
    print('\nSingle-electron resolution summary:\n')
    print(format_resolution_table(results_ext))

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        with save_path.open('w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=_RESOLUTION_FIELDS)
            writer.writeheader()
            for rows in results_ext:
                for res in rows:
                    if res is not None:
                        writer.writerow(_resolution_record(res))
        print(f'\nSaved resolution summary table (CSV) to {save_path}')

    return results_ext
