import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import biweight_location, biweight_midvariance
from scipy.optimize import curve_fit
from scipy.signal import find_peaks as scipy_find_peaks
import math
import csv

from pathlib import Path
from glob import glob
from tqdm import tqdm

# --- Pedestal subtraction + double-Gaussian zero/one fitting/plotting ---
# These now live in the standalone `pedestal_subtract` package. They are imported
# here (verbatim source) so this module's public names stay unchanged and the rest
# of the pipeline keeps working. Edit them in the pedestal_subtract package.
from pedestal_subtract.core import (
    convert_to_electrons,
    _smooth_counts,
    _make_histogram,
    _estimate_peak_width,
    _clip_to_bounds,
    _auto_zero_one_setup,
    calculate_noise_gain,
    double_gauss,
    pedestal_subtract,
    _PEDSUB_ALGO_VERSION,
    _PEDSUB_HEADER_KEYS,
    _pedsub_cache_path,
    _pedsub_header_matches,
    pedestal_subtract_ext_cached,
    _finish_fig,
    _fit_double_gauss_electrons,
    plot_zero_one_peaks,
    get_fits,
    _value_for_extension,
    get_zero_one_peaks_ext,
)

#---------------- ANALYSIS FUNCTIONS ----------------------------


#---------------- (2) Find peaks ----------------------------
# Finds all electron peaks in flattened pixel charge array per extension
# First converts data from ADU to electrons (if not already in electrons)
# Input is charge data (in ADU or electrons) from one extension
# bins is the number of bins given for initial charge histogram
# bin_factor is the multiple of the length of range used in fitting all peaks (number of bins per peak essentially)
# bin_factor is also used to define the distance parameter given to scipy_find_peaks 
# for the min distance between peaks (with buffer given by buffer) 
def find_all_peaks(data,
                   width,
                   buffer,
                   pedestal,
                   noise,
                   gain,
                   bins='default',
                   flatten=True,
                   do_convert_to_electrons=True,
                   range_left='left_of_zero',
                   range_right=2500,
                   bin_factor=8,
                   prominence=None):

    if flatten:
        data=np.array(data).flatten()

    if do_convert_to_electrons:
        data = convert_to_electrons(data, pedestal, gain, flatten=False)

    if range_left=='default':
        if do_convert_to_electrons:
            range_left = -0.5
        else:
            range_left=pedestal-3*noise-1 # If finding electron peaks in ADU, start fitting peaks at the left side of the zero electron peak

    hist_range = (range_left, range_right)

    if bins=='default':
       bins=math.floor((hist_range[1]-hist_range[0])*bin_factor)

    counts, edges = np.histogram(data, bins=bins,range=hist_range)
    centers = 0.5 * (edges[1:] + edges[:-1])
    peaks, properties = scipy_find_peaks(counts, height=0, width=width*bin_factor, distance=bin_factor-buffer, prominence=prominence)

    return counts, edges, peaks, centers, properties, hist_range


#---------------- (3) Fit nonlinearity ----------------------------
def fit_nonlinearity(peaks, centers, pedestal, gain, fit_range_right, do_convert_to_electrons=False, fit_bounds_low=-100, fit_bounds_high=100):
    # If specified, convert to electrons: subtract the pedestal (mean of zero electron peak) from all charge values, 
    # Then divide by the gain (difference between zero and 1 electron peak)
    if do_convert_to_electrons:
        peak_charge_e = np.array([(centers[p]-pedestal)/gain for p in peaks])

    # If peaks were given in electrons, do not convert to electrons. Just pick out electron values of peak locations.
    else:
        peak_charge_e = np.array([centers[p] for p in peaks])

    charge_minus_npeak = [(peak_charge_e[i] - i) for i in range(len(peaks))]
    fit_idx = int(np.searchsorted(peak_charge_e, fit_range_right))
    if fit_idx < 3:
        raise ValueError(
            f'Not enough peaks below fit_range_right={fit_range_right} to fit a parabola '
            f'(found {fit_idx}, need >=3). Total peaks detected: {len(peaks)}. '
            f'Loosen the peak-finder filters (lower peak_finder_widths, raise peak_finder_buffers, '
            f'lower peak_finder_prominences) or raise fit_range_right.'
        )
    parabola_coeff, parabola_pcov = curve_fit(parabola, peak_charge_e[:fit_idx], charge_minus_npeak[:fit_idx],
                           maxfev=2000, bounds=(fit_bounds_low, fit_bounds_high))
    return parabola_coeff, parabola_pcov, peak_charge_e, charge_minus_npeak

#---------------- (3b) Estimate optimal fit_range_right --------------------------
def estimate_optimal_fit_range_right(peaks, centers, pedestal, gain,
                                     min_charge=100, max_charge=None, n_candidates=40,
                                     do_convert_to_electrons=False,
                                     fit_bounds_low=-100, fit_bounds_high=100,
                                     tolerance=None, max_charge_percentile=None,
                                     min_peaks_in_fit=None):
    """Sweep candidate fit_range_right values and return the one that minimizes var(a).

    Score is pcov[0, 0] from the parabola curve_fit — the variance of the leading
    (curvature) coefficient. Small fit_range_right -> few points -> high var(a);
    large fit_range_right -> model breaks down at high charge -> high var(a). The
    minimum is the largest charge value where the parabola still fits cleanly.

    Robustness options (None = off, current behavior):
        tolerance (float):  Pick the smallest candidate whose score is within
                            (1 + tolerance) * min_score. Prevents overestimation
                            when the score curve is flat or noisy.
        max_charge_percentile (float): Cap max_charge at this percentile of the
                            detected peak charges (e.g. 90 keeps candidates out
                            of the noisy upper tail).
        min_peaks_in_fit (int):  Reject candidates that include fewer than this
                            many peaks below them. Prevents the estimator from
                            picking absurdly small ranges that minimize var(a)
                            with too few data points.
    """
    if do_convert_to_electrons:
        peak_charge_e = np.array([(centers[p] - pedestal) / gain for p in peaks])
    else:
        peak_charge_e = np.array([centers[p] for p in peaks])

    if len(peak_charge_e) < 5:
        raise ValueError(f'Need at least 5 peaks to estimate fit_range_right, got {len(peak_charge_e)}')

    if max_charge_percentile is not None:
        max_charge_pct = float(np.percentile(peak_charge_e, float(max_charge_percentile)))
        max_charge = max_charge_pct if max_charge is None else min(float(max_charge), max_charge_pct)
    elif max_charge is None:
        max_charge = float(peak_charge_e[-1])

    min_charge = max(float(min_charge), float(peak_charge_e[3]))
    if max_charge <= min_charge:
        raise ValueError(f'max_charge ({max_charge}) must exceed min_charge ({min_charge})')

    candidates = np.linspace(min_charge, max_charge, int(n_candidates))
    min_peaks_required = int(min_peaks_in_fit) if min_peaks_in_fit is not None else 0
    results = []
    for candidate in candidates:
        if min_peaks_required > 0:
            n_peaks_below = int(np.searchsorted(peak_charge_e, float(candidate)))
            if n_peaks_below < min_peaks_required:
                continue
        try:
            _, pcov, _, _ = fit_nonlinearity(peaks, centers, pedestal, gain, float(candidate),
                                              do_convert_to_electrons=do_convert_to_electrons,
                                              fit_bounds_low=fit_bounds_low,
                                              fit_bounds_high=fit_bounds_high)
        except (ValueError, RuntimeError, np.linalg.LinAlgError):
            continue
        score = pcov[0, 0]
        if not np.isfinite(score) or score <= 0:
            continue
        results.append((float(candidate), float(score)))

    if not results:
        raise ValueError('Could not find any valid fit_range_right candidate; try lowering min_charge or loosening peak finder settings')

    if tolerance is not None and tolerance > 0:
        # Among all candidates within tolerance of the min score, prefer the smallest charge.
        min_score = min(s for _, s in results)
        threshold = min_score * (1 + float(tolerance))
        acceptable = [(c, s) for c, s in results if s <= threshold]
        return min(acceptable, key=lambda x: x[0])[0]

    return min(results, key=lambda x: x[1])[0]


def estimate_optimal_fit_range_right_ext(peaks_ext, centers_ext, pedestals, gains,
                                          min_charge=100, max_charges=None, n_candidates=40,
                                          do_convert_to_electrons=False,
                                          fit_bounds_low=-100, fit_bounds_high=100,
                                          tolerance=None, max_charge_percentile=None,
                                          min_peaks_in_fit=None, verbose=True):
    n_ext = len(peaks_ext)
    if max_charges is None:
        max_charges = [None] * n_ext
    elif not isinstance(max_charges, (list, tuple)):
        max_charges = [max_charges] * n_ext

    result = []
    for i in range(n_ext):
        opt = estimate_optimal_fit_range_right(
            peaks_ext[i], centers_ext[i], pedestals[i], gains[i],
            min_charge=min_charge, max_charge=max_charges[i], n_candidates=n_candidates,
            do_convert_to_electrons=do_convert_to_electrons,
            fit_bounds_low=fit_bounds_low, fit_bounds_high=fit_bounds_high,
            tolerance=tolerance, max_charge_percentile=max_charge_percentile,
            min_peaks_in_fit=min_peaks_in_fit,
        )
        result.append(opt)
        if verbose:
            print(f'EXT {i+1}: estimated optimal fit_range_right = {opt:.0f} e-')
    return result


def estimate_fit_range_right_by_noise_onset(peaks, centers, pedestal, gain,
                                            min_charge=100, do_convert_to_electrons=False,
                                            window=30, factor=2.5):
    """Detect where the nonlinearity curve starts becoming noisy.

    For each peak i, fit a local quadratic in a window of ``window`` peaks centered
    at i and compute the residual MAD. Establish a baseline MAD from the early
    "safe" region (below the 30th percentile of charges above ``min_charge``).
    Walk forward; return the first peak charge where the rolling MAD exceeds
    ``factor * baseline``.

    This is an alternative to ``estimate_optimal_fit_range_right`` (which minimizes
    var(a)). It tries to reproduce a visual inspection: "where does the curve get
    jagged?"
    """
    if do_convert_to_electrons:
        peak_charge_e = np.array([(centers[p] - pedestal) / gain for p in peaks])
    else:
        peak_charge_e = np.array([centers[p] for p in peaks])
    charge_minus_npeak = np.array([peak_charge_e[i] - i for i in range(len(peaks))])

    n = len(peak_charge_e)
    if n < window + 10:
        raise ValueError(f'Need at least {window + 10} peaks for noise-onset detection, got {n}')

    half = int(window) // 2
    rolling_mad = np.full(n, np.nan)
    for i in range(half, n - half):
        x = peak_charge_e[i - half:i + half + 1]
        y = charge_minus_npeak[i - half:i + half + 1]
        coef = np.polyfit(x, y, 2)
        residuals = y - np.polyval(coef, x)
        rolling_mad[i] = np.median(np.abs(residuals - np.median(residuals)))

    above_min = peak_charge_e >= float(min_charge)
    if not above_min.any():
        raise ValueError(f'No peaks above min_charge={min_charge}')
    cutoff = np.percentile(peak_charge_e[above_min], 30)
    safe_mask = above_min & (peak_charge_e <= cutoff) & np.isfinite(rolling_mad)
    safe_idx = np.where(safe_mask)[0]
    if len(safe_idx) < 5:
        raise ValueError('Not enough peaks in safe region to establish baseline MAD')
    baseline_mad = float(np.median(rolling_mad[safe_idx]))
    if baseline_mad <= 0:
        baseline_mad = float(np.mean(rolling_mad[safe_idx]))
    threshold = float(factor) * baseline_mad

    start = int(safe_idx[-1]) + 1
    for i in range(start, n - half):
        if np.isfinite(rolling_mad[i]) and rolling_mad[i] > threshold:
            return float(peak_charge_e[i])

    return float(peak_charge_e[-1])


def estimate_fit_range_right_by_noise_onset_ext(peaks_ext, centers_ext, pedestals, gains,
                                                min_charge=100, do_convert_to_electrons=False,
                                                window=30, factor=2.5, verbose=True):
    result = []
    for i in range(len(peaks_ext)):
        try:
            opt = estimate_fit_range_right_by_noise_onset(
                peaks_ext[i], centers_ext[i], pedestals[i], gains[i],
                min_charge=min_charge, do_convert_to_electrons=do_convert_to_electrons,
                window=window, factor=factor,
            )
        except ValueError as exc:
            raise ValueError(f'EXT {i+1}: {exc}') from exc
        result.append(opt)
        if verbose:
            print(f'EXT {i+1}: noise-onset fit_range_right = {opt:.0f} e- (window={window}, factor={factor})')
    return result


#---------------- (3c) Estimate fit_range_right by changepoint (default) ----------
def _local_roughness(x, y, win):
    """Centered local robust scatter: 1.4826 * MAD of residuals from a local
    quadratic fit in a window of ``win`` points.

    Because the quadratic absorbs the local curvature, this is flat through the
    whole clean parabola (regardless of how curved it is) and steps up sharply
    where the data turns noisy. Being local, it never accumulates the
    global-extrapolation drift that fooled covariance-based scores.
    """
    n = len(x)
    half = int(win) // 2
    rough = np.full(n, np.nan)
    for i in range(half, n - half):
        sl = slice(i - half, i + half + 1)
        c = np.polyfit(x[sl], y[sl], 2)
        r = y[sl] - np.polyval(c, x[sl])
        rough[i] = 1.4826 * np.median(np.abs(r - np.median(r)))
    return rough


def estimate_fit_range_right_changepoint(peaks, centers, pedestal, gain,
                                         do_convert_to_electrons=False,
                                         win=25, baseline_pct=25, factor=4.0,
                                         k=4.0, persist=4, floor=0.15,
                                         return_debug=False):
    """Find the charge where the clean parabola gives way to noise.

    This is the default ``fit_range_right`` estimator. It is a two-stage,
    drift-free changepoint detector on the nonlinearity curve (peak charge vs.
    charge-minus-peak-number):

    Stage 1 (coarse, robust) -- centered LOCAL roughness (see _local_roughness)
        is flat through the whole clean parabola and steps up at the noise
        onset. A baseline is taken from the early clean region (charges below
        ``baseline_pct`` percentile); the first center index where roughness
        stays above ``max(factor * baseline, floor)`` for ``persist``
        consecutive points marks "noise is near here". The absolute ``floor``
        (e-) keeps an ultra-quiet baseline from setting an impossibly tight
        threshold that would trip on sub-noise waviness.

    Stage 2 (precise onset) -- everything up to ``coarse - win//2`` is
        guaranteed clean (the roughness window had not yet reached noise). A
        parabola is fit to those points, ``sigma_clean`` is estimated (robust),
        and the cut is placed at the first run of ``persist`` points whose
        absolute residual exceeds ``max(k * sigma_clean, floor)``. The short
        refinement span avoids both extrapolation drift and the half-window
        underestimate of a centered detector.

    Compared with ``estimate_optimal_fit_range_right`` (minimize var(a)), this
    is immune to the two var(a) failure modes: a flat score curve collapsing to
    the smallest range, and a near-linear curve whose score decreases
    monotonically into the noisy tail.
    """
    if do_convert_to_electrons:
        peak_charge_e = np.array([(centers[p] - pedestal) / gain for p in peaks], dtype=float)
    else:
        peak_charge_e = np.array([centers[p] for p in peaks], dtype=float)
    x = peak_charge_e
    y = np.array([peak_charge_e[i] - i for i in range(len(peaks))], dtype=float)

    n = len(x)
    half = int(win) // 2
    if n < win + persist + 4:
        raise ValueError(f'Need at least {win + persist + 4} peaks for changepoint '
                         f'detection, got {n}')

    rough = _local_roughness(x, y, win)

    qcut = np.percentile(x, baseline_pct)
    base_mask = (x <= qcut) & np.isfinite(rough)
    if base_mask.sum() < 3:
        base_mask = np.isfinite(rough)
    baseline = float(np.median(rough[base_mask]))
    if not np.isfinite(baseline) or baseline <= 0:
        baseline = float(np.nanmedian(rough)) or 1e-9
    thr = max(factor * baseline, float(floor))

    # ---- Stage 1: coarse crossing with persistence ----
    start = max(int(np.searchsorted(x, qcut)), half)
    coarse = None
    for i in range(start, n - half - persist):
        seg = rough[i:i + persist]
        if np.all(np.isfinite(seg)) and np.all(seg > thr):
            coarse = i
            break
    if coarse is None:
        # No sustained noise onset detected -> the curve is clean to the end.
        if return_debug:
            return float(x[-1]), dict(rough=rough, baseline=baseline, thr=thr,
                                      coarse_idx=None, clean_end=None,
                                      sigma_clean=None, res_thr=None, break_idx=n - 1)
        return float(x[-1])

    # ---- Stage 2: precise onset refinement on guaranteed-clean points ----
    clean_end = max(coarse - half, 8)
    coef = np.polyfit(x[:clean_end], y[:clean_end], 2)
    resid = y[:clean_end] - np.polyval(coef, x[:clean_end])
    sigma_clean = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    if not np.isfinite(sigma_clean) or sigma_clean <= 0:
        sigma_clean = float(np.std(resid)) or 1e-9

    res_thr = max(k * sigma_clean, float(floor))
    absres = np.abs(y - np.polyval(coef, x))
    break_idx = n - 1
    for i in range(clean_end, n - persist):
        if np.all(absres[i:i + persist] > res_thr):
            break_idx = i
            break

    if return_debug:
        return float(x[break_idx]), dict(rough=rough, baseline=baseline, thr=thr,
                                         coarse_idx=coarse, clean_end=clean_end,
                                         sigma_clean=sigma_clean, res_thr=res_thr,
                                         break_idx=break_idx)
    return float(x[break_idx])


def estimate_fit_range_right_changepoint_ext(peaks_ext, centers_ext, pedestals, gains,
                                             do_convert_to_electrons=False,
                                             win=25, baseline_pct=25, factor=4.0,
                                             k=4.0, persist=4, floor=0.15,
                                             cross_check=True, confidence_rel_tol=0.15,
                                             verbose=True):
    """Per-extension changepoint estimate of fit_range_right.

    Returns ``(result, diagnostics)`` where ``result`` is the list of
    fit_range_right values (one per extension) and ``diagnostics`` is a list of
    per-extension dicts. When ``cross_check`` is on, each diagnostic also holds
    the independent var(a) estimate and a confidence label: extensions where the
    two methods disagree by more than ``confidence_rel_tol`` (relative to the
    changepoint value) are flagged 'LOW' for review.
    """
    result = []
    diagnostics = []
    for i in range(len(peaks_ext)):
        try:
            cp = estimate_fit_range_right_changepoint(
                peaks_ext[i], centers_ext[i], pedestals[i], gains[i],
                do_convert_to_electrons=do_convert_to_electrons,
                win=win, baseline_pct=baseline_pct, factor=factor,
                k=k, persist=persist, floor=floor,
            )
        except ValueError as exc:
            raise ValueError(f'EXT {i+1}: {exc}') from exc

        d = {'ext': i + 1, 'fit_range_right': cp, 'method': 'changepoint'}
        if cross_check:
            try:
                va = estimate_optimal_fit_range_right(
                    peaks_ext[i], centers_ext[i], pedestals[i], gains[i],
                    do_convert_to_electrons=do_convert_to_electrons,
                )
                rel = abs(cp - va) / cp if cp else float('inf')
                d['var_a_cross_check'] = va
                d['rel_diff'] = rel
                d['confidence'] = 'LOW' if rel > confidence_rel_tol else 'OK'
            except Exception as exc:  # cross-check is best-effort
                d['var_a_cross_check'] = None
                d['rel_diff'] = None
                d['confidence'] = 'UNKNOWN'
                d['cross_check_error'] = str(exc)

        result.append(cp)
        diagnostics.append(d)
        if verbose:
            msg = f'EXT {i+1}: changepoint fit_range_right = {cp:.0f} e-'
            if cross_check:
                va = d.get('var_a_cross_check')
                if va is not None:
                    msg += (f'  (var(a) cross-check = {va:.0f} e-, '
                            f'rel diff = {d["rel_diff"]*100:.0f}%, '
                            f'confidence = {d["confidence"]})')
                    if d['confidence'] == 'LOW':
                        msg += '  <-- REVIEW'
                else:
                    msg += f'  (var(a) cross-check failed: {d.get("cross_check_error", "")})'
            print(msg)
    return result, diagnostics


#---------------- (4) Get nonlinearity at charge value ----------------------------
# use coefficients of parabola fit to nonlinearity curve from fit_nonlinearity to estimate the nonlinearity
# at a single value or list of values
# required arguments: q (charge value that parabola equation is evaluated at, can be int, float or list), parabola_coeff (the curve_fit optimal parabola coefficients found in fit_nonlinearity)
# optional arguments: parabola_pcov, fit_range_right (used to ascertain how confident we should be in the nonlinearity value); print_values (if you want to print what the function outputs)
def get_nonlinearity_at(q, parabola_coeff, parabola_pcov=None, fit_range_right=None, print_values=False, out_lines=None):
    a = parabola_coeff[0]
    b = parabola_coeff[1]
    c = parabola_coeff[2]

    def _emit(s):
        print(s)
        if out_lines is not None:
            out_lines.append(s)

    _emit(f'Parabola best fit: Nonlin(q) = {np.round(a,5)}*q^2 + {np.round(b,5)}*q + {np.round(c,5)}\n')

    if isinstance(q, (float, int, np.floating, np.integer)):
        nonlinearity_at_q = np.round(a*q**2 + b*q + c, 3)
    elif isinstance(q, (list, np.ndarray)):
        nonlinearity_at_q = [np.round((a*q_i**2 + b*q_i + c), 3) for q_i in q]
    else:
        raise TypeError(f'q must be a number or list, got {type(q)}')

    if isinstance(q, (float, int, np.floating, np.integer)):
        _emit(f'Interpolated nonlinearity at {q} e- = {nonlinearity_at_q} e-\n')
    else:
        for i in range(len(q)):
            _emit(f'Interpolated nonlinearity at {q[i]} e- = {nonlinearity_at_q[i]} e-\n')

    if print_values:
        if fit_range_right is not None:
            _emit(f'Parabola fit was performed up to {fit_range_right} e-\n')
        if parabola_pcov is not None:
            _emit(f'Covariance of fit = {parabola_pcov}\n')

    return nonlinearity_at_q


#---------------- (5) Single-electron resolution at a charge ----------------------------
# Quantify how well single-electron peaks are resolved near a charge q by fitting a
# constrained sum of Gaussians ("comb") to the windowed charge spectrum:
#   - means FIXED at the already-detected peak locations (absorbs the nonlinearity shift,
#     which moves absolute peak positions but barely changes the local ~1 e- spacing),
#   - a single SHARED sigma across components (single-electron readout noise is ~constant
#     across a small window); this sigma, in electrons, IS the resolution figure of merit,
#   - free amplitudes (they follow the charge-distribution envelope).
# Goodness of fit is reported as the comb's reduced chi^2 and as delta-AIC vs a single
# broad Gaussian "unresolved" null (positive favors the comb -> peaks are resolved).

def _gauss_single(x, amp, mu, sigma):
    return amp * np.exp(-0.5 * ((np.asarray(x, dtype=float) - mu) / sigma) ** 2)


def _make_comb_model(means):
    """Return a curve_fit-compatible model f(x, sigma, A_0, A_1, ...) that sums
    len(means) Gaussians with FIXED means and a single SHARED sigma."""
    means = np.asarray(means, dtype=float)

    def comb(x, sigma, *amps):
        x = np.asarray(x, dtype=float)
        y = np.zeros_like(x)
        for amp, mu in zip(amps, means):
            y += amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
        return y

    return comb, len(means)


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


_RESOLUTION_FIELDS = ['ext', 'q', 'window', 'n_peaks_fit', 'expected_peaks',
                      'sigma_e', 'sigma_e_err', 'reduced_chi2', 'delta_aic', 'verdict']


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


#---------------- PLOTTING FUNCTIONS ----------------------------


#---------------- Plot all electron peaks ----------------------------
# Input is list of data from each of four extensions
# ylim can be 'none' or tuple=(ylim_bottom, ylim_top)
def plot_all_peaks(counts_ext, 
                   peaks_ext, 
                   centers_ext, 
                   xlim, ylim='none', 
                   yscale='log', 
                   plot_individual=True, plot_together=False,
                   draw_lines=True, linecolor='r', linestyle='--',
                   peak_number_labels_individual=True, peak_number_labels_together=True, peak_number_label_size=8,
                   individual_figsize=(6,5), subplots_figsize=(9,7),
                   additional_title='',
                   suptitle='Peaks in Pixel Charge Distribution',
                   nimages=10,
                   sharex=True,
                   sharey=True,
                   show_titles=True,
                   save_plots=False,
                   show_plots=True,
                   fig_path='./', file='peak_finder',
                   dpi=350):

    fig_path = Path(fig_path)
    if file != 'peak_finder':
        base_name = file[:-5] + '_peak_finder'
    else:
        base_name = file
    fig_name = fig_path / base_name

    if plot_individual:
        for ext, counts in enumerate(counts_ext):
            peaks=peaks_ext[ext]
            centers=centers_ext[ext]
            bin_width = centers[1] - centers[0]
            
            fig, ax = plt.subplots(1, 1, figsize=individual_figsize, constrained_layout=True)
            if show_titles:
                ax.set_title(f'{additional_title}{suptitle} (Nimages = {nimages}): EXT {ext + 1}', fontsize=12, pad=10)
            ax.bar(centers, counts, align='center', edgecolor='none', linewidth=0, width=bin_width)
            ax.set_xlabel(r'Charge ($e^-$)')
            ax.set_ylabel('N')
            if yscale!='linear':
                ax.set_yscale(yscale)
            xlim_e = xlim[ext] if isinstance(xlim, list) else xlim
            ylim_e = ylim[ext] if isinstance(ylim, list) else ylim
            if xlim_e!='none':
                ax.set_xlim(xlim_e)
            if ylim_e!='none':
                ax.set_ylim(ylim_e)

            # draw vertical lines and labels at each peak
            if draw_lines:
                for i,p in enumerate(peaks):
                    peak_x = centers[p]

                    ax.axvline(peak_x, linestyle=linestyle, color=linecolor)
                    if peak_number_labels_individual:
                        label_y = 0.93 if (i % 2 == 0) else 0.83
                        ax.text(peak_x,
                            label_y,
                            f"{i}",
                            transform=ax.get_xaxis_transform(),
                            verticalalignment='top',
                            horizontalalignment='center',
                            color=linecolor,
                            fontsize=peak_number_label_size)

            if save_plots:
                output_path = fig_name.with_stem(fig_name.stem + f'_EXT{ext+1}').with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            _finish_fig(show_plots)

    if plot_together:
        fig, axs = plt.subplots(2,2,figsize=subplots_figsize,constrained_layout=True,sharex=sharex,sharey=sharey)
        axs=axs.flatten()
        if show_titles:
            fig.suptitle(f'{additional_title}{suptitle} (Nimages = {nimages})')

        for ext, counts in enumerate(counts_ext):
            peaks=peaks_ext[ext]
            centers=centers_ext[ext]
            bin_width = centers[1] - centers[0]
            ax = axs[ext]

            ax.bar(centers, counts, align='center', edgecolor='none', linewidth=0, width=bin_width)
            ax.set_xlabel(r'Charge ($e^-$)')
            ax.set_ylabel('N')
            if yscale!='linear':
                ax.set_yscale(yscale)
            xlim_e = xlim[ext] if isinstance(xlim, list) else xlim
            ylim_e = ylim[ext] if isinstance(ylim, list) else ylim
            if xlim_e!='none':
                ax.set_xlim(xlim_e)
            if ylim_e!='none':
                ax.set_ylim(ylim_e)
            ax.set_title(f'EXT {ext + 1}')

            # draw vertical lines and labels at each peak
            if draw_lines:
                for i,p in enumerate(peaks):
                    peak_x = centers[p]

                    ax.axvline(peak_x, linestyle=linestyle, color=linecolor)
                    if peak_number_labels_together:
                        label_y = 0.93 if (i % 2 == 0) else 0.83
                        ax.text(peak_x,
                            label_y,
                            f"{i}",
                            transform=ax.get_xaxis_transform(),
                            verticalalignment='top',
                            horizontalalignment='center',
                            color=linecolor,
                            fontsize=peak_number_label_size)

        for i in (0, 1):
            axs[i].set_xlabel('')
            axs[i].tick_params(labelbottom=True)
        for i in (1, 3):
            axs[i].set_ylabel('')
            axs[i].tick_params(labelleft=True)

        if save_plots:
            output_path = fig_name.with_suffix('.jpeg')
            plt.savefig(str(output_path), dpi=dpi)
            print(f'Saved plot to {output_path}')
        _finish_fig(show_plots)


#---------------- Plot nonlinearity ----------------------------
# xlim and ylim can be 'default', 'none', or tuple(ylim_bottom, ylim_top)
def plot_nonlinearity(peaks_ext,
                      parabola_coeffs, 
                      peak_charge_e_ext, 
                      charge_minus_npeak_ext,
                      fit_range_right_ext,
                      xlim='default', ylim='default',
                      individual_figsize=(6,5), subplots_figsize=(9,7),
                      additional_title='',
                      suptitle='Pixel Charge Nonlinearity',
                      nimages=10,
                      line_color='r', 
                      scatter_color='b', 
                      s=2, 
                      alpha=0.5,
                      plot_individual=False,
                      plot_together=True,
                      sharex=True,
                      sharey=True,
                      show_titles=True,
                      save_plots=False,
                      show_plots=True,
                      fig_path='./', file='nonlinearity_curve',
                      dpi=350):

    fig_path = Path(fig_path)
    if file != 'nonlinearity_curve':
        base_name = file[:-5] + '_nonlinearity'
    else:
        base_name = file
    fig_name = fig_path / base_name

    if type(fit_range_right_ext) is not list:
        fit_range_right_ext = [fit_range_right_ext]*len(peaks_ext)

    if plot_individual:
        for ext, peaks in enumerate(peaks_ext):
            fig, ax = plt.subplots(1, 1, figsize=individual_figsize, constrained_layout=True)
            if show_titles:
                ax.set_title(f'{additional_title}{suptitle} (Nimages = {nimages}): EXT {ext + 1}', fontsize=12, pad=10)
            ax.grid()

            parabola_coeff=parabola_coeffs[ext]
            peak_charge_e=peak_charge_e_ext[ext]
            charge_minus_npeak=charge_minus_npeak_ext[ext]
            fit_range_right=fit_range_right_ext[ext]
            xlim_e = xlim[ext] if isinstance(xlim, list) else xlim
            ylim_e = ylim[ext] if isinstance(ylim, list) else ylim

            ax.plot(peak_charge_e, parabola(peak_charge_e, *parabola_coeff), color=line_color,
                        label=r'$%5.6f x^2 + %5.3f x + %5.3f$' %tuple(parabola_coeff))
            ax.scatter(peak_charge_e, charge_minus_npeak, c=scatter_color, s=s, alpha=alpha)
            ax.legend(loc="upper right", fontsize=10)
            ax.set_xlabel(r'Measured Pixel Charge ($e^-$)')
            ax.set_ylabel(r'Measured Pixel Charge - Peak n. ($e^-$) ')

            if ylim_e=='default':
                fit_idx = int(np.searchsorted(peak_charge_e, fit_range_right))
                ax.set_ylim(min(charge_minus_npeak)-10, max(charge_minus_npeak[:fit_idx])+5)
            elif ylim_e!='none':
                ax.set_ylim(ylim_e)

            if xlim_e=='default':
                ax.set_xlim(-100, fit_range_right+500)
            elif xlim_e!='none':
                ax.set_xlim(xlim_e)

            if save_plots:
                output_path = fig_name.with_stem(fig_name.stem + f'_EXT{ext+1}').with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            _finish_fig(show_plots)

    if plot_together:
        fig, axs = plt.subplots(2, 2, figsize=subplots_figsize, constrained_layout=True, sharex=sharex, sharey=sharey)
        axs=axs.flatten()
        if show_titles:
            fig.suptitle(f'{additional_title}{suptitle} (Nimages = {nimages})')
        for ext, peak_charge_e in enumerate(peak_charge_e_ext):
            ax = axs[ext]
            ax.grid()

            parabola_coeff=parabola_coeffs[ext]
            charge_minus_npeak=charge_minus_npeak_ext[ext]
            fit_range_right=fit_range_right_ext[ext]
            xlim_e = xlim[ext] if isinstance(xlim, list) else xlim
            ylim_e = ylim[ext] if isinstance(ylim, list) else ylim

            ax.plot(peak_charge_e, parabola(peak_charge_e, *parabola_coeff), color=line_color,
                        label=r'$%5.6f x^2 + %5.3f x + %5.3f$' %tuple(parabola_coeff))
            ax.scatter(peak_charge_e, charge_minus_npeak, c=scatter_color, s=s, alpha=alpha)
            ax.legend(loc="upper right", fontsize=8)
            ax.set_title(f'EXT {ext + 1}')
            ax.set_xlabel(r'Measured Pixel Charge ($e^-$)')
            ax.set_ylabel(r'Measured Pixel Charge - Peak n. ($e^-$) ')

            if ylim_e=='default':
                fit_idx = int(np.searchsorted(peak_charge_e, fit_range_right))
                ax.set_ylim(min(charge_minus_npeak)-10, max(charge_minus_npeak[:fit_idx])+5)
            elif ylim_e!='none':
                ax.set_ylim(ylim_e)

            if xlim_e=='default':
                ax.set_xlim(-100, fit_range_right+500)
            elif xlim_e!='none':
                ax.set_xlim(xlim_e)

        for i in (0, 1):
            axs[i].set_xlabel('')
            axs[i].tick_params(labelbottom=True)
        for i in (1, 3):
            axs[i].set_ylabel('')
            axs[i].tick_params(labelleft=True)

        if save_plots:
            output_path = fig_name.with_suffix('.jpeg')
            plt.savefig(str(output_path), dpi=dpi)
            print(f'Saved plot to {output_path}')
        _finish_fig(show_plots)


#---------------- Plot single-electron resolution windows ----------------------------
_RESOLUTION_VERDICT_COLORS = {
    'well resolved': 'tab:green', 'marginal': 'tab:orange', 'unresolved': 'tab:red',
}


def _draw_resolution_window(ax, res, sigma_well=0.25, sigma_limit=0.5):
    """Draw one resolution window (data bars, comb components, comb sum, null,
    and a verdict annotation) onto an existing Axes."""
    xw, yw = res['xw'], res['yw']
    bw = res['bin_width']
    comb, _ = _make_comb_model(res['means'])
    xf = np.linspace(res['lo'], res['hi'], 2000)

    ax.bar(xw, yw, width=bw, align='center', color='0.8', edgecolor='none',
           label='data', zorder=1)

    sigma = res['comb_popt'][0]
    amps = res['comb_popt'][1:]
    for i, (amp, mu) in enumerate(zip(amps, res['means'])):
        ax.plot(xf, _gauss_single(xf, amp, mu, sigma), color='tab:blue', lw=0.8,
                alpha=0.7, zorder=2, label='components' if i == 0 else None)
    ax.plot(xf, comb(xf, *res['comb_popt']), color='tab:red', lw=2.0, zorder=4,
            label='n-Gaussian fit')
    if res['null_popt'] is not None:
        ax.plot(xf, _gauss_single(xf, *res['null_popt']), color='k', lw=1.3,
                ls='--', zorder=3, label='single-Gaussian null')

    verdict = res.get('verdict') or classify_resolution(res, sigma_well, sigma_limit)
    color = _RESOLUTION_VERDICT_COLORS.get(verdict, 'k')
    ax.set_xlabel(r'Charge ($e^-$)')
    ax.set_ylabel('N')
    txt = (rf'$\sigma$ = {res["sigma_e"]:.3f} $\pm$ {res["sigma_e_err"]:.3f} $e^-$ ($\sigma_0$ = {res["sigma_seed_e"]:.3f})'+'\n'
           f'peaks fit: {res["n_components"]} of {res["expected_peaks"]} expected\n'
           rf'reduced $\chi^2$ = {res["reduced_chi2"]:.2f}'
           f'   $\\Delta$AIC = {res["delta_aic"]:.0f}\n'
           f'verdict: {verdict.upper()}'
           rf'  ($\sigma$<{sigma_well} resolved), <{sigma_limit} marginal')
    ax.text(0.02, 0.97, txt, transform=ax.transAxes, va='top', ha='left',
            fontsize=9, zorder=11,
            bbox=dict(boxstyle='round', fc='white', ec=color, lw=1.5, alpha=0.85))
    # Draw the legend last and force it above all data bars, fit curves, and the
    # annotation box so it is never occluded.
    leg = ax.legend(loc='upper right', fontsize=8)
    leg.set_zorder(12)
    return verdict


def plot_resolution(results_ext, charges,
                    individual_figsize=(8, 6.5), subplots_figsize=(13, 9),
                    additional_title='',
                    suptitle='Single-Electron Resolution',
                    nimages=10,
                    sigma_well=0.25, sigma_limit=0.5,
                    plot_individual=False, plot_together=True,
                    sharex=False, sharey=False,
                    show_titles=True,
                    save_plots=False, show_plots=True,
                    fig_path='./', file='resolution', dpi=350):
    """Plot the resolution windows from resolution_at_charge_ext.

    ``plot_individual`` makes one figure per (extension, charge). ``plot_together``
    makes one 2x2 figure per charge (one subplot per extension).
    """
    charges = _normalize_charges(charges)
    fig_path = Path(fig_path)
    if file != 'resolution':
        base_name = file[:-5] + '_resolution'
    else:
        base_name = file
    fig_name = fig_path / base_name
    n_ext = len(results_ext)

    if plot_individual:
        for rows in results_ext:
            for res in rows:
                if res is None:
                    continue
                fig, ax = plt.subplots(figsize=individual_figsize,
                                       constrained_layout=True)
                _draw_resolution_window(ax, res, sigma_well, sigma_limit)
                if show_titles:
                    ax.set_title(
                        f'{additional_title}{suptitle} (Nimages = {nimages}): '
                        f'EXT {res["ext"]}, q = {res["q"]:.0f} $e^-$  '
                        f'(window {res["lo"]:.0f}-{res["hi"]:.0f}, '
                        f'{res["n_components"]}/{res["expected_peaks"]} peaks)',
                        fontsize=11, pad=10)
                if save_plots:
                    output_path = fig_name.with_stem(
                        fig_name.stem + f'_EXT{res["ext"]}_q{res["q"]:.0f}'
                    ).with_suffix('.jpeg')
                    plt.savefig(str(output_path), dpi=dpi)
                    print(f'Saved plot to {output_path}')
                _finish_fig(show_plots)

    if plot_together:
        for j, q in enumerate(charges):
            present = [results_ext[e][j] for e in range(n_ext)
                       if j < len(results_ext[e]) and results_ext[e][j] is not None]
            if not present:
                continue
            fig, axs = plt.subplots(2, 2, figsize=subplots_figsize,
                                    constrained_layout=True, sharex=sharex,
                                    sharey=sharey)
            axs = axs.flatten()
            if show_titles:
                fig.suptitle(f'{additional_title}{suptitle} at q = {q:g} $e^-$ '
                             f'(Nimages = {nimages})')
            for e in range(len(axs)):
                res = results_ext[e][j] if e < n_ext and j < len(results_ext[e]) else None
                if res is None:
                    axs[e].set_visible(False)
                    continue
                _draw_resolution_window(axs[e], res, sigma_well, sigma_limit)
                if show_titles:
                    axs[e].set_title(f'EXT {res["ext"]}')
            if save_plots:
                output_path = fig_name.with_stem(
                    fig_name.stem + f'_q{q:.0f}'
                ).with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            _finish_fig(show_plots)


#---------------- UTILITY FUNCTIONS ----------------------------


        

def get_all_peaks_ext(data_ext, widths, buffers, pedestals, double_gauss_popts, gains, bins='default', flatten=True, do_convert_to_electrons=True, range_left='default', range_right=2000, bin_factor=8, prominences=None, print_values=False):
    counts_ext = []
    edges_ext = []
    peaks_ext = []
    centers_ext = []
    hist_ranges = []

    if type(widths) is not list:
        widths = [widths] * len(data_ext)
    if type(buffers) is not list:
        buffers = [buffers] * len(data_ext)
    if prominences is None or not isinstance(prominences, list):
        prominences = [prominences] * len(data_ext)

    if print_values:
        print('\nFinding peaks for each extension with the following parameters:\n')
        print(f'Widths: {widths}')
        print(f'Buffers: {buffers}')
        print(f'Prominences: {prominences}')
        print(f'Pedestals: {pedestals}')
        print(f'Gains: {gains}')

    for ext, data in enumerate(data_ext):
        width = widths[ext]
        buffer = buffers[ext]
        prominence = prominences[ext]
        pedestal = pedestals[ext]
        noise = double_gauss_popts[ext][0]
        gain = gains[ext]

        counts, edges, peaks, centers, properties, hist_range = find_all_peaks(data, 
                                                                               width, 
                                                                               buffer, 
                                                                               pedestal,
                                                                               noise,
                                                                               gain,
                                                                               bins=bins,
                                                                               flatten=flatten,
                                                                               do_convert_to_electrons=do_convert_to_electrons,
                                                                               range_left=range_left,
                                                                               range_right=range_right,
                                                                               bin_factor=bin_factor,
                                                                               prominence=prominence)
    
        counts_ext.append(counts)
        edges_ext.append(edges)
        peaks_ext.append(peaks)
        centers_ext.append(centers)
        hist_ranges.append(hist_range)
    if print_values:
        print('***********************************************************')
    
    return counts_ext, edges_ext, peaks_ext, centers_ext, hist_ranges

def get_nonlinearity_ext(peaks_ext, centers_ext, pedestals, gains, fit_range_right_ext, do_convert_to_electrons=False, fit_bounds_low=-100, fit_bounds_high=100):
    
    peak_charge_e_ext = []
    charge_minus_npeak_ext = []
    parabola_coeffs = []
    parabola_pcovs = []

    if type(fit_range_right_ext) is not list:
        fit_range_right_ext = [fit_range_right_ext]*len(peaks_ext)

    for ext, peaks in enumerate(peaks_ext):

        centers=centers_ext[ext]
        pedestal=pedestals[ext]
        gain=gains[ext]
        fit_range_right=fit_range_right_ext[ext]

        try:
            parabola_coeff, parabola_pcov, peak_charge_e, charge_minus_npeak = fit_nonlinearity(peaks,
                                                                             centers,
                                                                             pedestal,
                                                                             gain,
                                                                             fit_range_right,
                                                                             do_convert_to_electrons,
                                                                             fit_bounds_low,
                                                                             fit_bounds_high)
        except ValueError as exc:
            raise ValueError(f'EXT {ext + 1}: {exc}') from exc

        peak_charge_e_ext.append(peak_charge_e)
        charge_minus_npeak_ext.append(charge_minus_npeak)
        parabola_coeffs.append(parabola_coeff)
        parabola_pcovs.append(parabola_pcov)

    return peak_charge_e_ext, charge_minus_npeak_ext, parabola_coeffs, parabola_pcovs

def get_nonlinearity_at_ext(q, parabola_coeffs, parabola_pcovs, fit_range_right_ext, save_path=None):
    out_lines = [] if save_path is not None else None

    nonlinearity_at_q_ext = []
    for ext, parabola_coeff in enumerate(parabola_coeffs):
        print('\n***********************************************************')
        header = f'EXT {ext+1}:'
        print(f'\n{header}\n')
        if out_lines is not None:
            out_lines.append('')
            out_lines.append(header)
            out_lines.append('')
        parabola_pcov = parabola_pcovs[ext]
        fit_range_right = fit_range_right_ext[ext]
        nonlinearity_at_q = get_nonlinearity_at(q, parabola_coeff, parabola_pcov, fit_range_right, out_lines=out_lines)
        nonlinearity_at_q_ext.append(nonlinearity_at_q)

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        with save_path.open('w', encoding='utf-8') as f:
            f.write('\n'.join(out_lines))
        print(f'Saved nonlinearity-at-charge text output to {save_path}')

    return nonlinearity_at_q_ext

def _normalize_charges(nonlinearity_charges):
    """Return nonlinearity_charges as a flat list of floats (accepts a scalar or a list)."""
    if nonlinearity_charges is None:
        return []
    if isinstance(nonlinearity_charges, (list, tuple, np.ndarray)):
        return [float(q) for q in nonlinearity_charges]
    return [float(nonlinearity_charges)]

def _nonlin_column_name(charge):
    """CSV/column name for the nonlinearity evaluated at a charge, e.g. 'nonlinearity_at_500_e'."""
    return f'nonlinearity_at_{charge:g}_e'

def build_extension_summary(gains, double_gauss_popts, parabola_coeffs, nonlinearity_charges=500):
    """Build a per-extension summary of gain, noise in electrons, and nonlinearity at charge(s).

    Args:
        gains: list of per-extension gains (ADU/e-), i.e. m1 - m0 from the double-Gaussian fit.
        double_gauss_popts: list of per-extension double-Gaussian popts (s0, m0, s1, m1, N0, N1).
            s0 (popt[0]) is the standard deviation of the zero-electron peak in ADU.
        parabola_coeffs: list of per-extension parabola coefficients (a, b, c) from the
            nonlinearity fit.
        nonlinearity_charges: charge or list of charges (in e-) at which to evaluate the
            nonlinearity (default 500).

    Returns:
        (rows, charges) where charges is the normalized list of charges, and rows is a list of
        dicts (one per extension) with flat keys:
          'ext', 'gain_adu_per_e', 'noise_e', and one 'nonlinearity_at_<charge>_e' per charge.
    """
    charges = _normalize_charges(nonlinearity_charges)
    rows = []
    for ext, gain in enumerate(gains):
        noise_adu = double_gauss_popts[ext][0]  # s0: std of zero-electron peak in ADU
        noise_e = noise_adu / gain  # convert ADU -> e- by dividing by gain (ADU/e-)
        a, b, c = parabola_coeffs[ext]
        row = {
            'ext': ext + 1,
            'gain_adu_per_e': gain,
            'noise_e': noise_e,
        }
        for q in charges:
            row[_nonlin_column_name(q)] = a * q**2 + b * q + c
        rows.append(row)
    return rows, charges

def _summary_fieldnames(charges):
    return ['ext', 'gain_adu_per_e', 'noise_e'] + [_nonlin_column_name(q) for q in charges]

def format_extension_summary(rows, charges):
    """Format per-extension summary rows (from build_extension_summary) as an aligned table string."""
    headers = ['EXT', 'Gain [ADU/e-]', 'Noise [e-]'] + [f'Nonlin @ {q:g} e- [e-]' for q in charges]
    fields = _summary_fieldnames(charges)
    cells = [[f'{r[f]:.4f}' if f != 'ext' else f'{r[f]}' for f in fields] for r in rows]
    widths = [max(len(h), *(len(row[i]) for row in cells)) if cells else len(h)
              for i, h in enumerate(headers)]
    fmt = lambda vals: '  '.join(v.rjust(widths[i]) for i, v in enumerate(vals))
    lines = [fmt(headers), '-' * (sum(widths) + 2 * (len(widths) - 1))]
    lines += [fmt(row) for row in cells]
    return '\n'.join(lines)

def summarize_extensions(gains, double_gauss_popts, parabola_coeffs, nonlinearity_charges=500,
                         save_path=None):
    """Print a per-extension table of gain, noise in e-, and nonlinearity, and optionally save it as CSV.

    Returns the list of summary rows from build_extension_summary.
    """
    rows, charges = build_extension_summary(gains, double_gauss_popts, parabola_coeffs, nonlinearity_charges)
    print('\n***********************************************************')
    print('\nPer-extension summary:\n')
    print(format_extension_summary(rows, charges))

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = _summary_fieldnames(charges)
        with save_path.open('w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f'\nSaved per-extension summary table (CSV) to {save_path}')

    return rows


def parabola(x, a, b, c):
    return a*x**2 + b*x + c
