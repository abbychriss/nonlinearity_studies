"""
fit_range — split out of nonlinearity_studies.py.
"""
import numpy as np

from .fit_data import fit_nonlinearity

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
