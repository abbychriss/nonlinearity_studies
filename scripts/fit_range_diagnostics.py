#!/usr/bin/env python3
"""Diagnostic harness for comparing fit_range_right estimators.

Reconstructs the per-extension nonlinearity curve (peak_charge_e vs
charge_minus_npeak) for one or more FITS files, then computes the right fit
bound several different ways and overlays them on the data so we can judge,
per extension, which estimator tracks the true "clean parabola -> noise" break.

Methods compared
----------------
  ground_truth      : hand-tuned value passed via --truth (one per extension)
  var_a             : current estimator, argmin of parabola pcov[0,0]
  var_a_tol         : var_a but smallest candidate within (1+tol)*min_score
  var_a_pct         : var_a capped at a charge percentile
  var_a_minpk       : var_a rejecting candidates with too few enclosed peaks
  noise_onset       : current experimental rolling-MAD onset detector
  anchored_cp       : PROPOSED anchored-changepoint with persistence

Nothing here is written back into the package; it only reads it. Plots are
saved to a diagnostics subfolder and nothing is ever cleared.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from nonlinearity_studies.nonlinearity_studies import (  # noqa: E402
    get_fits,
    pedestal_subtract_ext_cached,
    get_zero_one_peaks_ext,
    get_all_peaks_ext,
    estimate_optimal_fit_range_right,
    estimate_fit_range_right_by_noise_onset,
    parabola,
)
from scipy.optimize import curve_fit  # noqa: E402


# ----------------------------------------------------------------------------
# Proposed estimator: anchored forward changepoint with persistence
# ----------------------------------------------------------------------------
def _local_roughness(x, y, win):
    """Centered local robust scatter: MAD*1.4826 of residuals from a local
    quadratic fit in a window of `win` points. Removes curvature, leaves noise,
    and (being local) never accumulates global-extrapolation drift."""
    n = len(x)
    half = win // 2
    rough = np.full(n, np.nan)
    for i in range(half, n - half):
        sl = slice(i - half, i + half + 1)
        c = np.polyfit(x[sl], y[sl], 2)
        r = y[sl] - np.polyval(c, x[sl])
        rough[i] = 1.4826 * np.median(np.abs(r - np.median(r)))
    return rough


def anchored_changepoint(peak_charge_e, charge_minus_npeak,
                         win=25, baseline_pct=25, factor=4.0,
                         k=4.0, persist=4, floor=0.15, return_debug=False):
    """Find the charge where the clean parabola gives way to noise.

    Two-stage, drift-free:

    Stage 1 (coarse, robust) -- centered LOCAL roughness (see _local_roughness)
        is flat through the whole clean parabola (curvature removed) and steps
        up at the noise onset. Establish a baseline from the early clean region
        (charges below ``baseline_pct`` percentile) and find the first center
        index where roughness sustains above ``factor * baseline`` for
        ``persist`` consecutive points. This says "noise is near here" without
        being fooled by curvature or extrapolation.

    Stage 2 (precise onset) -- everything up to (coarse_center - half_window) is
        guaranteed clean (the roughness window had not yet touched noise). Fit a
        parabola to those clean points, estimate sigma_clean (floored), then walk
        forward and cut at the first run of ``persist`` points whose residual
        exceeds ``k * sigma_clean``. The refinement span is short, so there is no
        half-window underestimate and no drift.
    """
    x = np.asarray(peak_charge_e, dtype=float)
    y = np.asarray(charge_minus_npeak, dtype=float)
    n = len(x)
    half = win // 2
    if n < win + persist + 4:
        raise ValueError(f'Need at least {win + persist + 4} peaks, got {n}')

    rough = _local_roughness(x, y, win)

    qcut = np.percentile(x, baseline_pct)
    base_mask = (x <= qcut) & np.isfinite(rough)
    if base_mask.sum() < 3:
        base_mask = np.isfinite(rough)
    baseline = float(np.median(rough[base_mask]))
    if not np.isfinite(baseline) or baseline <= 0:
        baseline = float(np.nanmedian(rough)) or 1e-9
    # Absolute floor so an ultra-quiet early region cannot set an impossibly
    # tight threshold that trips on sub-noise waviness. `floor` is in e- and is
    # the minimum roughness that counts as "noise"; it sits between the clean
    # scatter (~0.02-0.05 e-) and the noisy-region roughness (~0.3-0.5 e-).
    thr = max(factor * baseline, float(floor))

    # ---- Stage 1: coarse crossing with persistence ----
    start = int(np.searchsorted(x, qcut))
    start = max(start, half)
    coarse = None
    for i in range(start, n - half - persist):
        seg = rough[i:i + persist]
        if np.all(np.isfinite(seg)) and np.all(seg > thr):
            coarse = i
            break
    if coarse is None:
        if return_debug:
            return float(x[-1]), dict(rough=rough, baseline=baseline, thr=thr,
                                      coarse_idx=None, clean_end=None,
                                      sigma_clean=None, break_idx=n - 1)
        return float(x[-1])

    # ---- Stage 2: precise onset refinement on guaranteed-clean points ----
    clean_end = max(coarse - half, 8)
    coef = np.polyfit(x[:clean_end], y[:clean_end], 2)
    resid = y[:clean_end] - np.polyval(coef, x[:clean_end])
    sigma_clean = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    if not np.isfinite(sigma_clean) or sigma_clean <= 0:
        sigma_clean = float(np.std(resid)) or 1e-9

    # Absolute residual threshold with the same floor logic as stage 1.
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


def var_a_score_curve(peak_charge_e, charge_minus_npeak, min_charge=100, n=80):
    """Return (candidates, pcov[0,0]) for plotting the var(a) score landscape."""
    x = np.asarray(peak_charge_e, dtype=float)
    y = np.asarray(charge_minus_npeak, dtype=float)
    lo = max(min_charge, x[3])
    cands = np.linspace(lo, x[-1], n)
    scores = []
    out = []
    for c in cands:
        fi = int(np.searchsorted(x, c))
        if fi < 4:
            continue
        try:
            _, pcov = curve_fit(parabola, x[:fi], y[:fi], maxfev=2000, bounds=(-100, 100))
            s = pcov[0, 0]
        except Exception:
            continue
        if np.isfinite(s) and s > 0:
            out.append(c)
            scores.append(s)
    return np.array(out), np.array(scores)


# ----------------------------------------------------------------------------
# Data reconstruction (mirrors run_nonlinearity_studies up through peak finding)
# ----------------------------------------------------------------------------
def build_curves(fits_path, do_pedsub=None):
    # Auto-skip pedestal subtraction for files that are already subtracted
    # (named *.pedsub.fits) to avoid double-subtracting.
    if do_pedsub is None:
        do_pedsub = not str(fits_path).endswith('.pedsub.fits')
    data_ext = get_fits(str(fits_path))
    if do_pedsub:
        data_ext = pedestal_subtract_ext_cached(
            data_ext, source_path=Path(fits_path),
            n_std_to_mask=1.5, axis='row',
            use_biweight_loc=True, use_biweight_midvar=True,
            cache_dir=None, force=False, verbose=True,
        )
    else:
        print(f'Skipping pedestal subtraction (already subtracted): {fits_path}')
    _, _, pedestals, gains, double_gauss_popts, _ = get_zero_one_peaks_ext(
        data_ext, n=100, fit_bounds='default',
    )
    _, _, peaks_ext, centers_ext, _ = get_all_peaks_ext(
        data_ext, widths=[0.1, 0.1, 0.1, 0.1], buffers=[3, 3, 3, 3],
        prominences=[None, None, None, None],
        pedestals=pedestals, double_gauss_popts=double_gauss_popts, gains=gains,
        bins='default', flatten=True, do_convert_to_electrons=True,
        range_left='default', range_right=1500, bin_factor=10, print_values=False,
    )
    curves = []
    for ext in range(len(peaks_ext)):
        peaks, centers = peaks_ext[ext], centers_ext[ext]
        x = np.array([centers[p] for p in peaks], dtype=float)
        y = np.array([x[i] - i for i in range(len(peaks))], dtype=float)
        curves.append((peaks, centers, x, y, pedestals[ext], gains[ext]))
    return curves


METHODS = ['ground_truth', 'var_a', 'var_a_tol', 'var_a_pct', 'var_a_minpk',
           'noise_onset', 'anchored_cp']
COLORS = {
    'ground_truth': 'k', 'var_a': 'tab:red', 'var_a_tol': 'tab:orange',
    'var_a_pct': 'tab:brown', 'var_a_minpk': 'tab:pink',
    'noise_onset': 'tab:green', 'anchored_cp': 'tab:blue',
}


def compute_methods(peaks, centers, x, y, pedestal, gain, truth=None):
    res = {}
    if truth is not None:
        res['ground_truth'] = float(truth)

    def safe(fn):
        try:
            return float(fn())
        except Exception as e:
            return ('ERR', str(e)[:60])

    res['var_a'] = safe(lambda: estimate_optimal_fit_range_right(
        peaks, centers, pedestal, gain, do_convert_to_electrons=False))
    res['var_a_tol'] = safe(lambda: estimate_optimal_fit_range_right(
        peaks, centers, pedestal, gain, do_convert_to_electrons=False, tolerance=0.10))
    res['var_a_pct'] = safe(lambda: estimate_optimal_fit_range_right(
        peaks, centers, pedestal, gain, do_convert_to_electrons=False,
        max_charge_percentile=90))
    res['var_a_minpk'] = safe(lambda: estimate_optimal_fit_range_right(
        peaks, centers, pedestal, gain, do_convert_to_electrons=False,
        min_peaks_in_fit=300))
    res['noise_onset'] = safe(lambda: estimate_fit_range_right_by_noise_onset(
        peaks, centers, pedestal, gain, do_convert_to_electrons=False,
        window=20, factor=1.8))
    res['anchored_cp'] = safe(lambda: anchored_changepoint(x, y))
    return res


def plot_file(fits_path, truths, out_dir, label):
    curves = build_curves(fits_path)
    n_ext = len(curves)
    summary = {}

    fig, axes = plt.subplots(n_ext, 2, figsize=(15, 4 * n_ext))
    if n_ext == 1:
        axes = axes[None, :]

    for ext, (peaks, centers, x, y, ped, gain) in enumerate(curves):
        truth = truths[ext] if truths and ext < len(truths) else None
        res = compute_methods(peaks, centers, x, y, ped, gain, truth)
        summary[ext] = res

        # ---- left: nonlinearity curve + breakpoints ----
        axL = axes[ext, 0]
        axL.scatter(x, y, s=3, color='b', alpha=0.6, label='data')
        # parabola fit to the anchored_cp result, for visual sanity
        cp = res.get('anchored_cp')
        if isinstance(cp, float):
            fi = int(np.searchsorted(x, cp))
            if fi >= 4:
                coef, _ = curve_fit(parabola, x[:fi], y[:fi], maxfev=2000, bounds=(-100, 100))
                xx = np.linspace(x[0], x[-1], 400)
                axL.plot(xx, parabola(xx, *coef), color='tab:blue', lw=1.2,
                         label='parabola @ anchored_cp')
        for mname in METHODS:
            v = res.get(mname)
            if isinstance(v, float):
                axL.axvline(v, color=COLORS[mname], ls='--', lw=1.4,
                            label=f'{mname}={v:.0f}')
        axL.set_title(f'{label} EXT {ext+1}: nonlinearity + breakpoints')
        axL.set_xlabel('charge (e-)'); axL.set_ylabel('charge - peak#')
        axL.set_ylim(min(-50, np.percentile(y, 1)), max(10, np.percentile(y, 99)))
        axL.legend(fontsize=7, ncol=2, loc='lower left')
        axL.grid(alpha=0.3)

        # ---- right: diagnostics (var_a score curve + anchored z trace) ----
        axR = axes[ext, 1]
        cands, scores = var_a_score_curve(x, y)
        if len(scores):
            axR.semilogy(cands, scores, color='tab:red', lw=1, label='var(a)=pcov[0,0]')
            axR.set_ylabel('var(a) [log]', color='tab:red')
            axR.tick_params(axis='y', labelcolor='tab:red')
        axR2 = axR.twinx()
        try:
            bp, dbg = anchored_changepoint(x, y, return_debug=True)
            axR2.plot(x, dbg['rough'], color='tab:blue', lw=0.9, alpha=0.8,
                      label='local roughness')
            axR2.axhline(dbg['thr'], color='tab:blue', ls=':', lw=1,
                         label=f"thr={dbg['thr']:.2f}")
            axR2.axhline(dbg['baseline'], color='tab:cyan', ls=':', lw=1,
                         label=f"baseline={dbg['baseline']:.2f}")
            axR2.axvline(bp, color='tab:blue', ls='--', lw=1.4)
            if dbg.get('clean_end') is not None:
                axR2.axvline(x[dbg['clean_end']], color='gray', ls='-', lw=1, alpha=0.5)
            finite = dbg['rough'][np.isfinite(dbg['rough'])]
            top = np.percentile(finite, 98) if len(finite) else 1
            axR2.set_ylim(0, max(top, dbg['thr'] * 1.5))
            axR2.set_ylabel('local roughness (e-)', color='tab:blue')
            axR2.tick_params(axis='y', labelcolor='tab:blue')
            axR2.legend(fontsize=6, loc='upper left')
        except Exception as e:
            axR2.text(0.5, 0.5, f'anchored err:\n{e}', transform=axR2.transAxes, fontsize=7)
        if isinstance(truth, (int, float)):
            axR.axvline(truth, color='k', ls='--', lw=1.4)
        axR.set_title(f'{label} EXT {ext+1}: diagnostics')
        axR.set_xlabel('charge (e-)')
        axR.grid(alpha=0.3)

    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f'{label}_{Path(fits_path).stem}_fitrange_diagnostics.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    print(f'\n===== {label} =====')
    hdr = f'{"EXT":>4} ' + ' '.join(f'{m:>12}' for m in METHODS)
    print(hdr)
    for ext in range(n_ext):
        row = f'{ext+1:>4} '
        for m in METHODS:
            v = summary[ext].get(m)
            if isinstance(v, float):
                row += f'{v:>12.0f}'
            elif isinstance(v, tuple):
                row += f'{"ERR":>12}'
            else:
                row += f'{"-":>12}'
        print(row)
    print(f'saved: {out_path}')
    return summary


def sweep(vr5, vr6, param_sets):
    """Print anchored_cp breakpoints for several parameter sets vs truth."""
    truth = [423, 745, 650, 674]
    files = [('VR-5', vr5), ('VR-6', vr6)]
    cache = {}
    for label, path in files:
        cache[label] = build_curves(path)

    for ps in param_sets:
        print(f"\n--- params: {ps} ---")
        for label, _ in files:
            curves = cache[label]
            vals = []
            for (peaks, centers, x, y, ped, gain) in curves:
                try:
                    vals.append(f'{anchored_changepoint(x, y, **ps):.0f}')
                except Exception as e:
                    vals.append(f'ERR({str(e)[:15]})')
            tline = ' '.join(f'{t:>6}' for t in truth)
            vline = ' '.join(f'{v:>6}' for v in vals)
            print(f'  {label}  truth: {tline}   anchored: {vline}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vr5', default='/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies/VR-5/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits')
    ap.add_argument('--vr6', default='/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies/VR-6/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits')
    ap.add_argument('--out', default='/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies/presentation_plots/fit_range_diagnostics')
    ap.add_argument('--sweep', action='store_true', help='Run parameter sweep instead of plotting')
    args = ap.parse_args()

    if args.sweep:
        param_sets = [
            dict(floor=0.0, factor=4.0, persist=4),
            dict(floor=0.10, factor=4.0, persist=4),
            dict(floor=0.15, factor=4.0, persist=4),
            dict(floor=0.20, factor=4.0, persist=4),
            dict(floor=0.25, factor=4.0, persist=4),
            dict(floor=0.15, factor=4.0, persist=6),
            dict(floor=0.20, factor=5.0, persist=5),
        ]
        sweep(args.vr5, args.vr6, param_sets)
        return

    out = Path(args.out)
    base = '/Users/abbychriss/Privitera_335/data/test_chamber/VR_studies'
    # (label, path, truths)  truths=None where no hand-tuned ground truth exists
    files = [
        ('VR-5', args.vr5, [423, 745, 650, 674]),
        ('VR-6', args.vr6, [423, 745, 650, 674]),
        ('VR-4', f'{base}/VR-4/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits', None),
        ('VR-6-Diego', f'{base}/VR_6_Diego_params/stitched-fits/cds_avg_img_L2_250x3500x500_1x1_L2_125_11_stitched.fits', None),
        ('VR-7', f'{base}/VR-7/stitched-fits/cds_avg_img_L2_250x3500x500_1x1_L2_124_12_stitched.fits', None),
        ('VR-8', f'{base}/VR-8/stitched-fits/cds_avg_img_L2_250x3500x500_1x1_L2_125_10_stitched.fits', None),
        ('VR-9', f'{base}/VR-9/stitched-fits/cds_avg_img_L2_250x3500x500_1x1_L2_125_10_stitched.fits', None),
    ]
    for label, path, truths in files:
        plot_file(path, truths=truths, out_dir=out, label=label)


if __name__ == '__main__':
    main()
