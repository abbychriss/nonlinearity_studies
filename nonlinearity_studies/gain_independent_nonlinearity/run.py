#!/usr/bin/env python3
"""
Gain-independent version of the nonlinearity-studies peak-finding pipeline.

This script follows the exact same steps as
``nonlinearity_studies/run_nonlinearity_studies.py`` -- FITS loading, pedestal
subtraction, and the zero/one double-Gaussian fit are called verbatim, with the
same functions and the same arguments. It deviates from the production
pipeline at exactly one point: the all-peaks spectrum is built with
``do_convert_to_electrons=False``, so ``find_all_peaks`` / ``get_all_peaks_ext``
returns and plots the pedestal-subtracted charge spectrum in ADU instead of
converting it to electrons via the zero/one gain. The zero/one fit is still
run (verbatim) because its pedestal and noise are needed to set the left edge
and peak-width filters of the ADU histogram; its gain is not used to rescale
anything here.

Peak finding is then two steps at two binnings. The coarse scipy.signal
step shares the production key ``nonlinearity_peakfinder_bin_factor`` (bins per
electron, converted to bins per ADU with the nominal ``ADU_PER_E`` so the
detection is binned as upstream bins it), and a sub-bin log-parabola refinement
with no upstream counterpart re-bins the same spectrum at its own
``sub_bin_peakfinder_bin_factor`` (bins per ADU) to localize each of those peaks
more precisely. The all-peaks legend reports both factors, each in both units.

Everything downstream of that -- the nonlinearity parabola fit, single-electron
resolution, zero/one-peak plots, and charge-per-column diagnostic -- is out of
scope for this script and is not run. In their place this package computes the
gain-independent observables: the differential gain delta_g(n) and its
cumulative sum, polynomial fits to the ADU charge spectrum q(n), the running
average gain g_fit(i), the direct gain g_dir(n), and a scan over calibration
gains (see the modules next to this file).

Usage:
As executable:
    ./nonlinearity_studies/gain_independent_nonlinearity/run.py [OPTIONS] <file_path>

As module:
    python -m nonlinearity_studies.gain_independent_nonlinearity.run [OPTIONS] <file_path>

Or after `pip install -e .`:
    gain-independent-nonlinearity [OPTIONS] <file_path>
"""

import numpy as np
import argparse
import glob
import json
from datetime import datetime
from pathlib import Path
import sys
import re
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 12  # bump default for axis labels, tick labels, legends

# Nominal ADU per electron. Used ONLY to convert the electron-valued peak-finder
# knobs -- nonlinearity_peakfinder_bin_factor (bins/e-) and peak_finder_widths
# (electrons) -- into the ADU units this module's histogram is in, so those knobs
# mean the same thing here as in the standard analysis. A fixed number rather than
# the per-extension zero/one gains, so the peak finder stays independent of the
# fit; the data itself is never rescaled by it, and the legend's bins/e- figures
# are computed from the real per-extension gains, not from this.
ADU_PER_E = 1.3

# Same plumbing reuse as run_nonlinearity_studies.py: overscan / fit-column / gain-seed
# helpers come from pedestal_subtract, not reimplemented here.
from pedestal_subtract.core import (
    overscan_cols_from_header,
    get_fits_header,
    _scalar_for_extension,
)
from pedestal_subtract.__main__ import (
    OVERSCAN_COLS,
    _overscan_ext_indices,
    _coerce_fit_cols,
    _normalize_fit_cols,
    _resolve_fit_col_masks,
    _coerce_gain_guess,
)
from pedestal_subtract.constants import _ZERO_ONE_N_BINS

# Absolute imports throughout, with the repo root put on sys.path first, so this
# file works identically whether it is imported as part of the installed package
# or executed as a bare script path (where Python would otherwise only know
# about this file's own directory and relative imports would fail).
_REPO_ROOT = str(Path(__file__).resolve().parents[2])
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from analysis_tools import stitch_fits
from nonlinearity_studies.cli_config import (
    _derive_data_path,
    _wildcard_free_base,
    _load_json_config,
    _config_default,
    _window_scale,
    _n_bins,
    _float_or_auto,
    _lim_token,
    _resolved_args_dict,
    _run_hash,
    _parse_lim,
)
from nonlinearity_studies.nonlinearity_studies import (
    get_fits,
    pedestal_subtract_ext_cached,
    get_zero_one_peaks_ext,
    get_all_peaks_ext,
    scale_buffers,
    plot_all_peaks,
)
from nonlinearity_studies.fit_range import (
    estimate_fit_range_right_changepoint_ext,
)
from nonlinearity_studies.gain_independent_nonlinearity.peaks import (
    sub_bin_log_parabola,
    gaussian_refit_failed_peaks,
    find_zero_peak_ext,
)
from nonlinearity_studies.gain_independent_nonlinearity.deltas import (
    get_standard_gain_ext,
    get_gain_indep_nonlinearity_ext,
    fit_delta_g_line_ext,
    get_cumulative_delta_ext,
    fit_alpha_from_S_ext,
)
from nonlinearity_studies.gain_independent_nonlinearity.spectrum_fit import (
    estimate_fit_n_max_ext,
    fit_charge_spectrum_ext,
    fit_charge_spectrum_cubic_ext,
    get_cubic_truncation_ext,
    get_gain_from_fit_ext,
    fit_charge_spectrum_running_ext,
    fit_running_gain_window_ext,
    fit_gain_fit_residual_cumsum_ext,
    get_gain_fit_residual_ext,
)
from nonlinearity_studies.gain_independent_nonlinearity.direct_gain import (
    get_direct_gain_ext,
    get_direct_delta_ext,
    fit_direct_delta_line_ext,
)
from nonlinearity_studies.gain_independent_nonlinearity.gain_scan import (
    scan_gain_nonlinearity_ext,
)
from nonlinearity_studies.gain_independent_nonlinearity import plots
from nonlinearity_studies.gain_independent_nonlinearity.summary import (
    build_gain_indep_summary,
    format_gain_indep_summary,
    write_gain_indep_summary,
)

from nonlinearity_studies.gain_independent_nonlinearity import (
    compare_nonlinearity_methods, 
    plot_nonlinearity_curves)

def main(args=None):
    """
    The main executable function of the script.

    Args:
        args: The Namespace object containing the parsed command-line arguments.
              If None, arguments are parsed from command line.
    """
    if args is None:
        args = init_argparse()

    if args.extra_plot_title and not args.extra_plot_title.endswith((' ', '\n')):
        args.extra_plot_title = f'{args.extra_plot_title}\n'

    file_path = Path(args.file_string)

    # --- file resolution: identical to run_nonlinearity_studies.py -------------
    if not args.stitch_fits and glob.has_magic(args.file_string):
        matches = sorted(glob.glob(args.file_string))
        if not matches:
            print(f'\nError: no file matches the pattern: {args.file_string}\n')
            sys.exit(1)
        if len(matches) > 1:
            print(f'\nError: file_string matched {len(matches)} files; narrow the pattern '
                  'or use --stitch_fits to combine them:\n  ' + '\n  '.join(matches) + '\n')
            sys.exit(1)
        file_path = Path(matches[0])
    elif not args.stitch_fits and not file_path.exists():
        is_bare_name = (not file_path.is_absolute()) and str(file_path.parent) == '.'
        if not is_bare_name:
            print(f'\nError: file not found: {args.file_string}\n')
            sys.exit(1)

        matches = sorted(set(Path('.').rglob(file_path.name)))
        if len(matches) == 1:
            print(f'Note: "{args.file_string}" not found in the current directory; '
                  f'using the only file matching by name: {matches[0]}')
            file_path = matches[0]
        elif not matches:
            print(f'\nError: file not found: {args.file_string}\n')
            sys.exit(1)
        else:
            listed = '\n  '.join(str(m) for m in matches)
            print(f'\nError: "{args.file_string}" not found, and multiple files are named '
                  f'"{file_path.name}":\n  {listed}\n'
                  f'Pass the exact path to the one you want.\n')
            sys.exit(1)

    do_stitch_images = args.stitch_fits
    do_plot_all_peaks = args.plot_all_peaks_individual or args.plot_all_peaks_together
    save_output = args.save_output
    output_dir = args.output_dir
    verbose = args.verbose

    stitched_dir_name = 'stitched-fits'

    if do_stitch_images:
        input_dir, input_pattern = _derive_data_path(args.file_string)

        if stitched_dir_name in input_dir.parts:
            fits_file_path = next(input_dir.glob(input_pattern), None)
            if fits_file_path is None:
                print('\nError: no files found matching the specified stitched FITS pattern.')
                sys.exit(1)
        else:
            out_path = (_wildcard_free_base(args.file_string) / stitched_dir_name).resolve()
            stitched_file = stitch_fits(
                input_dir.parent,
                directory=input_dir.name,
                image=input_pattern,
                out_path=out_path,
                print_header=verbose,
            )
            if stitched_file is None:
                sys.exit(1)
            fits_file_path = Path(stitched_file)
    else:
        fits_file_path = file_path

    if stitched_dir_name in fits_file_path.parts:
        stitched_fits_idx = fits_file_path.parts.index(stitched_dir_name)
        base_path = Path(*fits_file_path.parts[:stitched_fits_idx])
        default_fig_path = base_path / 'plots'
    else:
        default_fig_path = fits_file_path.parent / 'plots'

    if output_dir is not None:
        fig_path = Path(output_dir)
    else:
        fig_path = default_fig_path

    # Each run gets its own subfolder keyed by a hash of the resolved args, same
    # convention as run_nonlinearity_studies.py, but tagged so a gain-independent
    # (ADU) run never collides with a production (electron) run of the same image
    # even if every other option happens to match.
    run_hash = _run_hash(args)
    fig_path = fig_path / f'{fits_file_path.stem}_gain_independent_{run_hash}'

    if save_output:
        fig_path.mkdir(parents=True, exist_ok=True)
        config_snapshot_path = fig_path / 'config.json'
        snapshot = {
            'run_hash': run_hash,
            'saved_at': datetime.now().isoformat(timespec='seconds'),
            'args': _resolved_args_dict(args),
        }
        with config_snapshot_path.open('w', encoding='utf-8') as f:
            json.dump(snapshot, f, indent=2, sort_keys=True, default=str)
        print(f'Run hash: {run_hash}')
        print(f'Output (plots) and config snapshot will be saved to {fig_path}')

    image_name = fits_file_path.name

    fits_path = str(fits_file_path)
    print(f'Analyzing image: {fits_path}\n')

    # --- load data ---------------------------------------------------------
    data_ext = get_fits(fits_path)
    raw_data_ext = data_ext

    # --- fit columns / overscan: identical plumbing to run_nonlinearity_studies.py
    fit_cols_spec = _normalize_fit_cols(args.fit_cols, len(data_ext))
    fit_cols_ext = _resolve_fit_col_masks(fit_cols_spec, raw_data_ext, args.n_std_to_mask)

    overscan_exts = _overscan_ext_indices(args.use_overscan_only, len(data_ext))
    if overscan_exts:
        overscan_range = overscan_cols_from_header(get_fits_header(fits_path))
        if overscan_range is None:
            overscan_range = tuple(args.overscan_cols)
            if verbose:
                print(f'Overscan columns: header keys not found; using configured {overscan_range}')
        elif verbose:
            print(f'Overscan columns from header: {overscan_range} (last {-overscan_range[0]} columns)')
    else:
        overscan_range = tuple(args.overscan_cols)
    overscan_cols = [overscan_range if i in overscan_exts else None
                     for i in range(len(data_ext))]

    # --- pedestal subtraction: identical call to run_nonlinearity_studies.py ---
    # Keeps the data in ADU -- this is where the production pipeline and this
    # script are still exactly the same.
    if args.do_pedestal_subtraction:
        data_ext = pedestal_subtract_ext_cached(
            data_ext,
            source_path=fits_file_path,
            n_std_to_mask=args.n_std_to_mask,
            axis=args.pedestal_subtraction_axis,
            use_biweight_loc=args.use_biweight_loc,
            use_biweight_midvar=args.use_biweight_midvar,
            max_iter=args.pedsub_max_iter,
            cache_dir=args.pedsub_cache_dir,
            force=args.force_pedsub,
            verbose=verbose,
            overscan_cols=overscan_cols,
        )

    # Full pedestal-subtracted frame, kept aside before the fit_cols restriction
    # below rebinds data_ext -- same reasoning as run_nonlinearity_studies.py: the
    # all-peaks finder scans every pixel, not just the zero/one calibration columns.
    full_data_ext = data_ext

    if any(fc is not None for fc in fit_cols_ext):
        data_ext = [
            np.asarray(d)[:, mask] if mask is not None else d
            for d, mask in zip(data_ext, fit_cols_ext)
        ]
        if verbose:
            kept = [int(m.sum()) if m is not None else None for m in fit_cols_ext]
            print(f'Restricting the fit columns per extension: spec={fit_cols_spec}, '
                  f'kept_columns={kept}')

    # Extract number of stitched images from filename; --nimages arg or config overrides.
    nimages = None
    match = re.search(r'_(\d+)_stitched', fits_path)
    if match:
        nimages = int(match.group(1))
    elif stitched_dir_name not in fits_file_path.parts:
        nimages = 1
    if args.nimages is not None:
        nimages = args.nimages

    if verbose and args.zero_one_gain_guess is not None:
        gain_seeds_ext = [_scalar_for_extension(args.zero_one_gain_guess, ext, len(data_ext))
                          for ext in range(len(data_ext))]
        print(f'Using guess for gains: {gain_seeds_ext}')

    # --- zero/one fit: identical call to run_nonlinearity_studies.py -----------
    # Run verbatim, for the same reason production runs it: it is what yields the
    # per-extension pedestal (needed as the zero point / left-edge anchor of the
    # ADU histogram) and noise (needed for the ADU histogram's left edge and the
    # peak-width filter). The gain it returns is NOT used to rescale the spectrum
    # here -- that omission is the one deviation this script makes from the
    # production pipeline.
    zero_one_counts_ext, zero_one_edges_ext, pedestals, gains, double_gauss_popts, zero_one_ranges = get_zero_one_peaks_ext(
        data_ext,
        n_bins=args.zero_one_n_bins,
        fit_bounds='default',
        window_left_scale=args.zero_one_window_left_scale,
        window_right_scale=args.zero_one_window_right_scale,
        gain_seed=args.zero_one_gain_guess,
    )

    # Right edge of the ADU histogram. In the production pipeline this is
    # expressed directly in electrons (default 1500) because the data has
    # already been gain-divided by the time it is histogrammed. Here the data
    # stays in ADU, so the same "how many electrons out do we want to see"
    # intent is converted to ADU per extension using the zero/one gain solely
    # to size the histogram window -- the gain still never rescales the data
    # itself. Any extension whose fit failed (non-finite gain) falls back to
    # --all_peaks_range_right_adu_fallback so a single bad extension can't
    # crash the whole run.
    if args.all_peaks_range_right_adu is not None:
        all_peaks_range_right = float(args.all_peaks_range_right_adu)
    else:
        candidates = [
            pedestals[e] + args.all_peaks_n_max_electrons * gains[e]
            for e in range(len(data_ext))
            if np.isfinite(gains[e])
        ]
        all_peaks_range_right = (
            float(max(candidates)) if candidates else float(args.all_peaks_range_right_adu_fallback)
        )

    # Peak finding here is two steps at two binnings, one key each:
    #
    #   coarse   scipy.signal.find_peaks, keyed on the production pipeline's
    #            nonlinearity_peakfinder_bin_factor -- bins per ELECTRON, so the
    #            detection step is binned exactly as run_nonlinearity_studies.py
    #            bins it. This spectrum is in ADU, so the number is converted to
    #            bins per ADU (below) for the histogram find_peaks actually runs
    #            on; the physical bin width is the same either way.
    #   sub-bin  the log-parabola refinement, keyed on this module's own
    #            sub_bin_peakfinder_bin_factor -- bins per ADU, since it is a
    #            localization step on the ADU spectrum with no counterpart in
    #            the production pipeline.
    coarse_bin_factor_e = args.nonlinearity_peakfinder_bin_factor
    coarse_bin_factor_adu = coarse_bin_factor_e / ADU_PER_E
    sub_bin_factor = args.sub_bin_peakfinder_bin_factor

    # peak_finder_widths are in electrons, and find_all_peaks converts them with
    # whatever bin factor it is handed -- bins per ADU here. Pre-multiplying by
    # ADU_PER_E cancels that, so width * ADU_PER_E * coarse_bin_factor_adu is
    # width * coarse_bin_factor_e: the same number of bins the production
    # pipeline would use, and still an electron-valued knob.
    _widths = args.peak_finder_widths
    peak_finder_widths_adu = ([w * ADU_PER_E for w in _widths]
                              if isinstance(_widths, list) else _widths * ADU_PER_E)

    # peak_finder_buffers is quoted normalized to 10 bins/e-, so scale_buffers
    # gives the same number of bins as in the production pipeline. find_all_peaks
    # then computes distance = bin_factor - buffer, and its bin_factor is bins per
    # ADU (one ADU's worth of bins) where the production formula subtracts from
    # one ELECTRON's worth; the offset below restores that, leaving
    # distance = coarse_bin_factor_e - buffer bins exactly as upstream.
    peak_finder_buffers_bins = [
        coarse_bin_factor_adu - (coarse_bin_factor_e - b)
        for b in scale_buffers(args.peak_finder_buffers, coarse_bin_factor_e, len(data_ext))
    ]

    # --- THE ONE DEVIATION FROM run_nonlinearity_studies.py --------------------
    # do_convert_to_electrons=False: find_all_peaks (called once per extension by
    # get_all_peaks_ext) histograms the pedestal-subtracted ADU data directly,
    # instead of first dividing by the zero/one gain. Everything else about the
    # call -- widths, buffers, prominences, pedestals, double_gauss_popts, gains --
    # is passed exactly as run_nonlinearity_studies.py passes it; the bin factor
    # is the coarse bins-per-electron key converted to bins per ADU (see above).
    # `data_adu_ext` returned here is the pedestal-subtracted charge spectrum in
    # ADU that this script exists to produce.
    counts_ext, edges_ext, peaks_ext, centers_ext, hist_ranges, data_adu_ext = get_all_peaks_ext(
        full_data_ext,
        widths=peak_finder_widths_adu,
        buffers=peak_finder_buffers_bins,
        prominences=args.peak_finder_prominences,
        pedestals=pedestals,
        double_gauss_popts=double_gauss_popts,
        gains=gains,
        bins='default',
        flatten=True,
        do_convert_to_electrons=False,
        range_left='default',
        range_right=all_peaks_range_right,
        nonlinearity_peakfinder_bin_factor=coarse_bin_factor_adu,
        print_values=verbose,
    )

    if verbose:
        for ext in range(len(data_ext)):
            print(f'EXT {ext + 1}: {len(peaks_ext[ext])} peaks found in ADU spectrum '
                  f'(pedestal={pedestals[ext]:.3f} ADU, histogram range={hist_ranges[ext]}, '
                  f'coarse bin factor={coarse_bin_factor_adu:.4g} bins/ADU '
                  f'= {coarse_bin_factor_e:g} bins/e-, '
                  f'buffer={peak_finder_buffers_bins[ext]:.3g} bins)')

    # --- sub-bin peak localization -----------------------------------------
    # Refines each scipy-found peak (a bin index) to a sub-bin location via the
    # log(counts) parabola vertex; see sub_bin_log_parabola's docstring.
    #
    # The refinement runs on its own histogram: the same ADU spectrum over the
    # same range, re-binned at sub_bin_peakfinder_bin_factor (bins/ADU, finer
    # than the coarse detection binning), with each scipy peak re-located to the
    # fine bin holding its ADU position. Only the sub-bin inputs change --
    # counts_ext and peaks_ext stay on the coarse grid, which is what the peak
    # numbering (find_zero_peak_ext) and the all-peaks histogram are built from.
    fine_counts_ext, fine_centers_ext, fine_peaks_ext = [], [], []
    for ext in range(len(data_ext)):
        lo, hi = hist_ranges[ext]
        n_bins = int(np.floor((hi - lo) * sub_bin_factor))
        counts, edges = np.histogram(data_adu_ext[ext], bins=n_bins, range=(lo, hi))
        centers = 0.5 * (edges[1:] + edges[:-1])
        # bins are uniform, so the fine bin holding each scipy peak's ADU
        # position is pure arithmetic -- no search needed
        peak_adu = centers_ext[ext][peaks_ext[ext]]
        idx = np.floor((peak_adu - lo) / ((hi - lo) / n_bins)).astype(int)
        fine_counts_ext.append(counts)
        fine_centers_ext.append(centers)
        fine_peaks_ext.append(np.clip(idx, 0, n_bins - 1))

    sub_bin_ext = [
        sub_bin_log_parabola(
            fine_counts_ext[ext], fine_centers_ext[ext], fine_peaks_ext[ext],
            half_width=args.sub_bin_half_width,
            min_good_bins=args.sub_bin_min_good_bins,
        )
        for ext in range(len(data_ext))
    ]
    if verbose:
        for ext in range(len(data_ext)):
            n_bad = int((~sub_bin_ext[ext]['ok']).sum())
            if n_bad:
                next_step = ('refitting with a seeded Gaussian' if args.gaussian_peakfind_fallback
                             else 'falling back to the raw bin center (Gaussian fallback disabled)')
                print(f'EXT {ext + 1}: {n_bad}/{len(peaks_ext[ext])} peaks failed the '
                      f'log-parabola fit; {next_step}')

    # Fallback for peaks the log-parabola window couldn't localize: a seeded
    # single-Gaussian curve_fit over a window sized by the zero/one sigma. Only
    # the failed subset (typically a few percent) is refit; see
    # gaussian_refit_failed_peaks's docstring. Gated behind a config/CLI flag
    # since curve_fit's per-peak Python loop is meaningfully slower than the
    # batched log-parabola solve, and may not be worth it in practice.
    if args.gaussian_peakfind_fallback:
        for ext in range(len(data_ext)):
            gaussian_refit_failed_peaks(
                fine_counts_ext[ext], fine_centers_ext[ext], fine_peaks_ext[ext], sub_bin_ext[ext],
                sigma_seed=double_gauss_popts[ext][0],
                k_sigma=args.gaussian_refit_k_sigma,
                min_good_bins=args.sub_bin_min_good_bins,
            )

    mu_ext = [sub_bin['mu'] for sub_bin in sub_bin_ext]

    if verbose and args.gaussian_peakfind_fallback:
        for ext in range(len(data_ext)):
            n_bad = int((~sub_bin_ext[ext]['ok']).sum())
            if n_bad:
                print(f'EXT {ext + 1}: {n_bad}/{len(peaks_ext[ext])} peaks still fell back '
                      f'to the raw bin center (Gaussian refit also failed)')

    # --- peak numbering ----------------------------------------------------
    # q_ext[ext][n] is the ADU charge of peak number n, with n = 0 the
    # zero-electron peak; see find_zero_peak_ext's docstring for why the peak
    # finder's own index 0 cannot be trusted to be that peak.
    # g_ext is the standard gain (= q(1) - q(0)) but run on the sub-bin localized peak locations, 
    # which is not necessarily the same as the double guassian fit to the zero-one peaks
    q_ext, zero_idx_ext = find_zero_peak_ext(counts_ext, peaks_ext, mu_ext)
    g_ext = get_standard_gain_ext(q_ext)

    if verbose:
        for ext in range(len(data_ext)):
            if zero_idx_ext[ext]:
                print(f'EXT {ext + 1}: dropped {zero_idx_ext[ext]} peak(s) left of the '
                      f'zero-electron peak before numbering')
            print(f'EXT {ext + 1}: standard gain g = q(1) - q(0) = {g_ext[ext]:.6f} ADU/e-')

    # No explicit --gain_indep_fit_n_max means "read the fit domain off the
    # data": the changepoint detector walks each extension's charge spectrum
    # and stops where the clean parabola gives way to noise, giving one n_max
    # per extension that every fit below then respects.
    n_max = args.gain_indep_fit_n_max
    if n_max is None:
        if verbose:
            print('\nEstimating the fit domain per extension via changepoint detection...')
        n_max = estimate_fit_n_max_ext(q_ext, g_ext, verbose=verbose)
    n_min = args.gain_indep_fit_n_min

    # --- differential gain and everything derived from it --------------------
    gain_indep_nonlin_ext = get_gain_indep_nonlinearity_ext(q_ext, g_ext)
    delta_g_line_fits_ext = fit_delta_g_line_ext(gain_indep_nonlin_ext, n_max=n_max)
    I_ext, S_ext = get_cumulative_delta_ext(gain_indep_nonlin_ext)
    alpha_fits_ext = fit_alpha_from_S_ext(S_ext, n_max=n_max)

    # --- polynomial fits to the ADU charge spectrum q(n) ---------------------
    quad_fits_ext = fit_charge_spectrum_ext(q_ext, n_max=n_max)
    gain_from_fit_ext = get_gain_from_fit_ext(quad_fits_ext, g_ext)
    cubic_fits_ext = None
    cubic_trunc_ext = None
    if args.gain_indep_fit_cubic:
        cubic_fits_ext = fit_charge_spectrum_cubic_ext(q_ext, n_max=n_max)
        cubic_trunc_ext = get_cubic_truncation_ext(quad_fits_ext, cubic_fits_ext, q_ext, n_max=n_max)

    # --- running fit: the average gain over the first i peaks ----------------
    running_ext = fit_charge_spectrum_running_ext(q_ext, n_min=n_min, n_max=n_max)
    gain_residual_ext, gain_residual_cumsum_ext = get_gain_fit_residual_ext(running_ext, g_ext)
    # The stretch of that curve the average-gain plot reports on: its mean and
    # its trend line, over a window of fit ranges set independently of n_min/n_max
    # so narrowing it leaves every other fit in this module alone.
    running_window_ext = fit_running_gain_window_ext(
        running_ext, n_min=args.gain_indep_avg_gain_n_min,
        n_max=args.gain_indep_avg_gain_n_max)
    gain_residual_cumsum_fits_ext = fit_gain_fit_residual_cumsum_ext(
        running_ext, gain_residual_cumsum_ext)

    # --- the gain read directly off the raw peaks ----------------------------
    g_dir_ext = get_direct_gain_ext(q_ext)
    delta_dir_ext, I_dir_ext, S_dir_ext = get_direct_delta_ext(g_dir_ext, g_ext)
    delta_dir_line_fits_ext = fit_direct_delta_line_ext(delta_dir_ext)
    alpha_fits_dir_ext = fit_alpha_from_S_ext(S_dir_ext, n_max=n_max)

    # --- scan the calibration gain and refit the standard nonlinearity -------
    # 'auto' (the default) places the parabola's right edge per extension with the production
    # changepoint detector, run on the zero/one-gain nonlinearity curve: the
    # scan and the nonlinearity-curve plot below then stop where that
    # extension's clean parabola gives way to noise instead of at one number
    # shared by all four. A plain number is used for every extension.
    scan_fit_range_right = args.gain_indep_scan_fit_range_right
    scan_ext = None
    if args.gain_indep_gain_scan:
        if isinstance(scan_fit_range_right, str):
            if verbose:
                print('\nEstimating the gain-scan fit_range_right per extension '
                      'via changepoint detection...')
            scan_fit_range_right, scan_fit_range_diagnostics = estimate_fit_range_right_changepoint_ext(
                peaks_ext, centers_ext, pedestals, gains,
                do_convert_to_electrons=True,
                fix_intercept=args.nonlinearity_fit_fix_intercept, verbose=verbose,
            )
            low = [d['ext'] for d in scan_fit_range_diagnostics if d.get('confidence') == 'LOW']
            if low and verbose:
                print(f'  WARNING: low-confidence fit_range_right on EXT {low} '
                      f'(changepoint and var(a) disagree); inspect the nonlinearity curve plot(s).')
        if verbose:
            print(f'\nScanning {sum(len(run["i"]) for run in running_ext)} candidate gains...')
        scan_ext = scan_gain_nonlinearity_ext(
            q_ext, running_ext,
            fit_range_right=scan_fit_range_right,
            eval_charge=args.gain_indep_scan_eval_charge,
            fix_intercept=args.nonlinearity_fit_fix_intercept,
        )

    # nonlinearity_at_gains_ext is a list of four lists (one per ext) containing the nonlinearity at q for the gains tested
    # gain_range is a list of the gains tested
    # q is the charge that the nonlinearity was evaluated at
    # The sweep is 400 candidate gains x 4 extensions of parabola fits and nothing
    # but its own plot consumes it, so it is skipped when that plot is switched off.
    nonlin_vs_gain_out = None
    if args.gain_indep_gain_scan and (args.plot_nonlin_vs_gain_individual
                                      or args.plot_nonlin_vs_gain_together):
        nonlin_vs_gain_out = compare_nonlinearity_methods.nonlin_at_q_vs_gain(peaks_ext, centers_ext, pedestals, q=500, gain_range=np.linspace(1.0, 1.4, 400), fit_range_right_ext=1000.0, fit_bounds_low=-100, fit_bounds_high=100, fix_intercept=args.nonlinearity_fit_fix_intercept)

    # alpha_dict_ext is list of 4 dicts (one per ext) containing ``alpha_at_max_i`, ``max_i``, ``alpha_at_best_i``, ``best_i``,``alpha_at_test_i`, ``test_i`, ``std_alphas``, ``max_gain_diff``.
    # index_alpha_gfit_sorted_ext is a list of dicts (one per extension) ``alphas_gains_sorted_ext`` containing all the indices, alphas, and gains sampled, sorted by index.
    test_i = args.gain_indep_compare_test_i
    # Both of these read best_i out of the gain scan, so they only run when the
    # scan did; alpha_dict_ext does not exist otherwise.
    alpha_dict_ext = index_alpha_gfit_sorted_ext = index_alpha_prime_ext = None
    if args.gain_indep_gain_scan:
        alpha_dict_ext, index_alpha_gfit_sorted_ext = compare_nonlinearity_methods.alpha_vs_index_ext(fit_ADU_spectrum_coeffs=running_ext, gain_optimization_out=scan_ext, test_i=test_i)

        # List of 4 dicts (one per ext) containing ``i`` and ``alpha``
        index_alpha_prime_ext = compare_nonlinearity_methods.alpha_prime_vs_index_ext(S_ext, alpha_dict_ext, test_i=test_i)

    # --- summary ------------------------------------------------------------
    summary_rows = build_gain_indep_summary(
        q_ext, zero_idx_ext, g_ext, quad_fits_ext, gain_from_fit_ext,
        delta_g_line_fits_ext, alpha_fits_ext, gain_residual_cumsum_ext,
        cubic_trunc_ext=cubic_trunc_ext, scan_ext=scan_ext, n_max=n_max)
    print(format_gain_indep_summary(summary_rows, eval_charge=args.gain_indep_scan_eval_charge))
    if save_output:
        write_gain_indep_summary(summary_rows, fig_path / 'gain_independent_summary.csv')

    # --- plot: same plot_all_peaks function as run_nonlinearity_studies.py -----
    # Same vertical-line/peak-number-label behavior. It is fed the ADU
    # peaks/centers/data computed above instead of the electron-converted ones,
    # so the histogram it draws is the ADU charge spectrum with a line and
    # number at every found peak. plot_all_peaks' xlabel defaults to
    # "Charge (e-)" for every other caller; here it is overridden to
    # "Charge (ADU)" via the xlabel kwarg added for this purpose.
    if do_plot_all_peaks:
        all_peaks_limit_kwargs = {}
        if args.plot_all_peaks_xlim is not None:
            all_peaks_limit_kwargs['xlim'] = _parse_lim(args.plot_all_peaks_xlim, n_ext=len(data_ext))
        if args.plot_all_peaks_ylim is not None:
            all_peaks_limit_kwargs['ylim'] = _parse_lim(args.plot_all_peaks_ylim, n_ext=len(data_ext))

        plot_all_peaks(counts_ext,
                    peaks_ext,
                    centers_ext,
                    **all_peaks_limit_kwargs,
                    yscale=args.plot_all_peaks_yscale,
                    data_ext=data_adu_ext,
                    window_bin_factor=args.plot_all_peaks_window_bin_factor,
                    peakfinder_bin_factor=coarse_bin_factor_adu,
                    sub_bin_bin_factor=sub_bin_factor,
                    # this spectrum is binned in ADU, so the legend says bins/ADU and
                    # adds the bins/e- equivalent (factor * gain, display only) to make
                    # the binning comparable to the standard bins-per-electron plot
                    bin_factor_unit='bins/ADU',
                    bin_factor_gains=gains,
                    plot_individual=args.plot_all_peaks_individual,
                    plot_together=args.plot_all_peaks_together,
                    sharex=args.plot_all_peaks_sharex,
                    sharey=args.plot_all_peaks_sharey,
                    show_titles=args.show_titles,
                    raw_peak_location_vertical_line=args.raw_peak_location_vertical_line,
                    linecolor='r',
                    linestyle='--',
                    raw_peak_number=args.raw_peak_number,
                    sub_bin_peaks_ext=mu_ext,
                    sub_bin_linecolor='tab:orange',
                    sub_bin_linestyle='--',
                    sub_bin_peak_location_vertical_line=args.sub_bin_peak_location_vertical_line,
                    sub_bin_peak_number=args.sub_bin_peak_number,
                    # the fitted curves themselves, one per successfully localized
                    # peak, colored by which method produced it
                    sub_bin_fits_ext=sub_bin_ext,
                    peak_number_labels_size=args.peak_number_labels_size,
                    individual_figsize=tuple(args.plot_all_peaks_individual_figsize),
                    subplots_figsize=tuple(args.plot_all_peaks_subplots_figsize),
                    additional_title=args.extra_plot_title,
                    suptitle='Peaks in Pixel Charge Distribution (ADU, gain-independent)',
                    nimages=nimages,
                    save_plots=save_output,
                    show_plots=args.show_plots,
                    fig_path=str(fig_path),
                    file=image_name,
                    dpi=350,
                    xlabel='Charge (ADU)')

    # --- gain-independent plots ---------------------------------------------
    # All of them share one set of figure-size / shared-axis settings; only the
    # individual/together switches are per plot.
    plot_kwargs = dict(
        individual_figsize=tuple(args.plot_gain_indep_nonlin_individual_figsize),
        subplots_figsize=tuple(args.plot_gain_indep_nonlin_together_figsize),
        sharex=args.plot_gain_indep_nonlin_sharex,
        sharey=args.plot_gain_indep_nonlin_sharey,
        show_titles=args.show_titles,
        additional_title=args.extra_plot_title,
        nimages=nimages,
        save_plots=save_output,
        show_plots=args.show_plots,
        fig_path=str(fig_path),
        file=image_name,
    )

    def _switches(name):
        return dict(plot_individual=getattr(args, f'plot_{name}_individual'),
                    plot_together=getattr(args, f'plot_{name}_together'))

    plots.plot_gain_indep_nonlinearity_ext(
        gain_indep_nonlin_ext, delta_g_line_fits_ext,
        **_switches('gain_indep_nonlin'), **plot_kwargs)
    plots.plot_charge_spectrum_fit_ext(
        q_ext, quad_fits_ext, cubic_fits_ext,
        **_switches('charge_spectrum_fit'), **plot_kwargs)
    plots.plot_charge_spectrum_ext(
        # q(n) spans a different range in every extension, so this one plot
        # overrides the shared sharey setting the rest of them use.
        q_ext, quad_fits_ext, cubic_fits_ext,
        **_switches('charge_spectrum'), **{**plot_kwargs, 'sharey': False})
    plots.plot_cumulative_delta_ext(
        S_ext, alpha_fits_ext,
        **_switches('cumulative_delta'), **plot_kwargs)
    plots.plot_cumulative_delta_ext(
            S_dir_ext, alpha_fits_dir_ext,
            **_switches('cumulative_dir_delta'),
            **{**plot_kwargs, 'suptitle': 'Cumulative Delta Direct Gain',
               'base_suffix': 'cumulative_dir_delta'})
    plots.plot_running_gain_ext(
        running_ext, g_ext, running_window_ext, gain_residual_ext,
        **_switches('running_gain'), **plot_kwargs)
    plots.plot_gain_fit_residual_ext(
        running_ext, gain_residual_ext,
        **_switches('gain_fit_residual'), **plot_kwargs)
    plots.plot_gain_fit_residual_cumsum_ext(
        running_ext, gain_residual_cumsum_ext, gain_residual_cumsum_fits_ext,
        **_switches('gain_fit_residual_cumsum'), **plot_kwargs)
    plots.plot_direct_gain_ext(
        g_dir_ext, g_ext, running_ext,
        **_switches('direct_gain'), **plot_kwargs)
    plots.plot_direct_delta_ext(
        delta_dir_ext, delta_dir_line_fits_ext,
        **_switches('direct_delta'), **plot_kwargs)
    if scan_ext is not None:
        plots.plot_gain_scan_ext(
            scan_ext, eval_charge=args.gain_indep_scan_eval_charge,
            **_switches('gain_scan'), **plot_kwargs)

    # Compare nonlinearity curves calculated from different gains. Needs the
    # gain scan for best_gain / best_i, so it is gated on it like the plots above.
    if scan_ext is not None:
        # One calibration gain per extension, not a curve: g_dir_ext is
        # g_dir(n), a value at every peak, and handing that to
        # get_nonlinearity_ext would divide the spectrum by an array. The scan
        # already computes the direct gain out to the peak eval_charge sits at
        # and treats it as the comparison point for the scanned candidates.
        g_dirs = [scan['g_dir_eval'] for scan in scan_ext]
        print(g_dirs)
        plot_nonlinearity_curves.plot_nonlinearity_curves(zero_one_gains=gains, sub_bin_gains=g_ext,
                                                          fit_ADU_spectrum_coeffs=running_ext,
                                                          gain_optimization_out=scan_ext, test_i=test_i,
                                                          g_dirs=g_dirs, peaks_ext=peaks_ext, centers_ext=centers_ext,
                                                          pedestals=pedestals,
                                                          fit_range_right_ext=scan_fit_range_right,
                                                          fix_intercept=args.nonlinearity_fit_fix_intercept,
                                                          **_switches('compare_nonlinearity_curves'),
                                                          # Each curve is converted to e- with a different
                                                          # gain, so the panels do not span a common y range;
                                                          # same override charge_spectrum uses above.
                                                          **{**plot_kwargs, 'sharey': False})

    # alpha / alpha' vs peak number. Both inputs come out of the gain scan, so
    # they are None without it and there is nothing to draw.
    if index_alpha_gfit_sorted_ext is not None or index_alpha_prime_ext is not None:
        # This one takes only the figure sizes, not the shared save/show/title
        # kwargs the *_ext plots above accept.
        plots.plot_alpha_alpha_prime_vs_gain_and_index_ext(
            index_alpha_gfit_sorted_ext, index_alpha_prime_ext,
            individual_figsize=tuple(args.plot_gain_indep_nonlin_individual_figsize),
            together_figsize=tuple(args.plot_gain_indep_nonlin_together_figsize),
            **_switches('alpha_alpha_prime'))

    if nonlin_vs_gain_out is not None:
        nonlinearity_at_gains_ext, gain_range, q = nonlin_vs_gain_out
        plots.plot_nonlin_vs_gain(
            nonlinearity_at_gains_ext, gain_range, q,
            **_switches('nonlin_vs_gain'), **plot_kwargs)

    return dict(
        pedestals=pedestals,
        gains=gains,
        double_gauss_popts=double_gauss_popts,
        counts_ext=counts_ext,
        edges_ext=edges_ext,
        peaks_ext=peaks_ext,
        centers_ext=centers_ext,
        hist_ranges=hist_ranges,
        data_adu_ext=data_adu_ext,
        sub_bin_ext=sub_bin_ext,
        q_ext=q_ext,
        zero_idx_ext=zero_idx_ext,
        g_ext=g_ext,
        gain_indep_nonlin_ext=gain_indep_nonlin_ext,
        delta_g_line_fits_ext=delta_g_line_fits_ext,
        I_ext=I_ext,
        S_ext=S_ext,
        alpha_fits_ext=alpha_fits_ext,
        quad_fits_ext=quad_fits_ext,
        cubic_fits_ext=cubic_fits_ext,
        cubic_trunc_ext=cubic_trunc_ext,
        gain_from_fit_ext=gain_from_fit_ext,
        running_ext=running_ext,
        gain_residual_ext=gain_residual_ext,
        gain_residual_cumsum_ext=gain_residual_cumsum_ext,
        g_dir_ext=g_dir_ext,
        delta_dir_ext=delta_dir_ext,
        scan_ext=scan_ext,
        summary_rows=summary_rows,
        fig_path=fig_path,
    )


def _add_plot_switches(parser, config, name, help_text, default_individual=False,
                       default_together=True):
    """Add the --plot_<name>_individual / _together pair (and their --no- forms)."""
    parser.add_argument(f"--plot_{name}_individual", action="store_true",
                        default=_config_default(config, f'plot_{name}_individual', default_individual),
                        help=f"Plot {help_text} as one figure per extension")
    parser.add_argument(f"--no-plot_{name}_individual", dest=f"plot_{name}_individual",
                        action="store_false",
                        help=f"Disable the individual {help_text} plots")
    parser.add_argument(f"--plot_{name}_together", action="store_true",
                        default=_config_default(config, f'plot_{name}_together', default_together),
                        help=f"Plot {help_text} as a combined 2x2 subplot")
    parser.add_argument(f"--no-plot_{name}_together", dest=f"plot_{name}_together",
                        action="store_false",
                        help=f"Disable the combined 2x2 {help_text} subplot")


def init_argparse(args=None):
    """
    Initializes the ArgumentParser object and defines arguments.

    Accepts the same JSON-config mechanism and the same subset of config keys
    used by run_nonlinearity_studies.py that are relevant to this script's
    scope (file loading/stitching, pedestal subtraction, overscan/fit_cols, the
    zero/one fit, peak finding, and the plot_all_peaks display options), plus
    the gain-independent calculation and plot options defined here. A config
    file written for run_nonlinearity_studies.py (e.g. containing
    nonlinearity-fit or resolution keys) can be passed in as-is; keys outside
    this scope are simply not read.
    """
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument("-j", "--json", dest="config", type=str, default=None,
                            help="Path to JSON file containing command-line arguments")
    pre_args, _ = pre_parser.parse_known_args(args)

    config = {}
    if pre_args.config is not None:
        try:
            config = _load_json_config(pre_args.config)
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            pre_parser.error(str(exc))

    parser = argparse.ArgumentParser(
        description="""Gain-independent nonlinearity-studies peak finder.

Runs the same FITS loading, pedestal subtraction, and zero/one double-Gaussian
fit as run_nonlinearity_studies.py, then finds every electron peak directly in
the pedestal-subtracted ADU spectrum instead of converting to electrons via the
zero/one gain, and computes the gain-independent observables from it.""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-j", "--json", type=str, default=pre_args.config,
                        help="Path to JSON file containing command-line arguments. Explicit CLI arguments override JSON values.")
    parser.add_argument('file_string', type=str, nargs='?',
                       default=_config_default(config, 'file_string', None),
                       help='Absolute or relative path to image file (.fz or .fits accepted) if not stitching, or to directory with images if stitching (can include glob pattern like "*"). May also be set in JSON config.')
    parser.add_argument("-f", "--stitch_fits", action="store_true",
                       default=_config_default(config, 'stitch_fits', False),
                       help="Stitch FITS files by extension")
    parser.add_argument("--no-stitch_fits", dest="stitch_fits", action="store_false",
                       help="Disable FITS stitching when enabled by JSON config")

    parser.add_argument("-s", "--save_output", action="store_true",
                       default=_config_default(config, 'save_output', False),
                       help="Save the plots (jpeg), the summary CSV and a config snapshot to the output directory")
    parser.add_argument("--no-save_output", dest="save_output", action="store_false",
                       help="Disable output saving when enabled by JSON config")
    parser.add_argument("-o", "--output_dir", type=str,
                       default=_config_default(config, 'output_dir', None),
                       help="Directory for all saved output (plots, summary CSV, config snapshot). Defaults to a 'plots/' folder alongside the source FITS.")
    parser.add_argument("--show_plots", action="store_true",
                       default=_config_default(config, 'show_plots', True),
                       help="Display plots interactively via plt.show() (default: True)")
    parser.add_argument("--no-show_plots", dest="show_plots", action="store_false",
                       help="Suppress interactive display; useful for batch/headless runs")
    parser.add_argument("-v", "--verbose", action="store_true",
                       default=_config_default(config, 'verbose', False),
                       help="Print verbose output")
    parser.add_argument("--no-verbose", dest="verbose", action="store_false",
                       help="Disable verbose output when enabled by JSON config")
    parser.add_argument("--extra_plot_title", type=str,
                       default=_config_default(config, 'extra_plot_title', '') or '',
                        help="Additional title added at beginning of all plot titles ('<Additional title>+<Default title>)")
    parser.add_argument("--nimages", type=int,
                       default=_config_default(config, 'nimages', None),
                        help="Number of images used in the stack. If not provided, extracted from filename for stitched files.")

    # --- peak finding (same semantics as run_nonlinearity_studies.py: the coarse
    # scipy step shares its --nonlinearity_peakfinder_bin_factor, and
    # peak_finder_widths/buffers and the fallback range are given in electrons.
    # The sub-bin refinement that follows it has no counterpart upstream and gets
    # its own bins-per-ADU key, --sub_bin_peakfinder_bin_factor) ---
    parser.add_argument("--peak_finder_widths", nargs='+', type=float,
                       default=_config_default(config, 'peak_finder_widths', [0.1, 0.1, 0.1, 0.1]),
                        help="Minimum peak width required by scipy.signal.find_peaks, in units of electrons (internally converted to bins at nonlinearity_peakfinder_bin_factor, so it means the same width as in the standard analysis even though this histogram is binned in ADU). Scalar (applied to all extensions) or one value per extension.")
    parser.add_argument("--peak_finder_buffers", nargs='+', type=float,
                       default=_config_default(config, 'peak_finder_buffers', [3.0, 3.0, 3.0, 3.0]),
                        help="Buffer, SUBTRACTED from nonlinearity_peakfinder_bin_factor to compute the minimum neighbor-peak distance: d = nonlinearity_peakfinder_bin_factor - buffer, in bins. Quoted NORMALIZED to 10 bins/e-, exactly as in the standard analysis, so the same number means the same peak separation in both pipelines. Applies to the coarse scipy step only -- the sub-bin refinement has no distance filter. Fractional values are allowed.")
    parser.add_argument("--peak_finder_prominences", nargs='+', type=float,
                       default=_config_default(config, 'peak_finder_prominences', None),
                        help="Minimum peak prominence required by scipy.signal.find_peaks, in histogram counts. Scalar or one value per extension. None disables the filter.")
    parser.add_argument("--nonlinearity_peakfinder_bin_factor", type=int,
                       default=_config_default(config, 'nonlinearity_peakfinder_bin_factor', 10),
                        help="Number of histogram bins per ELECTRON for the coarse scipy.signal.find_peaks step, and the knob for its filters (width = peak_finder_widths * this, distance = this - peak_finder_buffers). Shared with the standard analysis, and quoted in the same units, so the detection step is binned identically in both; this spectrum is in ADU, so it is converted to bins per ADU (divided by ADU_PER_E) for the histogram find_peaks runs on. The all-peaks legend reports it both ways. Higher = finer histogram, but more sensitive to noise. See --sub_bin_peakfinder_bin_factor for the refinement step's binning.")
    parser.add_argument("--sub_bin_peakfinder_bin_factor", type=float,
                       default=_config_default(config, 'sub_bin_peakfinder_bin_factor', 40),
                        help="Number of histogram bins per ADU for the sub-bin localization step, which re-bins the same ADU spectrum and refines each scipy-found peak on it. This step is unique to this module, so it is keyed in ADU (the units the spectrum is actually in) rather than electrons; the all-peaks legend reports it both ways. Normally finer than --nonlinearity_peakfinder_bin_factor to buy localization precision -- if the log-parabola fits start failing, that is the binning that got too fine for --sub_bin_half_width.")
    parser.add_argument("--all_peaks_n_max_electrons", type=float, default=1500,
                        help="Right edge of the ADU histogram, expressed as 'how many electrons out' we want to see -- matches the production pipeline's all_peaks_range_right default (1500 e-). Converted to ADU per extension via pedestal + n_max*gain (using the zero/one fit purely to size the window, never to rescale the data). Ignored if --all_peaks_range_right_adu is set.")
    parser.add_argument("--all_peaks_range_right_adu", type=float, default=None,
                        help="Right edge of the ADU histogram, given directly in ADU. Overrides --all_peaks_n_max_electrons.")
    parser.add_argument("--all_peaks_range_right_adu_fallback", type=float, default=20000.0,
                        help="Right edge of the ADU histogram used only if every extension's zero/one gain fit failed (non-finite), so the run can still proceed.")

    # --- sub-bin peak localization: refines each scipy-found peak (a bin index)
    # to a sub-bin location via the log(counts) parabola vertex -------------
    parser.add_argument("--sub_bin_half_width", type=int,
                        default=_config_default(config, 'sub_bin_half_width', 3),
                        help="Number of bins on each side of a detected peak included in the sub-bin log-parabola fit window (window is 2*half_width+1 bins wide). Bins of the SUB-BIN histogram, so its width in ADU is set together with --sub_bin_peakfinder_bin_factor: raise both to keep the same window while localizing more finely.")
    parser.add_argument("--sub_bin_min_good_bins", type=int,
                        default=_config_default(config, 'sub_bin_min_good_bins', 5),
                        help="Minimum number of in-range, positive-count bins required within the fit window for the sub-bin fit to be trusted; below this the peak falls back to its raw bin center (or the Gaussian refit, if that succeeds).")
    parser.add_argument("--gaussian_peakfind_fallback", action="store_true",
                        default=_config_default(config, 'gaussian_peakfind_fallback', True),
                        help="For peaks where the sub-bin log-parabola fit fails: refit with a seeded single-Gaussian curve_fit instead of leaving them at the raw bin center. Adds a per-peak Python loop over the (typically few-percent) failed subset, so disable if that's too slow to be worth it.")
    parser.add_argument("--no-gaussian_peakfind_fallback", dest="gaussian_peakfind_fallback", action="store_false",
                        help="Disable the Gaussian refit fallback; peaks where the sub-bin log-parabola fit fails stay at the raw bin center.")
    parser.add_argument("--gaussian_refit_k_sigma", type=float,
                        default=_config_default(config, 'gaussian_refit_k_sigma', 3.0),
                        help="For peaks where the sub-bin log-parabola fit fails: half-width of the fallback single-Gaussian curve_fit window, in units of the zero/one-electron sigma (window is +/- k_sigma * sigma_seed).")

    # --- pedestal subtraction: same flags as run_nonlinearity_studies.py -------
    parser.add_argument("--do_pedestal_subtraction", action="store_true",
                        default=_config_default(config, 'do_pedestal_subtraction', True),
                        help="Whether to do pedestal subtraction along one or two axes as preprocess step")
    parser.add_argument("--n_std_to_mask", type=float,
                        default=_config_default(config, 'n_std_to_mask', 1.5),
                        help="Number of standard deviations from mean charge along axis to mask for pedestal subtraction")
    parser.add_argument("--pedsub_max_iter", type=int,
                        default=_config_default(config, 'pedsub_max_iter', 5),
                        help="Max sigma-clip iterations for pedestal subtraction.")
    parser.add_argument("--pedestal_subtraction_axis", type=str,
                        default=_config_default(config, 'pedestal_subtraction_axis', 'row'),
                        help="Which axis to compute pedestals across. Options: 'row', 'col', 'row_then_col','col_then_row'")
    parser.add_argument("--pedsub_cache_dir", type=str,
                        default=_config_default(config, 'pedsub_cache_dir', None),
                        help="Directory to store pedestal-subtracted FITS cache. Defaults to same directory as the source FITS.")
    parser.add_argument("--force_pedsub", action="store_true",
                        default=_config_default(config, 'force_pedsub', False),
                        help="Force recomputation of pedestal subtraction, ignoring any existing cache file.")
    parser.add_argument("--no-force_pedsub", dest="force_pedsub", action="store_false",
                        help="Use cached pedestal-subtracted data when params match (default)")
    parser.add_argument("--use_overscan_only", nargs='*', type=int, metavar='EXT',
                        default=_config_default(config, 'use_overscan_only', False),
                        help="Estimate the per-row pedestal from the overscan columns only, then subtract it from the full frame. Give extension numbers (1-4) to apply to only those extensions, or pass the flag alone to apply to all.")
    parser.add_argument("--no-use_overscan_only", dest="use_overscan_only",
                        action="store_const", const=False,
                        help="Estimate the pedestal from the full frame for all extensions (default).")
    parser.add_argument("--overscan_cols", nargs=2, type=int,
                        default=_config_default(config, 'overscan_cols', list(OVERSCAN_COLS)),
                        metavar=('START', 'STOP'),
                        help="Column range (Python half-open slice START:STOP) the per-row pedestal is estimated from when --use_overscan_only is set.")
    parser.add_argument("--use_biweight_loc", action="store_true",
                        default=_config_default(config, 'use_biweight_loc', True),
                        help="If true, uses Tukey biweight location. If false, uses simple average.")
    parser.add_argument("--use_biweight_midvar", action="store_true",
                        default=_config_default(config, 'use_biweight_midvar', True),
                        help="If true, uses Tukey biweight midvariance. If false, uses standard deviation.")

    # --- zero/one fit: same flags as run_nonlinearity_studies.py (used here only
    # to seed the ADU histogram's pedestal/noise, not to rescale the spectrum) ---
    parser.add_argument("--zero_one_n_bins", type=_n_bins,
                        default=_config_default(config, 'zero_one_n_bins', _ZERO_ONE_N_BINS),
                        help="Strict number of bins spanning each internal zero/one histogram. Integer >= 10.")
    parser.add_argument("--zero_one_window_left_scale", type=_window_scale,
                        default=_config_default(config, 'zero_one_window_left_scale', 1.0),
                        help="Scale the left half-width of the auto-computed zero/one fit window (>=1.0).")
    parser.add_argument("--zero_one_window_right_scale", type=_window_scale,
                        default=_config_default(config, 'zero_one_window_right_scale', 1.0),
                        help="Scale the right half-width of the auto-computed zero/one fit window (>=1.0).")
    parser.add_argument("--fit_cols", nargs='+', type=int,
                        default=_config_default(config, 'fit_cols', None),
                        metavar='COL',
                        help="Restrict the zero/one fit (and its histograms) to image columns (a Python half-open slice START:STOP). Two ints for all extensions, or two per extension. In the JSON config use a [START, STOP] pair or a per-extension list.")
    parser.add_argument("--zero_one_gain_guess", nargs='+', type=float,
                        default=_config_default(config, 'zero_one_gain_guess', None),
                        metavar='GAIN',
                        help="Seed for the one-electron peak location, used to initialize the double-Gaussian fit. One value for all extensions, or one per extension.")

    # --- gain-independent calculations ---------------------------------------
    parser.add_argument("--gain_indep_fit_n_min", type=int,
                        default=_config_default(config, 'gain_indep_fit_n_min', 2),
                        help="Smallest number of peaks (counting the zero-electron peak) used by the running fit g_fit(i). 2 starts the scan at the two-point line fit, which reproduces the standard gain exactly, so the scan covers the whole space including the standard method.")
    parser.add_argument("--gain_indep_fit_n_max", type=int,
                        default=_config_default(config, 'gain_indep_fit_n_max', None),
                        help="Largest number of peaks (counting the zero-electron peak) any fit may use, i.e. fits stop at peak number n_max-1. None (the default) finds it per extension with the same changepoint detector the production pipeline uses for fit_range_right: the fit stops where the clean parabola in q(n) gives way to noise.")
    parser.add_argument("--gain_indep_avg_gain_n_min", type=int,
                        default=_config_default(config, 'gain_indep_avg_gain_n_min', 100),
                        help="Smallest fit range i included when the average gain plot averages and fits a line to g_fit(i). Only affects that plot's reported average and trend line, not the running fit itself or any other fit.")
    parser.add_argument("--gain_indep_avg_gain_n_max", type=int,
                        default=_config_default(config, 'gain_indep_avg_gain_n_max', 1400),
                        help="Largest fit range i included in the average gain plot's average and trend line. None uses every fit range the running fit produced.")
    parser.add_argument("--gain_indep_fit_cubic", action="store_true",
                        default=_config_default(config, 'gain_indep_fit_cubic', True),
                        help="Also fit the ADU charge spectrum to a cubic, and report how far it departs from the quadratic fit (the error of truncating the series at n^2).")
    parser.add_argument("--no-gain_indep_fit_cubic", dest="gain_indep_fit_cubic", action="store_false",
                        help="Skip the cubic fit and the truncation-error estimate.")
    parser.add_argument("--gain_indep_gain_scan", action="store_true",
                        default=_config_default(config, 'gain_indep_gain_scan', True),
                        help="Calibrate the spectrum with every candidate gain g_fit(i) and recompute the standard nonlinearity for each, to find the gain that minimizes it.")
    parser.add_argument("--no-gain_indep_gain_scan", dest="gain_indep_gain_scan", action="store_false",
                        help="Skip the calibration-gain scan (it is the slowest step: one parabola fit per candidate gain).")
    parser.add_argument("--gain_indep_scan_fit_range_right", type=_float_or_auto,
                        default=_config_default(config, 'gain_indep_scan_fit_range_right', 'auto'),
                        help="Right edge (in electrons) of the standard nonlinearity parabola fit performed at each candidate gain during the scan. 'auto' (the default) places it per extension with the changepoint estimator (the same one --fit_range_right auto uses in the production pipeline); a number is applied to every extension.")
    parser.add_argument("--nonlinearity_fit_fix_intercept", action="store_true",
                        default=_config_default(config, 'nonlinearity_fit_fix_intercept', True),
                        help="Hold the nonlinearity parabola's constant term at 0, which Nonlin(0) = 0 requires by construction. Left free it absorbs the misfit of a parabola against a curve that stops being quadratic before fit_range_right, offsetting the fitted curve from the data at low peak numbers.")
    parser.add_argument("--no-nonlinearity_fit_fix_intercept", dest="nonlinearity_fit_fix_intercept", action="store_false",
                        help="Fit the constant term as a free parameter (the pre-fix behavior).")
    parser.add_argument("--gain_indep_compare_test_i", nargs='+', type=int,
                        default=_config_default(config, 'gain_indep_compare_test_i', [10, 500, 1000]),
                        help="List of indices to evaluate the g_fit and alphas at in comparison.")
    parser.add_argument("--gain_indep_scan_eval_charge", type=float,
                        default=_config_default(config, 'gain_indep_scan_eval_charge', 500.0),
                        help="Charge (in electrons) at which the scanned nonlinearity is evaluated and minimized.")

    # --- plot_all_peaks display options: same flags as run_nonlinearity_studies.py
    parser.add_argument("--plot_all_peaks_individual", action="store_true",
                        default=_config_default(config, 'plot_all_peaks_individual', True),
                        help="Plot the ADU all-peaks spectrum as one figure per extension")
    parser.add_argument("--no-plot_all_peaks_individual", dest="plot_all_peaks_individual", action="store_false",
                        help="Disable individual all-peaks plots when enabled by JSON config")
    parser.add_argument("--plot_all_peaks_together", action="store_true",
                        default=_config_default(config, 'plot_all_peaks_together', True),
                        help="Plot the ADU all-peaks spectrum as a combined 2x2 subplot")
    parser.add_argument("--no-plot_all_peaks_together", dest="plot_all_peaks_together", action="store_false",
                        help="Disable combined all-peaks subplot")
    parser.add_argument("--plot_all_peaks_xlim", nargs='+', type=_lim_token,
                        default=_config_default(config, 'plot_all_peaks_xlim', None),
                        metavar='VAL',
                        help="X-axis limits (in ADU) for the all-peaks plot. Provide 2 values (LEFT RIGHT) to apply to all extensions, or 8 values for per-extension limits. Defaults to full histogram range if not set.")
    parser.add_argument("--plot_all_peaks_ylim", nargs='+', type=_lim_token,
                        default=_config_default(config, 'plot_all_peaks_ylim', None),
                        metavar='VAL',
                        help="Y-axis limits for the all-peaks plot. Provide 2 values (BOTTOM TOP) to apply to all extensions, or 8 values for per-extension limits. Defaults to auto if not set.")
    parser.add_argument("--plot_all_peaks_yscale", type=str,
                        default=_config_default(config, 'plot_all_peaks_yscale', 'linear'),
                        help="Y-axis scale for the all-peaks plot. Options: 'linear', 'log'")
    parser.add_argument("--plot_all_peaks_window_bin_factor", type=float,
                        default=_config_default(config, 'plot_all_peaks_window_bin_factor', None),
                        help="If set, re-histogram the all-peaks plot within the displayed --plot_all_peaks_xlim window at this many bins per ADU-equivalent-of-an-electron, instead of reusing the full-spectrum histogram. Requires a concrete --plot_all_peaks_xlim.")
    parser.add_argument("--plot_all_peaks_individual_figsize", nargs=2, type=float,
                        default=_config_default(config, 'plot_all_peaks_individual_figsize', [9, 6.5]),
                        metavar=('W', 'H'),
                        help="Figure size (width height) for individual all-peaks plots")
    parser.add_argument("--plot_all_peaks_subplots_figsize", nargs=2, type=float,
                        default=_config_default(config, 'plot_all_peaks_subplots_figsize', [13, 9]),
                        metavar=('W', 'H'),
                        help="Figure size (width height) for combined 2x2 all-peaks subplot")
    parser.add_argument("--show_titles", action="store_true",
                        default=_config_default(config, 'show_titles', True),
                        help="Show suptitles and per-axis titles on the plot")
    parser.add_argument("--no-show_titles", dest="show_titles", action="store_false",
                        help="Hide all suptitles and per-axis titles")
    parser.add_argument("--plot_all_peaks_sharex", action="store_true",
                        default=_config_default(config, 'plot_all_peaks_sharex', True),
                        help="Share x-axis range across the 2x2 all-peaks subplots")
    parser.add_argument("--no-plot_all_peaks_sharex", dest="plot_all_peaks_sharex", action="store_false",
                        help="Allow independent x-axis ranges per all-peaks subplot")
    parser.add_argument("--plot_all_peaks_sharey", action="store_true",
                        default=_config_default(config, 'plot_all_peaks_sharey', True),
                        help="Share y-axis range across the 2x2 all-peaks subplots")
    parser.add_argument("--no-plot_all_peaks_sharey", dest="plot_all_peaks_sharey", action="store_false",
                        help="Allow independent y-axis ranges per all-peaks subplot")
    parser.add_argument("--raw_peak_number", action="store_true",
                        default=_config_default(config, 'raw_peak_number', True),
                        help="Annotate each detected peak with its index (0, 1, 2, ...) at its raw bin-center location, on every all-peaks plot")
    parser.add_argument("--no-raw_peak_number", dest="raw_peak_number", action="store_false",
                        help="Hide the raw bin-center per-peak index annotations")
    parser.add_argument("--raw_peak_location_vertical_line", action="store_true",
                        default=_config_default(config, 'raw_peak_location_vertical_line', True),
                        help="Draw a vertical line at each detected peak's raw bin-center location, on every all-peaks plot")
    parser.add_argument("--no-raw_peak_location_vertical_line", dest="raw_peak_location_vertical_line", action="store_false",
                        help="Hide the raw bin-center vertical lines")
    parser.add_argument("--sub_bin_peak_number", action="store_true",
                        default=_config_default(config, 'sub_bin_peak_number', True),
                        help="Annotate each detected peak with its index (0, 1, 2, ...) at its sub-bin (log-parabola / Gaussian-refit) location, on every all-peaks plot")
    parser.add_argument("--no-sub_bin_peak_number", dest="sub_bin_peak_number", action="store_false",
                        help="Hide the sub-bin per-peak index annotations")
    parser.add_argument("--sub_bin_peak_location_vertical_line", action="store_true",
                        default=_config_default(config, 'sub_bin_peak_location_vertical_line', True),
                        help="Draw a vertical line at each detected peak's sub-bin location, on every all-peaks plot")
    parser.add_argument("--no-sub_bin_peak_location_vertical_line", dest="sub_bin_peak_location_vertical_line", action="store_false",
                        help="Hide the sub-bin vertical lines")
    parser.add_argument("--peak_number_labels_size", type=float,
                        default=_config_default(config, 'peak_number_labels_size', 8),
                        help="Font size for the per-peak index annotations on the all-peaks plot (applies to both raw and sub-bin labels)")

    # --- gain-independent plots: one switch pair each, shared figure settings ---
    _add_plot_switches(parser, config, 'gain_indep_nonlin',
                       "the differential gain delta_g(n) vs n, with its fitted line")
    _add_plot_switches(parser, config, 'charge_spectrum_fit',
                       "the ADU charge spectrum's departure from a straight line, with the fitted nonlinear term")
    _add_plot_switches(parser, config, 'charge_spectrum',
                       "the ADU charge spectrum q(n) itself, with the parabola fit Q(n) over it "
                       "and a q(n) - Q(n) residual panel")
    _add_plot_switches(parser, config, 'cumulative_delta',
                       "S(n) = I(n)/n vs n, with the fitted line whose slope is alpha'")
    _add_plot_switches(parser, config, 'cumulative_dir_delta',
                       "S(n) built from the direct gain differential, with the alpha' fit line over it")
    _add_plot_switches(parser, config, 'running_gain',
                       "the running average gain g_fit(i) vs i, against the standard gain, "
                       "with its windowed average and trend line and a g_fit(i) - g residual panel")
    _add_plot_switches(parser, config, 'gain_fit_residual',
                       "the running gain residual g_fit(i) - g vs i")
    _add_plot_switches(parser, config, 'gain_fit_residual_cumsum',
                       "the integrated running gain residual vs i, with a parabola fit through it")
    _add_plot_switches(parser, config, 'direct_gain',
                       "the direct gain g_dir(n) vs n, against g and g_fit(i)")
    _add_plot_switches(parser, config, 'direct_delta',
                       "the direct gain differential g - g_dir(n) vs n, with the fitted line through it")
    _add_plot_switches(parser, config, 'gain_scan',
                       "the standard nonlinearity vs the gain used to calibrate")
    _add_plot_switches(parser, config, 'compare_nonlinearity_curves',
                       "the standard nonlinearity curve computed with each of the gains "
                       "(zero/one, sub-bin, g_fit(i), optimal g_fit, direct), overlaid")
    _add_plot_switches(parser, config, 'alpha_alpha_prime',
                       "alpha and alpha' vs the peak number they were fit up to, "
                       "with g_fit(i) on a twinned axis")
    _add_plot_switches(parser, config, 'nonlin_vs_gain',
                       "the standard nonlinearity at a fixed charge vs a swept calibration gain")

    parser.add_argument("--plot_gain_indep_nonlin_individual_figsize", nargs=2, type=float,
                        default=_config_default(config, 'plot_gain_indep_nonlin_individual_figsize', [7, 6]),
                        metavar=('W', 'H'),
                        help="Figure size (width height) for every individual gain-independent plot")
    parser.add_argument("--plot_gain_indep_nonlin_together_figsize", nargs=2, type=float,
                        default=_config_default(config, 'plot_gain_indep_nonlin_together_figsize', [13, 9]),
                        metavar=('W', 'H'),
                        help="Figure size (width height) for every combined 2x2 gain-independent subplot")
    parser.add_argument("--plot_gain_indep_nonlin_sharex", action="store_true",
                        default=_config_default(config, 'plot_gain_indep_nonlin_sharex', True),
                        help="Share x-axis range across the 2x2 subplots of every gain-independent plot")
    parser.add_argument("--no-plot_gain_indep_nonlin_sharex", dest="plot_gain_indep_nonlin_sharex", action="store_false",
                        help="Allow independent x-axis ranges per gain-independent subplot")
    parser.add_argument("--plot_gain_indep_nonlin_sharey", action="store_true",
                        default=_config_default(config, 'plot_gain_indep_nonlin_sharey', True),
                        help="Share y-axis range across the 2x2 subplots of every gain-independent plot")
    parser.add_argument("--no-plot_gain_indep_nonlin_sharey", dest="plot_gain_indep_nonlin_sharey", action="store_false",
                        help="Allow independent y-axis ranges per gain-independent subplot")

    parsed_args = parser.parse_args(args)

    if parsed_args.file_string is None:
        parser.error('file_string is required unless it is provided in the JSON config')

    for _scale_arg in ('zero_one_window_left_scale', 'zero_one_window_right_scale'):
        if getattr(parsed_args, _scale_arg) < 1.0:
            parser.error(f"{_scale_arg} must be >= 1.0 (got {getattr(parsed_args, _scale_arg)})")
    if int(parsed_args.zero_one_n_bins) < 10:
        parser.error(f"zero_one_n_bins must be an integer >= 10 (got {parsed_args.zero_one_n_bins})")
    parsed_args.zero_one_n_bins = int(parsed_args.zero_one_n_bins)

    if parsed_args.gain_indep_fit_n_min < 2:
        parser.error(f"gain_indep_fit_n_min must be >= 2 (got {parsed_args.gain_indep_fit_n_min})")
    if (parsed_args.gain_indep_fit_n_max is not None
            and parsed_args.gain_indep_fit_n_max <= parsed_args.gain_indep_fit_n_min):
        parser.error(f"gain_indep_fit_n_max ({parsed_args.gain_indep_fit_n_max}) must exceed "
                     f"gain_indep_fit_n_min ({parsed_args.gain_indep_fit_n_min})")
    if (parsed_args.gain_indep_avg_gain_n_min is not None
            and parsed_args.gain_indep_avg_gain_n_max is not None
            and parsed_args.gain_indep_avg_gain_n_max <= parsed_args.gain_indep_avg_gain_n_min):
        parser.error(f"gain_indep_avg_gain_n_max ({parsed_args.gain_indep_avg_gain_n_max}) must exceed "
                     f"gain_indep_avg_gain_n_min ({parsed_args.gain_indep_avg_gain_n_min})")

    try:
        parsed_args.fit_cols = _coerce_fit_cols(parsed_args.fit_cols)
    except ValueError as e:
        parser.error(str(e))

    parsed_args.zero_one_gain_guess = _coerce_gain_guess(parsed_args.zero_one_gain_guess)
    _gg = parsed_args.zero_one_gain_guess
    _gg_values = _gg if isinstance(_gg, list) else ([] if _gg is None else [_gg])
    for v in _gg_values:
        if v is not None and float(v) <= 0:
            parser.error(f"zero_one_gain_guess values must be > 0 (got {v})")

    return parsed_args


if __name__ == '__main__':
    args = init_argparse()
    main(args)
