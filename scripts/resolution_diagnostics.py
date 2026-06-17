#!/usr/bin/env python3
"""Quick standalone harness for single-electron *resolution* at a charge q.

This is a thin convenience wrapper for experimenting on one FITS file without a
full JSON config. All the actual work -- the constrained n-Gaussian comb fit,
the goodness-of-fit / verdict, the summary table, and the plots -- now lives in
the package (nonlinearity_studies.nonlinearity_studies); see resolution_at_charge,
resolution_at_charge_ext, classify_resolution, summarize_resolution, and
plot_resolution. This script only loads the data (mirroring the pipeline up
through peak finding) and calls those functions, so there is a single source of
truth for the algorithm.

For production runs use the package CLI instead:
    run-nonlinearity-studies <fits> -r 1000 --resolution_window 10 --save_plots --save_csv
    run-nonlinearity-studies -j config/resolution_study.json

Usage
-----
    python scripts/resolution_diagnostics.py <fits_path> -q 1000 -n 10
    python scripts/resolution_diagnostics.py <fits_path> -q 200 600 1000 --ext 1 2 --together
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from nonlinearity_studies.nonlinearity_studies import (  # noqa: E402
    get_fits,
    pedestal_subtract_ext_cached,
    get_zero_one_peaks_ext,
    get_all_peaks_ext,
    resolution_at_charge_ext,
    summarize_resolution,
    plot_resolution,
)


def load_extensions(fits_path, do_pedsub=None, bin_factor=10, range_right=1500):
    """Load a FITS file and run pedestal subtraction + zero/one + all-peaks, the
    same way the pipeline does. Returns the inputs resolution_at_charge_ext needs."""
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
        print(f'Skipping pedestal subtraction: {fits_path}')

    _, _, pedestals, gains, double_gauss_popts, _ = get_zero_one_peaks_ext(
        data_ext, n=100, fit_bounds='default',
    )
    counts_ext, _, peaks_ext, centers_ext, _ = get_all_peaks_ext(
        data_ext, widths=[0.1, 0.1, 0.1, 0.1], buffers=[3, 3, 3, 3],
        prominences=[None, None, None, None],
        pedestals=pedestals, double_gauss_popts=double_gauss_popts, gains=gains,
        bins='default', flatten=True, do_convert_to_electrons=True,
        range_left='default', range_right=range_right, bin_factor=bin_factor,
        print_values=False,
    )
    return counts_ext, centers_ext, peaks_ext, gains, double_gauss_popts


def main(argv=None):
    p = argparse.ArgumentParser(
        description='Quantify single-electron resolution at charge q via a '
                    'constrained n-Gaussian comb fit (thin wrapper over the package).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('fits_path', type=str, help='Path to FITS/.fz file (or *.pedsub.fits).')
    p.add_argument('-q', '--resolution-at', nargs='+', type=float, required=True,
                   help='Charge value(s) in electrons to evaluate resolution at.')
    p.add_argument('-n', '--window', type=float, default=10.0,
                   help='Window width in electrons (~ number of peaks). Default 10.')
    p.add_argument('--ext', nargs='+', type=int, default=None,
                   help='1-based extension numbers to report/plot (default: all).')
    p.add_argument('--sigma-well', type=float, default=0.25,
                   help='sigma (e-) below which peaks are "well resolved". Default 0.25.')
    p.add_argument('--sigma-limit', type=float, default=0.5,
                   help='sigma (e-) at/above which peaks are "unresolved". Default 0.5.')
    p.add_argument('--min-peak-frac', type=float, default=0.6,
                   help='Detected/expected peak fraction below which the window is unresolved. Default 0.6.')
    p.add_argument('--bin-factor', type=int, default=10,
                   help='Histogram bins per electron (matches package default).')
    p.add_argument('--range-right', type=int, default=1500,
                   help='Right edge (e-) of the all-peaks histogram. Auto-extended past max q.')
    p.add_argument('--no-pedestal-subtraction', dest='pedsub', action='store_false',
                   help='Skip pedestal subtraction (default: auto).')
    p.add_argument('--individual', dest='individual', action='store_true',
                   help='One figure per (extension, charge) (default).')
    p.add_argument('--no-individual', dest='individual', action='store_false',
                   help='Suppress the per-(extension, charge) figures.')
    p.add_argument('--together', dest='together', action='store_true',
                   help='Also produce a combined 2x2 figure per charge.')
    p.add_argument('-o', '--output-dir', type=str, default=None,
                   help='Directory for saved plots. Default: ./plots/resolution.')
    p.add_argument('--no-save', dest='save', action='store_false',
                   help='Do not save plots to disk.')
    p.add_argument('--no-show', dest='show', action='store_false',
                   help='Do not display plots interactively (batch/headless).')
    p.set_defaults(pedsub=None, individual=True, together=False, save=True, show=True)
    args = p.parse_args(argv)

    fits_path = Path(args.fits_path)
    out_dir = Path(args.output_dir) if args.output_dir else Path('plots') / 'resolution'

    max_q = max(args.resolution_at) + args.window / 2.0
    range_right = max(args.range_right, int(np.ceil(max_q)) + 50)

    counts_ext, centers_ext, peaks_ext, gains, double_gauss_popts = load_extensions(
        fits_path, do_pedsub=args.pedsub, bin_factor=args.bin_factor,
        range_right=range_right)

    results_ext = resolution_at_charge_ext(
        counts_ext, centers_ext, peaks_ext, gains, double_gauss_popts,
        charges=args.resolution_at, window=args.window,
        sigma_well=args.sigma_well, sigma_limit=args.sigma_limit,
        min_peak_frac=args.min_peak_frac, verbose=True)

    if args.ext:
        n_ext = len(results_ext)
        results_ext = [results_ext[i - 1] for i in args.ext if 1 <= i <= n_ext]

    summarize_resolution(results_ext)

    if args.individual or args.together:
        plot_resolution(
            results_ext, args.resolution_at,
            sigma_well=args.sigma_well, sigma_limit=args.sigma_limit,
            plot_individual=args.individual, plot_together=args.together,
            save_plots=args.save, show_plots=args.show,
            fig_path=str(out_dir), file=fits_path.name)


if __name__ == '__main__':
    main()
