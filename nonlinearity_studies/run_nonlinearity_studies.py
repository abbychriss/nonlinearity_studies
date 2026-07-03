#!/usr/bin/env python3
"""
Command-line interface for nonlinearity studies analysis.

This script provides functionality to:
- Stitch FITS images by extension
- Fit zero and one electron peaks
- Analyze all electron peaks
- Calculate and visualize nonlinearity

Usage:
As executable:
    ./nonlinearity_studies/run_nonlinearity_studies.py [OPTIONS] <file_path>

As module:
    python -m nonlinearity_studies.run_nonlinearity_studies [OPTIONS] <file_path>
    
Or after pip installation:
    run-nonlinearity-studies [OPTIONS] <file_path>
"""

import numpy as np
import argparse
import glob
import hashlib
import json
from datetime import datetime
from pathlib import Path
import sys
import re
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 12  # bump default for axis labels, tick labels, legends

# Overscan / fit-column / gain-seed plumbing is reused from pedestal_subtract rather
# than re-implemented here: the coercion + per-extension normalization helpers live in
# its CLI module, the FITS-header geometry readers in its core. Importing __main__ as a
# submodule (not as the entry point) leaves its main() guarded and unexecuted.
from pedestal_subtract.core import (
    overscan_cols_from_header,
    get_fits_header,
    _scalar_for_extension,
    plot_charge_per_column,
)
from pedestal_subtract.__main__ import (
    OVERSCAN_COLS,
    _overscan_ext_indices,
    _coerce_fit_cols,
    _normalize_fit_cols,
    _coerce_gain_guess,
)

# Handle imports for both direct execution and module import
if __name__ == "__main__":
    # When run as script, add parent directory to path
    sys.path.insert(0, str(Path(__file__).parent))
    from stitch_fits import stitch_fits
    from nonlinearity_studies import (
        get_fits,
        pedestal_subtract_ext_cached,
        get_zero_one_peaks_ext,
        get_all_peaks_ext,
        get_nonlinearity_ext,
        estimate_optimal_fit_range_right_ext,
        estimate_fit_range_right_by_noise_onset_ext,
        estimate_fit_range_right_changepoint_ext,
        summarize_extensions,
        resolution_at_charge_ext,
        summarize_resolution,
        build_resolution_summary_wide,
        plot_zero_one_peaks,
        plot_all_peaks,
        plot_nonlinearity,
        plot_resolution
    )
else:
    # When imported as module, use relative imports
    from .stitch_fits import stitch_fits
    from .nonlinearity_studies import (
        get_fits,
        pedestal_subtract_ext_cached,
        get_zero_one_peaks_ext,
        get_all_peaks_ext,
        get_nonlinearity_ext,
        estimate_optimal_fit_range_right_ext,
        estimate_fit_range_right_by_noise_onset_ext,
        estimate_fit_range_right_changepoint_ext,
        summarize_extensions,
        resolution_at_charge_ext,
        summarize_resolution,
        build_resolution_summary_wide,
        plot_zero_one_peaks,
        plot_all_peaks,
        plot_nonlinearity,
        plot_resolution
    )

def _derive_data_path(file_path_str):
    """
    Derive the data path from the file path.
    
    Args:
        file_path_str: The file path string (can contain glob patterns like '*')
        
    Returns:
        tuple: (directory_path, image_pattern) where directory_path is the path to search
               and image_pattern is the file pattern to match
    """
    # Remove trailing slashes
    clean_path_str = file_path_str.rstrip('/')
    
    # Split path and extract the image pattern (last component)
    parts = clean_path_str.split('/')
    image_pattern = parts[-1]  # e.g., '*' or '*.fz'
    
    # Directory path is everything before the pattern
    directory_path = '/'.join(parts[:-1])  # e.g., 'examples/images/ten-images'
    
    # Convert to Path and make absolute if relative
    dir_path = Path(directory_path)
    if not dir_path.is_absolute():
        dir_path = Path.cwd() / dir_path

    return dir_path, image_pattern


def _wildcard_free_base(pattern_str):
    """Deepest leading directory of ``pattern_str`` that contains no glob wildcard.

    When the file pattern spans multiple directories (e.g. ``.../PD07*/cds_avg*.fz``),
    the stitched output must be anchored somewhere deterministic rather than inside a
    literal wildcard-named folder (the old ``out_path`` reused ``PD07*`` verbatim, which
    created an actual directory called ``PD07*``). This walks the leading path segments
    up to the first one containing a glob magic character, so a concrete single
    directory still gets its own ``stitched-fits`` subfolder, while a wildcard pattern
    falls back to the common parent above the wildcard. Always returns an absolute path.
    """
    base_parts = []
    for part in Path(pattern_str).parts:
        if glob.has_magic(part):
            break
        base_parts.append(part)
    base = Path(*base_parts) if base_parts else Path('.')
    if not base.is_absolute():
        base = Path.cwd() / base
    return base


CONFIG_KEYS = {
    'file_string',
    'stitch_fits',
    'plot_charge_per_column',
    'plot_charge_per_column_figsize',
    'plot_zero_one_adu',
    'get_nonlinearity_at',
    'resolution_at',
    'resolution_window',
    'resolution_sigma_well',
    'resolution_sigma_limit',
    'resolution_min_peak_frac',
    'plot_resolution_individual',
    'plot_resolution_together',
    'plot_resolution_individual_figsize',
    'plot_resolution_subplots_figsize',
    'plot_resolution_sharex',
    'plot_resolution_sharey',
    'save_plots',
    'save_csv',
    'output_dir',
    'show_plots',
    'verbose',
    'nimages',
    'extra_plot_title',
    'peak_finder_widths',
    'peak_finder_buffers',
    'peak_finder_prominences',
    'bin_factor',
    'fit_range_right',
    'auto_fit_range_tolerance',
    'auto_fit_range_max_charge_percentile',
    'auto_fit_range_min_peaks_in_fit',
    'auto_fit_range_method',
    'changepoint_window',
    'changepoint_factor',
    'changepoint_floor',
    'changepoint_persist',
    'fit_range_confidence_tol',
    'do_pedestal_subtraction',
    'n_std_to_mask',
    'pedsub_max_iter',
    'pedestal_subtraction_axis',
    'pedsub_cache_dir',
    'force_pedsub',
    'use_overscan_only',
    'overscan_cols',
    'use_biweight_loc',
    'use_biweight_midvar',
    'zero_one_n_bins',
    'zero_one_window_left_scale',
    'zero_one_window_right_scale',
    'zero_one_peakfind_density',
    'fit_cols',
    'zero_one_gain_guess',
    'plot_zero_one_individual_figsize',
    'plot_zero_one_subplots_figsize',
    'plot_all_peaks_xlim',
    'plot_all_peaks_ylim',
    'plot_all_peaks_yscale',
    'plot_all_peaks_individual_figsize',
    'plot_all_peaks_subplots_figsize',
    'plot_nonlinearity_individual_figsize',
    'plot_nonlinearity_subplots_figsize',
    'plot_nonlinearity_xlim',
    'plot_nonlinearity_ylim',
    'plot_zero_one_individual',
    'plot_zero_one_together',
    'plot_all_peaks_individual',
    'plot_all_peaks_together',
    'plot_nonlinearity_individual',
    'plot_nonlinearity_together',
    'plot_zero_one_electrons',
    'plot_zero_one_yscale',
    'plot_zero_one_xlim',
    'plot_zero_one_ylim',
    'plot_zero_one_ylim_electrons',
    'electron_fit_mode',
    'show_titles',
    'plot_zero_one_sharex',
    'plot_all_peaks_sharex',
    'plot_nonlinearity_sharex',
    'plot_zero_one_sharey',
    'plot_all_peaks_sharey',
    'peak_number_labels_individual',
    'peak_number_labels_together',
    'peak_number_label_size',
    'plot_nonlinearity_sharey',
}

def _load_json_config(config_path):
    config_path = Path(config_path)
    with config_path.open('r', encoding='utf-8') as config_file:
        config = json.load(config_file)

    if not isinstance(config, dict):
        raise ValueError(f'JSON config must contain an object at the top level: {config_path}')

    config = {key.replace('-', '_'): value for key, value in config.items()}
    unknown_keys = sorted(set(config) - CONFIG_KEYS)
    if unknown_keys:
        raise ValueError(
            f'Unknown JSON config option(s): {", ".join(unknown_keys)}. '
            f'Allowed options are: {", ".join(sorted(CONFIG_KEYS))}'
        )

    return config


def _config_default(config, key, fallback):
    return config[key] if key in config else fallback


def _window_scale(value):
    """argparse type for the zero/one fit-window scale factors: a float >= 1.0.

    Values below 1 shrink the window inside the one-electron peak, leaving the
    fit nothing to fit, so they are rejected up front."""
    f = float(value)
    if f < 1.0:
        raise argparse.ArgumentTypeError(f"must be >= 1.0 (got {f})")
    return f


def _n_bins(value):
    """argparse type for the zero/one histogram bin count: an integer >= 10
    (the double-Gaussian has 6 free parameters, so fewer bins is ill-posed)."""
    f = int(value)
    if f < 10:
        raise argparse.ArgumentTypeError(f"must be an integer >= 10 (got {f})")
    return f


def _peakfind_density(value):
    """argparse type for the peak-finding histogram density (bins per ADU): a
    number >= 1, used only to locate the zero/one peaks, independent of the
    fit/plot bin count set by --zero_one_n_bins."""
    f = float(value)
    if f < 1.0:
        raise argparse.ArgumentTypeError(f"must be >= 1 (got {f})")
    return f


def _normalize_scalar_or_list(value):
    if isinstance(value, list) and len(value) == 1:
        return value[0]
    return value


def _int_or_auto(s):
    """argparse type that accepts an int or the literal string 'auto'."""
    if isinstance(s, str) and s.lower() == 'auto':
        return 'auto'
    return int(s)


def _lim_token(s):
    """argparse type for an xlim/ylim element: a float, or a sentinel string.

    Lets the limits accept the string sentinels 'auto'/'default'/'none' (e.g. from
    a JSON config, where argparse applies this type to a string default) alongside
    numeric bounds.
    """
    if isinstance(s, str) and s.lower() in ('auto', 'default', 'none'):
        return s.lower()
    return float(s)


# Args that should NOT influence the run-identity hash (operational/output flags only).
_RUN_HASH_EXCLUDE = {
    'config', 'json', 'verbose', 'save_plots', 'save_csv', 'output_dir', 'show_plots',
    'pedsub_cache_dir', 'force_pedsub',
    'plot_zero_one_adu', 'plot_zero_one_electrons',
    'plot_zero_one_individual', 'plot_zero_one_together',
    'plot_all_peaks_individual', 'plot_all_peaks_together',
    'plot_nonlinearity_individual', 'plot_nonlinearity_together',
}


def _resolved_args_dict(args):
    """Return a JSON-serializable dict of the args namespace."""
    out = {}
    for k, v in vars(args).items():
        if isinstance(v, Path):
            out[k] = str(v)
        else:
            out[k] = v
    return out


def _run_hash(args, length=8):
    """Compute a short stable hash from the analysis-meaningful args."""
    d = _resolved_args_dict(args)
    for k in _RUN_HASH_EXCLUDE:
        d.pop(k, None)
    payload = json.dumps(d, sort_keys=True, default=str).encode('utf-8')
    return hashlib.sha1(payload).hexdigest()[:length]


def _parse_lim(raw, n_ext=4):
    """Convert a raw xlim/ylim arg to 'default', a single tuple, or a list of per-extension tuples.

    Accepts:
      None                          -> 'default'
      'auto' / 'default' / 'none'   -> passed through (string sentinels)
      [l, r]                        -> (l, r)  applied to all extensions
      [[l1,r1],[l2,r2],...]         -> [(l1,r1), ...] one per extension (JSON format)
      [[l1,r1], null, ...]          -> [(l1,r1), 'none', ...]  null = auto for that extension
      [l1, r1, l2, r2, ...]         -> [(l1,r1), ...] one per extension (CLI flat format)
    """
    if raw is None:
        return 'default'
    # A lone sentinel string from CLI arrives wrapped in a list (nargs='+'); unwrap it.
    if isinstance(raw, (list, tuple)) and len(raw) == 1 and isinstance(raw[0], str):
        raw = raw[0]
    if isinstance(raw, str):
        allowed = {'auto', 'default', 'none'}
        if raw not in allowed:
            raise ValueError(f'xlim/ylim string must be one of {sorted(allowed)}, got {raw!r}')
        return raw
    # Per-extension form: a length-n_ext list whose entries are each a [left, right]
    # pair or null. null means "auto for that extension" (-> 'none'). Detect it from
    # any pair/null entry, so a leading null (e.g. [null, [0,40], ...]) still counts.
    if any(isinstance(item, (list, tuple)) or item is None for item in raw):
        if len(raw) != n_ext:
            raise ValueError(
                f'Expected {n_ext} per-extension [left, right] pairs (null allowed), got {len(raw)}'
            )
        return ['none' if item is None else tuple(item) for item in raw]
    if len(raw) == 2:
        return tuple(raw)
    if len(raw) == 2 * n_ext:
        return [tuple(raw[i*2:(i+1)*2]) for i in range(n_ext)]
    raise ValueError(
        f'xlim/ylim must be 2 values (all extensions) or {2*n_ext} values '
        f'(one [left, right] pair per extension), got {len(raw)}'
    )


def main(args=None):
    """
    The main executable function of the script.
    
    Args:
        args: The Namespace object containing the parsed command-line arguments.
              If None, arguments are parsed from command line.
    """
    if args is None:
        args = init_argparse()

    # If a non-empty extra_plot_title doesn't already end with a separator,
    # append a newline so it sits on its own line above the default title.
    if args.extra_plot_title and not args.extra_plot_title.endswith((' ', '\n')):
        args.extra_plot_title = f'{args.extra_plot_title}\n'

    file_path = Path(args.file_string)

    if not args.stitch_fits and glob.has_magic(args.file_string):
        # A glob pattern (e.g. '..._VDD-20p5_*.fz') must be expanded -- otherwise the
        # literal '*' path never exists and the check below reports "file not found".
        # The single-file (non-stitch) path analyzes one image, so an ambiguous match
        # is an error: narrow the pattern or use --stitch_fits to combine the files.
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
        # The path doesn't exist as typed. Only fall back to a tree-wide search when the
        # user passed a BARE filename (no directory component) -- that is the intended
        # convenience. A path that includes directories (or is absolute) is taken
        # literally: if it doesn't exist we error out rather than silently analyzing a
        # different same-named file elsewhere in the tree. The old code grabbed the FIRST
        # rglob match of the basename, so a wrong or nonexistent path such as
        # 'examples/images/VR-3/stitched-fits/avg_..._stitched.fits' was silently
        # resolved to the identically named file under VR-6 and the wrong image was
        # analyzed.
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

    # Get values from argparse arguments
    do_stitch_images = args.stitch_fits
    do_plot_charge_per_column = args.plot_charge_per_column
    do_plot_zero_one_peaks = args.plot_zero_one_adu or args.plot_zero_one_electrons
    do_plot_all_peaks = args.plot_all_peaks_individual or args.plot_all_peaks_together
    get_nonlinearity_at_charges = args.get_nonlinearity_at

    # Unpack single value from list for cleaner interface
    if get_nonlinearity_at_charges is not None:
        get_nonlinearity_at_charges = _normalize_scalar_or_list(get_nonlinearity_at_charges)
    
    # Right bound for parabolic fit: accepts an int, a list of ints (one per extension),
    # or the literal 'auto' to invoke the data-driven estimator (see estimate_optimal_fit_range_right_ext).
    def _contains_auto(v):
        if isinstance(v, str):
            return v.lower() == 'auto'
        if isinstance(v, list):
            return any(_contains_auto(x) for x in v)
        return False

    if _contains_auto(args.fit_range_right):
        fit_range_right_ext = 'auto'
    else:
        fit_range_right_ext = _normalize_scalar_or_list(args.fit_range_right)

    do_plot_nonlinearity = args.plot_nonlinearity_individual or args.plot_nonlinearity_together
    save_plots = args.save_plots
    save_csv = args.save_csv
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
            # Anchor the output at the deepest wildcard-free directory of the pattern, so a
            # multi-directory glob (e.g. PD07*/...) writes one combined file to a single
            # 'stitched-fits' folder instead of a literal 'PD07*' directory. out_path is
            # absolute, so it overrides stitch_fits's own file_path/out_path join.
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

    # Each run gets its own subfolder keyed by a hash of the resolved args, with a
    # snapshot of the config alongside the plots so the run is fully reproducible.
    run_hash = _run_hash(args)
    fig_path = fig_path / f'{fits_file_path.stem}_{run_hash}'

    if save_plots or save_csv:
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
        outputs = 'Plots' if save_plots else 'CSV output'
        print(f'{outputs} and config snapshot will be saved to {fig_path}')

    image_name = fits_file_path.name

    fits_path = str(fits_file_path)
    print(f'Analyzing image: {fits_path}\n')
    
    # Load data from FITS file
    data_ext = get_fits(fits_path)

    # Charge per column is a RAW-data diagnostic (median charge per column before pedestal
    # subtraction reveals anomalous columns); keep this reference since pedestal subtraction
    # and the fit-column restriction below both rebind data_ext.
    raw_data_ext = data_ext

    # Resolve the per-extension fit-column slice up front (before pedestal subtraction,
    # which rebinds data_ext), so it is available to the fit-column restriction below.
    fit_cols_ext = _normalize_fit_cols(args.fit_cols, len(data_ext))

    # Per-extension overscan setting: each selected extension estimates its per-row
    # pedestal from the overscan columns only (still subtracted from the full frame);
    # the rest estimate from the full frame.
    overscan_exts = _overscan_ext_indices(args.use_overscan_only, len(data_ext))
    if overscan_exts:
        # Prefer the overscan columns computed from the CCD-geometry header keys
        # (PRESCAN, PHYSCOL, NCOL, NCOLPRE, NSBIN); fall back to --overscan_cols /
        # the config value only when the header lacks those keys.
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

    # If user specifies pedestal subtraction process, apply to all data before doing analysis
    # This will subtract the average baseline charge of specified axis (row/column) from that axis but keeps the data in ADU
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

    # Restrict the columns used for the zero/one fit (and the plotted histograms), if
    # configured. Done after pedestal subtraction so the per-row pedestal (and any
    # overscan estimate) still sees the full frame. The fit flattens each extension, so
    # this only changes which pixels enter the zero/one histogram, not the geometry.
    if any(fc is not None for fc in fit_cols_ext):
        data_ext = [
            np.asarray(d)[:, fc[0]:fc[1]] if fc is not None else d
            for d, fc in zip(data_ext, fit_cols_ext)
        ]
        if verbose:
            print(f'Restricting the fit columns per extension: {fit_cols_ext}')

    if fit_range_right_ext != 'auto':
        if not isinstance(fit_range_right_ext, list):
            fit_range_right_ext = [fit_range_right_ext] * len(data_ext)
        elif len(fit_range_right_ext) != len(data_ext):
            print(
                f'\nError: fit_range_right must be a single value, one value per extension '
                f'({len(data_ext)} values for this file), or the literal "auto".\n'
            )
            sys.exit(1)

    # Extract number of stitched images from filename; --nimages arg or config overrides.
    # A single, unstitched image carries no count in its name -> it is exactly one image.
    nimages = None
    match = re.search(r'_(\d+)_stitched', fits_path)
    if match:
        nimages = int(match.group(1))
    elif stitched_dir_name not in fits_file_path.parts:
        nimages = 1
    if args.nimages is not None:
        nimages = args.nimages

    if do_plot_charge_per_column:
        plot_charge_per_column(
            raw_data_ext,
            n_std_to_mask=args.n_std_to_mask,
            fit_cols_ext=fit_cols_ext,
            figsize=tuple(args.plot_charge_per_column_figsize),
            additional_title=args.extra_plot_title,
            show_titles=args.show_titles,
            nimages=nimages,
            verbose=verbose,
            save_plots=save_plots,
            show_plots=args.show_plots,
            fig_path=str(fig_path),
            file=image_name,
            dpi=350,
        )

    if verbose and args.zero_one_gain_guess is not None:
        # Show the per-extension gain seeds actually applied (resolved with the same
        # helper the fit uses), so a single value is reported as the list it broadcasts to.
        gain_seeds_ext = [_scalar_for_extension(args.zero_one_gain_guess, ext, len(data_ext))
                          for ext in range(len(data_ext))]
        print(f'Using guess for gains: {gain_seeds_ext}')

    # Fit zeroth and first electron peaks to double gaussians
    zero_one_counts_ext, zero_one_edges_ext, pedestals, gains, double_gauss_popts, zero_one_ranges = get_zero_one_peaks_ext(
        data_ext,
        n=args.zero_one_n_bins,
        fit_bounds='default',
        window_left_scale=args.zero_one_window_left_scale,
        window_right_scale=args.zero_one_window_right_scale,
        peakfind_density=args.zero_one_peakfind_density,
        gain_seed=args.zero_one_gain_guess,
    )

    # Extend the all-peaks histogram if a resolution charge sits near/above the
    # default right edge, so the [q - n/2, q + n/2] window is fully covered.
    all_peaks_range_right = 1500
    if args.resolution_at is not None:
        res_charges = args.resolution_at if isinstance(args.resolution_at, list) else [args.resolution_at]
        needed = max(res_charges) + args.resolution_window / 2.0 + 50
        all_peaks_range_right = max(all_peaks_range_right, int(np.ceil(needed)))

    # Apply scipy peak finder to find location of every electron peak
    counts_ext, edges_ext, peaks_ext, centers_ext, hist_ranges = get_all_peaks_ext(data_ext,
                                                                                widths=args.peak_finder_widths,
                                                                                buffers=args.peak_finder_buffers,
                                                                                prominences=args.peak_finder_prominences,
                                                                                pedestals=pedestals,
                                                                                double_gauss_popts=double_gauss_popts,
                                                                                gains=gains,
                                                                                bins='default',
                                                                                flatten=True,
                                                                                do_convert_to_electrons=True,
                                                                                range_left='default',
                                                                                range_right=all_peaks_range_right,
                                                                                bin_factor=args.bin_factor,
                                                                                print_values=verbose)

    fit_range_diagnostics = None
    if fit_range_right_ext == 'auto':
        if args.auto_fit_range_method == 'noise_onset':
            if verbose:
                print(f'\nEstimating fit_range_right via noise-onset detection (window={args.noise_onset_window}, factor={args.noise_onset_factor})...')
            fit_range_right_ext = estimate_fit_range_right_by_noise_onset_ext(
                peaks_ext, centers_ext, pedestals, gains,
                do_convert_to_electrons=False,
                window=args.noise_onset_window, factor=args.noise_onset_factor,
                verbose=verbose,
            )
        elif args.auto_fit_range_method == 'var_a':
            if verbose:
                print('\nEstimating optimal fit_range_right per extension (var_a)...')
                if args.auto_fit_range_tolerance is not None:
                    print(f'  tolerance = {args.auto_fit_range_tolerance}')
                if args.auto_fit_range_max_charge_percentile is not None:
                    print(f'  max_charge_percentile = {args.auto_fit_range_max_charge_percentile}')
                if args.auto_fit_range_min_peaks_in_fit is not None:
                    print(f'  min_peaks_in_fit = {args.auto_fit_range_min_peaks_in_fit}')
            fit_range_right_ext = estimate_optimal_fit_range_right_ext(
                peaks_ext, centers_ext, pedestals, gains,
                do_convert_to_electrons=False,
                tolerance=args.auto_fit_range_tolerance,
                max_charge_percentile=args.auto_fit_range_max_charge_percentile,
                min_peaks_in_fit=args.auto_fit_range_min_peaks_in_fit,
                verbose=verbose,
            )
        else:  # 'changepoint' (default)
            if verbose:
                print(f'\nEstimating fit_range_right via changepoint detection '
                      f'(window={args.changepoint_window}, factor={args.changepoint_factor}, '
                      f'floor={args.changepoint_floor}, persist={args.changepoint_persist})...')
            fit_range_right_ext, fit_range_diagnostics = estimate_fit_range_right_changepoint_ext(
                peaks_ext, centers_ext, pedestals, gains,
                do_convert_to_electrons=False,
                win=args.changepoint_window, factor=args.changepoint_factor,
                floor=args.changepoint_floor, persist=args.changepoint_persist,
                confidence_rel_tol=args.fit_range_confidence_tol,
                verbose=verbose,
            )
            low = [d['ext'] for d in fit_range_diagnostics if d.get('confidence') == 'LOW']
            if low and verbose:
                print(f'  WARNING: low-confidence fit_range_right on EXT {low} '
                      f'(changepoint and var(a) disagree by > {args.fit_range_confidence_tol*100:.0f}%); '
                      f'inspect the nonlinearity plot(s).')

    # Persist the auto-estimate diagnostics (incl. per-extension confidence) so a
    # saved run records which extensions need review.
    if save_plots and fit_range_diagnostics is not None:
        diag_path = fig_path / 'fit_range_estimate.json'
        with diag_path.open('w', encoding='utf-8') as f:
            json.dump({
                'method': args.auto_fit_range_method,
                'confidence_rel_tol': args.fit_range_confidence_tol,
                'params': {
                    'window': args.changepoint_window,
                    'factor': args.changepoint_factor,
                    'floor': args.changepoint_floor,
                    'persist': args.changepoint_persist,
                },
                'per_extension': fit_range_diagnostics,
            }, f, indent=2, default=str)
        print(f'Saved fit_range estimate diagnostics to {diag_path}')

    # Fit parabola to nonlinearity curve
    peak_charge_e_ext, charge_minus_npeak_ext, parabola_coeffs, _ = get_nonlinearity_ext(peaks_ext,
                                                                                        centers_ext, 
                                                                                        pedestals, 
                                                                                        gains, 
                                                                                        fit_range_right_ext, 
                                                                                        do_convert_to_electrons=False,
                                                                                        fit_bounds_low=-100, 
                                                                                        fit_bounds_high=100)

    # Single-electron resolution at the requested charge(s): fit a constrained
    # n-Gaussian comb in a window around each q and report sigma (e-), reduced
    # chi^2, delta-AIC vs a single-Gaussian null, and a 3-tier verdict. Computed
    # before the per-extension summary table below so its results (pivoted to one
    # row per extension) can be merged into the same extension_summary.csv.
    resolution_rows, resolution_charges = None, None
    if args.resolution_at is not None:
        resolution_charges = _normalize_scalar_or_list(args.resolution_at)
        resolution_results_ext = resolution_at_charge_ext(
            counts_ext, centers_ext, peaks_ext, gains, double_gauss_popts,
            charges=resolution_charges,
            window=args.resolution_window,
            sigma_well=args.resolution_sigma_well,
            sigma_limit=args.resolution_sigma_limit,
            min_peak_frac=args.resolution_min_peak_frac,
            verbose=verbose,
        )
        # Console table only -- resolution_summary.csv (one row per (ext, q)) is no
        # longer written; its columns are merged into extension_summary.csv instead
        # (see build_resolution_summary_wide / summarize_extensions below).
        summarize_resolution(resolution_results_ext)
        resolution_rows, resolution_charges = build_resolution_summary_wide(
            resolution_results_ext, resolution_charges
        )
        if args.plot_resolution_individual or args.plot_resolution_together:
            plot_resolution(
                resolution_results_ext,
                resolution_charges,
                individual_figsize=tuple(args.plot_resolution_individual_figsize),
                subplots_figsize=tuple(args.plot_resolution_subplots_figsize),
                additional_title=args.extra_plot_title,
                suptitle='Single-Electron Resolution',
                nimages=nimages,
                sigma_well=args.resolution_sigma_well,
                sigma_limit=args.resolution_sigma_limit,
                plot_individual=args.plot_resolution_individual,
                plot_together=args.plot_resolution_together,
                sharex=args.plot_resolution_sharex,
                sharey=args.plot_resolution_sharey,
                show_titles=args.show_titles,
                save_plots=save_plots,
                show_plots=args.show_plots,
                fig_path=str(fig_path),
                file=image_name,
                dpi=350,
            )

    # Per-extension summary table: gain (ADU/e-), noise (e-), nonlinearity at the same
    # charge(s) requested via --get_nonlinearity_at (falling back to 500 e- when none given),
    # and -- when --resolution_at was given -- the resolution columns computed above, merged
    # in by extension. One combined extension_summary.csv covers everything.
    summary_charges = get_nonlinearity_at_charges if get_nonlinearity_at_charges is not None else 500
    summarize_extensions(
        gains,
        double_gauss_popts,
        parabola_coeffs,
        nonlinearity_charges=summary_charges,
        resolution_rows=resolution_rows,
        resolution_charges=resolution_charges,
        verbose=verbose,
        save_path=(fig_path / 'extension_summary.csv') if save_csv else None,
    )

    # Fit a double gaussian to zero + 1 electron peak in each extension
    if do_plot_zero_one_peaks:
        plot_zero_one_peaks(data_ext,
                            zero_one_counts_ext,
                            zero_one_edges_ext,
                            pedestals,
                            gains,
                            double_gauss_popts,
                            zero_one_ranges,
                            individual_figsize=tuple(args.plot_zero_one_individual_figsize),
                            subplots_figsize=tuple(args.plot_zero_one_subplots_figsize),
                            xlim=args.plot_zero_one_xlim,
                            ylim=args.plot_zero_one_ylim,
                            ylim_electrons=args.plot_zero_one_ylim_electrons,
                            electron_fit_mode=args.electron_fit_mode or 'transform',
                            additional_title=f'{args.extra_plot_title}' if args.extra_plot_title else '',
                            suptitle='Double-Gaussian Fit to Zero-One Electron Peaks',
                            nimages=nimages,
                            yscale=args.plot_zero_one_yscale,
                            fontsize=12,
                            n=100,
                            do_plot_adu=args.plot_zero_one_adu,
                            do_convert_to_electrons=args.plot_zero_one_electrons,
                            plot_individual=args.plot_zero_one_individual,
                            plot_together=args.plot_zero_one_together,
                            sharex=args.plot_zero_one_sharex,
                            sharey=args.plot_zero_one_sharey,
                            show_titles=args.show_titles,
                            save_plots=save_plots,
                            show_plots=args.show_plots,
                            fig_path=str(fig_path),
                            file=image_name,
                            dpi=350)

    if do_plot_all_peaks:
        # Only override the all-peaks x/y limits when the user passes them on the CLI.
        # When a limit flag is absent, omit it so plot_all_peaks' own default applies
        # (its xlim default is the 500-510 e- window) instead of forcing the full
        # histogram range here.
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
                    plot_individual=args.plot_all_peaks_individual,
                    plot_together=args.plot_all_peaks_together,
                    sharex=args.plot_all_peaks_sharex,
                    sharey=args.plot_all_peaks_sharey,
                    show_titles=args.show_titles,
                    draw_lines=True,
                    linecolor='r',
                    linestyle='--',
                    peak_number_labels_individual=args.peak_number_labels_individual,
                    peak_number_labels_together=args.peak_number_labels_together,
                    peak_number_label_size=args.peak_number_label_size,
                    individual_figsize=tuple(args.plot_all_peaks_individual_figsize),
                    subplots_figsize=tuple(args.plot_all_peaks_subplots_figsize),
                    additional_title=args.extra_plot_title,
                    suptitle='Peaks in Pixel Charge Distribution',
                    nimages=nimages,
                    save_plots=save_plots,
                    show_plots=args.show_plots,
                    fig_path=str(fig_path),
                    file=image_name,
                    dpi=350)

    if do_plot_nonlinearity:
        plot_nonlinearity(peaks_ext,
                        parabola_coeffs,
                        peak_charge_e_ext,
                        charge_minus_npeak_ext,
                        fit_range_right_ext,
                        xlim=_parse_lim(args.plot_nonlinearity_xlim, n_ext=len(data_ext)),
                        ylim=_parse_lim(args.plot_nonlinearity_ylim, n_ext=len(data_ext)),
                        individual_figsize=tuple(args.plot_nonlinearity_individual_figsize),
                        subplots_figsize=tuple(args.plot_nonlinearity_subplots_figsize),
                        additional_title=args.extra_plot_title,
                        suptitle='Pixel Charge Nonlinearity',
                        nimages=nimages,
                        line_color='r',
                        scatter_color='b',
                        s=2,
                        alpha=1,
                        plot_individual=args.plot_nonlinearity_individual,
                        plot_together=args.plot_nonlinearity_together,
                        sharex=args.plot_nonlinearity_sharex,
                        sharey=args.plot_nonlinearity_sharey,
                        show_titles=args.show_titles,
                        save_plots=save_plots,
                        show_plots=args.show_plots,
                        fig_path=str(fig_path),
                        file=image_name,
                        dpi=350)
        

def init_argparse(args=None):
    """
    Initializes the ArgumentParser object and defines arguments.
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
        description="""Run nonlinearity analysis pipeline.

This script can:
- Stitch FITS images by extension
- Fit to zeroth/first electron peaks to compute pedestal, noise, gain
- Fit and plot all electron peaks
- Compute and plot nonlinearity
                                    
You can enable any combination of steps using flags below.""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-j", "--json", type=str, default=pre_args.config,
                        help="Path to JSON file containing command-line arguments. Explicit CLI arguments override JSON values.")
    parser.add_argument('file_string', type=str, nargs='?',
                       default=_config_default(config, 'file_string', None),
                       help='Absolute or relative path to image file (.fz or .fits accepted) if not stitching (ex: data/avg_img.fz), or to directory with images if stitching (can include glob pattern like "*", ex: data/03-12-2026/avg*.fz or data/*). May also be set in JSON config.')
    parser.add_argument("-f", "--stitch_fits", action="store_true",
                       default=_config_default(config, 'stitch_fits', False),
                       help="Stitch FITS files by extension")
    parser.add_argument("--no-stitch_fits", dest="stitch_fits", action="store_false",
                       help="Disable FITS stitching when enabled by JSON config")
    parser.add_argument("--plot_charge_per_column", action="store_true",
                       default=_config_default(config, 'plot_charge_per_column', False),
                       help="Plot the median charge per column for each extension (2x2 grid) on the "
                            "raw, pre-pedestal-subtraction data -- a diagnostic for anomalous columns. "
                            "A column whose median is >= --n_std_to_mask biweight-SDs from the "
                            "extension's biweight location is flagged red. Columns excluded by "
                            "--fit_cols are shaded light grey.")
    parser.add_argument("--no-plot_charge_per_column", dest="plot_charge_per_column", action="store_false",
                       help="Disable the charge-per-column plot when enabled by JSON config")
    parser.add_argument("-z", "--plot_zero_one_adu", action="store_true",
                       default=_config_default(config, 'plot_zero_one_adu', False),
                       help="Plot fits to zero/one electron peaks in ADU. Combined with --plot_zero_one_electrons to also (or only) produce the electron-units version.")
    parser.add_argument("--no-plot_zero_one_adu", dest="plot_zero_one_adu", action="store_false",
                       help="Disable ADU zero/one peak plotting when enabled by JSON config")
    parser.add_argument("-g", "--get_nonlinearity_at", nargs='+', type=float,
                       default=_config_default(config, 'get_nonlinearity_at', None),
                       help="Estimate nonlinearity at specified charge value(s) using parabolic fit")
    parser.add_argument("-r", "--resolution_at", nargs='+', type=float,
                       default=_config_default(config, 'resolution_at', None),
                       help="Quantify single-electron resolution at the given charge value(s) in electrons. "
                            "For each q, fits a constrained n-Gaussian comb in a window [q - n/2, q + n/2] "
                            "and reports sigma (e-), reduced chi^2, delta-AIC vs a single-Gaussian null, and a verdict.")
    parser.add_argument("--resolution_window", type=float,
                       default=_config_default(config, 'resolution_window', 10.0),
                       help="Window width n (in electrons, ~ number of peaks) for the resolution comb fit. Default 10.")
    parser.add_argument("--resolution_sigma_well", type=float,
                       default=_config_default(config, 'resolution_sigma_well', 0.25),
                       help="sigma (e-) below which peaks are 'well resolved' (separation > 4 sigma). Default 0.25.")
    parser.add_argument("--resolution_sigma_limit", type=float,
                       default=_config_default(config, 'resolution_sigma_limit', 0.5),
                       help="sigma (e-) at/above which peaks are 'unresolved' (no central valley, separation < 2 sigma). Default 0.5.")
    parser.add_argument("--resolution_min_peak_frac", type=float,
                       default=_config_default(config, 'resolution_min_peak_frac', 0.6),
                       help="Detected/expected peak fraction below which the window is 'unresolved' regardless of sigma. Default 0.6.")
    parser.add_argument("--plot_resolution_individual", action="store_true",
                       default=_config_default(config, 'plot_resolution_individual', False),
                       help="Plot one resolution-window figure per (extension, charge)")
    parser.add_argument("--no-plot_resolution_individual", dest="plot_resolution_individual", action="store_false",
                       help="Disable individual resolution-window plots when enabled by JSON config")
    parser.add_argument("--plot_resolution_together", action="store_true",
                       default=_config_default(config, 'plot_resolution_together', True),
                       help="Plot resolution windows as a combined 2x2 subplot (one per extension) per charge")
    parser.add_argument("--no-plot_resolution_together", dest="plot_resolution_together", action="store_false",
                       help="Disable the combined resolution-window subplot")
    parser.add_argument("--plot_resolution_individual_figsize", nargs=2, type=float,
                       default=_config_default(config, 'plot_resolution_individual_figsize', [8, 6.5]),
                       metavar=('W', 'H'),
                       help="Figure size (width height) for individual resolution-window plots")
    parser.add_argument("--plot_resolution_subplots_figsize", nargs=2, type=float,
                       default=_config_default(config, 'plot_resolution_subplots_figsize', [13, 9]),
                       metavar=('W', 'H'),
                       help="Figure size (width height) for the combined 2x2 resolution-window subplot")
    parser.add_argument("--plot_resolution_sharex", action="store_true",
                       default=_config_default(config, 'plot_resolution_sharex', False),
                       help="Share x-axis range across the 2x2 resolution-window subplots")
    parser.add_argument("--no-plot_resolution_sharex", dest="plot_resolution_sharex", action="store_false",
                       help="Allow independent x-axis ranges per resolution-window subplot (default)")
    parser.add_argument("--plot_resolution_sharey", action="store_true",
                       default=_config_default(config, 'plot_resolution_sharey', False),
                       help="Share y-axis range across the 2x2 resolution-window subplots")
    parser.add_argument("--no-plot_resolution_sharey", dest="plot_resolution_sharey", action="store_false",
                       help="Allow independent y-axis ranges per resolution-window subplot (default)")
    parser.add_argument("-s", "--save_plots", action="store_true",
                       default=_config_default(config, 'save_plots', False),
                       help="Save all plots as jpeg images")
    parser.add_argument("--no-save_plots", dest="save_plots", action="store_false",
                       help="Disable plot saving when enabled by JSON config")
    parser.add_argument("--save_csv", action="store_true",
                       default=_config_default(config, 'save_csv', False),
                       help="Save the per-extension summary table as a CSV (extension_summary.csv) in the output directory")
    parser.add_argument("--no-save_csv", dest="save_csv", action="store_false",
                       help="Disable summary CSV saving when enabled by JSON config")
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
    parser.add_argument("--peak_finder_widths", nargs='+', type=float,
                       default=_config_default(config, 'peak_finder_widths', [0.1, 0.1, 0.1, 0.1]),
                        help="Minimum peak width required by scipy.signal.find_peaks, in units of electrons (internally multiplied by bin_factor to convert to bins). Scalar (applied to all extensions) or one value per extension. Higher = stricter filter on narrow noise spikes.")
    parser.add_argument("--peak_finder_buffers", nargs='+', type=int,
                       default=_config_default(config, 'peak_finder_buffers', [3, 3, 3, 3]),
                        help="Buffer (in bins), SUBTRACTED from bin_factor to compute the minimum neighbor-peak distance: d = bin_factor - buffer. With bin_factor=10: buffer=0 -> 10 bins (1 e- spacing, physical), buffer=3 -> 7 bins (~0.7 e-, loose), buffer=-2 -> 12 bins (1.2 e-, strict). Larger buffer = looser, smaller/negative = stricter.")
    parser.add_argument("--peak_finder_prominences", nargs='+', type=float,
                       default=_config_default(config, 'peak_finder_prominences', None),
                        help="Minimum peak prominence required by scipy.signal.find_peaks, in histogram counts (same units as the y-axis of plot_all_peaks). Scalar or one value per extension. None disables the filter. Larger = stricter. Often the most robust filter for separating real electron peaks from noise.")
    parser.add_argument("--bin_factor", type=int,
                       default=_config_default(config, 'bin_factor', 10),
                        help="Number of histogram bins per electron used in the all-peaks histogram. Also controls peak_finder_widths conversion (electrons -> bins) and the buffer math (distance = bin_factor - buffer). Higher = finer resolution but more sensitive to noise.")
    parser.add_argument("--fit_range_right", nargs='+', type=_int_or_auto,
                       default=_config_default(config, 'fit_range_right', 'auto'),
                        help="Right charge bound (in electrons) for the parabolic nonlinearity fit. Accepts: a single int applied to all extensions (e.g. 500), one int per extension (e.g. 600 850 750 1050), or the literal 'auto' to enable the data-driven estimator that picks the value minimizing var(a) of the parabola fit per extension.")
    parser.add_argument("--auto_fit_range_tolerance", type=float,
                       default=_config_default(config, 'auto_fit_range_tolerance', None),
                        help="(Only with --fit_range_right auto) If set, among candidates whose score is within (1+tolerance)*min_score, pick the smallest charge. Prevents overestimation when the covariance curve is flat/noisy. E.g. 0.10 = 10%%.")
    parser.add_argument("--auto_fit_range_max_charge_percentile", type=float,
                       default=_config_default(config, 'auto_fit_range_max_charge_percentile', None),
                        help="(Only with --fit_range_right auto) If set, cap max_charge at this percentile of detected peak charges (e.g. 90). Keeps candidates out of the noisy upper tail.")
    parser.add_argument("--auto_fit_range_min_peaks_in_fit", type=int,
                       default=_config_default(config, 'auto_fit_range_min_peaks_in_fit', None),
                        help="(Only with --fit_range_right auto) If set, reject candidate fit_range_right values that include fewer than this many peaks below them. Prevents picking absurdly small ranges.")
    parser.add_argument("--auto_fit_range_method", type=str,
                        default=_config_default(config, 'auto_fit_range_method', 'changepoint'),
                        choices=['changepoint', 'var_a', 'noise_onset'],
                        help="Estimator to use when --fit_range_right auto. 'changepoint' (default) is a two-stage local-roughness changepoint detector that finds where the clean parabola gives way to noise and is robust to both var_a failure modes. 'var_a' minimizes parabola pcov[0,0]. 'noise_onset' (experimental) walks forward and returns the first charge where rolling MAD exceeds a clean-region baseline.")
    parser.add_argument("--changepoint_window", type=int,
                        default=_config_default(config, 'changepoint_window', 25),
                        help="(Changepoint method) window size in peaks for the local-roughness quadratic fit / MAD. Odd values recommended.")
    parser.add_argument("--changepoint_factor", type=float,
                        default=_config_default(config, 'changepoint_factor', 4.0),
                        help="(Changepoint method) local roughness must exceed factor * baseline roughness (or the absolute floor, whichever is larger) to mark the noise onset. Lower = more sensitive.")
    parser.add_argument("--changepoint_floor", type=float,
                        default=_config_default(config, 'changepoint_floor', 0.15),
                        help="(Changepoint method) absolute roughness floor in electrons. Prevents an ultra-quiet early region from setting an impossibly tight threshold. Sits between clean scatter (~0.02-0.05) and noisy roughness (~0.3-0.5).")
    parser.add_argument("--changepoint_persist", type=int,
                        default=_config_default(config, 'changepoint_persist', 4),
                        help="(Changepoint method) number of consecutive points that must exceed the threshold to count as the noise onset. Higher = more robust to isolated outliers / precursor bumps.")
    parser.add_argument("--fit_range_confidence_tol", type=float,
                        default=_config_default(config, 'fit_range_confidence_tol', 0.15),
                        help="(Changepoint method) relative tolerance for the var(a) cross-check. Extensions where changepoint and var(a) disagree by more than this fraction are flagged 'LOW' confidence for review. E.g. 0.15 = 15%%.")
    parser.add_argument("--noise_onset_window", type=int, default=30,
                        help="(Noise-onset method) sliding window size in peaks for the local quadratic fit / MAD.")
    parser.add_argument("--noise_onset_factor", type=float, default=2.5,
                        help="(Noise-onset method) rolling MAD must exceed factor * baseline MAD to count as noise onset. Lower = more sensitive.")
    parser.add_argument("--do_pedestal_subtraction", action="store_true",
                        default=_config_default(config, 'do_pedestal_subtraction', True),
                        help="Whether to do pedestal subtraction along one or two axes as preprocess step")
    parser.add_argument("--n_std_to_mask", type=float,
                        default=_config_default(config, 'n_std_to_mask', 1.5),
                        help="Number of standard deviations from mean charge along axis to mask for pedestal subtraction")
    parser.add_argument("--pedsub_max_iter", type=int,
                        default=_config_default(config, 'pedsub_max_iter', 5),
                        help="Max sigma-clip iterations for pedestal subtraction. Each pass re-estimates the pedestal from the clipped zero-peak core; stops early once it converges.")
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
                        help="Estimate the per-row pedestal from the overscan columns "
                             "(see --overscan_cols) only, then subtract it from the full "
                             "frame. Give extension numbers (1-4) to apply to only those "
                             "extensions, or pass the flag alone to apply to all. In the "
                             "JSON config use true/false or a list like [1, 3]. "
                             "(Put the FITS path before this flag, or in the config, so it "
                             "isn't read as an extension number.)")
    parser.add_argument("--no-use_overscan_only", dest="use_overscan_only",
                        action="store_const", const=False,
                        help="Estimate the pedestal from the full frame for all extensions (default).")
    parser.add_argument("--overscan_cols", nargs=2, type=int,
                        default=_config_default(config, 'overscan_cols', list(OVERSCAN_COLS)),
                        metavar=('START', 'STOP'),
                        help="Column range (Python half-open slice START:STOP) the per-row "
                             "pedestal is estimated from when --use_overscan_only is set. "
                             "Negative endpoints count from the right. Used only as a fallback "
                             "when the FITS header lacks the CCD-geometry keys. Default: the "
                             "last 147 columns ([-147:]).")
    parser.add_argument("--use_biweight_loc", action="store_true",
                        default=_config_default(config, 'use_biweight_loc', True),
                        help="If true, uses Tukey biweight location (more robust - iteratively gets rid of outliers). If false, uses simple average.")
    parser.add_argument("--use_biweight_midvar", action="store_true",
                        default=_config_default(config, 'use_biweight_midvar', True),
                        help="If true, uses Tukey biweight midvariance. If false, uses standard deviation.")
    parser.add_argument("--zero_one_n_bins", type=_n_bins,
                        default=_config_default(config, 'zero_one_n_bins', 100),
                        help="Number of bins spanning the zero/one fit window at window "
                             "scale 1.0 (the count scales up automatically when the window "
                             "is widened, keeping bin width constant). Integer >= 10. Default 100.")
    parser.add_argument("--zero_one_window_left_scale", type=_window_scale,
                        default=_config_default(config, 'zero_one_window_left_scale', 1.0),
                        help="Scale the left half-width of the auto-computed zero/one fit "
                             "window (>=1.0; >1 widens toward lower charge). Default 1.0.")
    parser.add_argument("--zero_one_window_right_scale", type=_window_scale,
                        default=_config_default(config, 'zero_one_window_right_scale', 1.0),
                        help="Scale the right half-width of the auto-computed zero/one fit "
                             "window (>=1.0; >1 widens past the one-electron peak). Default 1.0.")
    parser.add_argument("--zero_one_peakfind_density", type=_peakfind_density,
                        default=_config_default(config, 'zero_one_peakfind_density', 10),
                        help="Bins-per-ADU of the internal histograms used to LOCATE the zero/one "
                             "peaks (separate from --zero_one_n_bins, which sets the fit/plot bins). "
                             "Raise for finer detection, lower to aggregate sparse low-statistics "
                             "hits. Number >= 1. Default 10.")
    parser.add_argument("--fit_cols", nargs='+', type=int,
                        default=_config_default(config, 'fit_cols', None),
                        metavar='COL',
                        help="Restrict the zero/one fit (and the plotted histograms) to image "
                             "columns (a Python half-open slice START:STOP). Pass two ints "
                             "(START STOP) to apply one range to every extension, or two per "
                             "extension (e.g. 8 ints for 4 extensions) for per-extension ranges. "
                             "Negative endpoints count from the right. Applied after pedestal "
                             "subtraction, so the pedestal still uses the full frame/overscan. Omit "
                             "to use all columns. In the JSON config use a [START, STOP] pair, or a "
                             "per-extension list of [START, STOP]/null (null = all columns for that "
                             "extension; null endpoints like [256, null] give an open-ended slice).")
    parser.add_argument("--zero_one_gain_guess", nargs='+', type=float,
                        default=_config_default(config, 'zero_one_gain_guess', None),
                        metavar='GAIN',
                        help="Seed for the one-electron peak location -- a guess for the gain "
                             "(ADU/e-) -- used to initialize the double-Gaussian fit instead of "
                             "auto-detecting the one-electron bump. Pass one value to apply it to "
                             "every extension, or one per extension (e.g. 4 values for 4 "
                             "extensions). In the JSON config use a single number or a "
                             "per-extension list of numbers/null (null = auto-detect that "
                             "extension). Omit to auto-detect for all. Values must be > 0.")
    parser.add_argument("--plot_zero_one_individual_figsize", nargs=2, type=float,
                        default=_config_default(config, 'plot_zero_one_individual_figsize', [7, 5]),
                        metavar=('W', 'H'),
                        help="Figure size (width height) for individual zero/one peak plots")
    parser.add_argument("--plot_zero_one_subplots_figsize", nargs=2, type=float,
                        default=_config_default(config, 'plot_zero_one_subplots_figsize', [13, 9]),
                        metavar=('W', 'H'),
                        help="Figure size (width height) for combined 2x2 zero/one peak subplot")
    parser.add_argument("--plot_charge_per_column_figsize", nargs=2, type=float,
                        default=_config_default(config, 'plot_charge_per_column_figsize', [13, 9]),
                        metavar=('W', 'H'),
                        help="Figure size (width height) for the charge-per-column plot")
    parser.add_argument("--plot_all_peaks_xlim", nargs='+', type=_lim_token,
                        default=_config_default(config, 'plot_all_peaks_xlim', None),
                        metavar='VAL',
                        help="X-axis limits for all-peaks plot. Provide 2 values (LEFT RIGHT) to apply to all extensions, or 8 values (L1 R1 L2 R2 L3 R3 L4 R4) for per-extension limits. Defaults to full histogram range if not set.")
    parser.add_argument("--plot_all_peaks_ylim", nargs='+', type=_lim_token,
                        default=_config_default(config, 'plot_all_peaks_ylim', None),
                        metavar='VAL',
                        help="Y-axis limits for all-peaks plot. Provide 2 values (BOTTOM TOP) to apply to all extensions, or 8 values (B1 T1 B2 T2 B3 T3 B4 T4) for per-extension limits. Defaults to auto if not set.")
    parser.add_argument("--plot_all_peaks_yscale", type=str,
                        default=_config_default(config, 'plot_all_peaks_yscale', 'linear'),
                        help="Y-axis scale for all-peaks plot. Options: 'linear', 'log'")
    parser.add_argument("--plot_all_peaks_individual_figsize", nargs=2, type=float,
                        default=_config_default(config, 'plot_all_peaks_individual_figsize', [7, 6]),
                        metavar=('W', 'H'),
                        help="Figure size (width height) for individual all-peaks plots")
    parser.add_argument("--plot_all_peaks_subplots_figsize", nargs=2, type=float,
                        default=_config_default(config, 'plot_all_peaks_subplots_figsize', [13, 9]),
                        metavar=('W', 'H'),
                        help="Figure size (width height) for combined 2x2 all-peaks subplot")
    parser.add_argument("--plot_nonlinearity_individual_figsize", nargs=2, type=float,
                        default=_config_default(config, 'plot_nonlinearity_individual_figsize', [6, 5]),
                        metavar=('W', 'H'),
                        help="Figure size (width height) for individual nonlinearity plots")
    parser.add_argument("--plot_nonlinearity_subplots_figsize", nargs=2, type=float,
                        default=_config_default(config, 'plot_nonlinearity_subplots_figsize', [13, 9]),
                        metavar=('W', 'H'),
                        help="Figure size (width height) for combined 2x2 nonlinearity subplot")
    parser.add_argument("--plot_nonlinearity_xlim", nargs='+', type=_lim_token,
                        default=_config_default(config, 'plot_nonlinearity_xlim', None),
                        metavar='VAL',
                        help="X-axis limits for nonlinearity plot. Provide 2 values (LEFT RIGHT) to apply to all extensions, or 8 values (L1 R1 L2 R2 L3 R3 L4 R4) for per-extension limits. Defaults to auto if not set.")
    parser.add_argument("--plot_nonlinearity_ylim", nargs='+', type=_lim_token,
                        default=_config_default(config, 'plot_nonlinearity_ylim', None),
                        metavar='VAL',
                        help="Y-axis limits for nonlinearity plot. Provide 2 values (BOTTOM TOP) to apply to all extensions, or 8 values (B1 T1 B2 T2 B3 T3 B4 T4) for per-extension limits. Defaults to auto if not set.")
    parser.add_argument("--plot_zero_one_individual", action="store_true",
                        default=_config_default(config, 'plot_zero_one_individual', False),
                        help="Plot zero/one peaks as one figure per extension")
    parser.add_argument("--no-plot_zero_one_individual", dest="plot_zero_one_individual", action="store_false",
                        help="Disable individual zero/one peak plots when enabled by JSON config")
    parser.add_argument("--plot_zero_one_together", action="store_true",
                        default=_config_default(config, 'plot_zero_one_together', True),
                        help="Plot zero/one peaks as a combined 2x2 subplot")
    parser.add_argument("--no-plot_zero_one_together", dest="plot_zero_one_together", action="store_false",
                        help="Disable combined zero/one peak subplot")
    parser.add_argument("--plot_all_peaks_individual", action="store_true",
                        default=_config_default(config, 'plot_all_peaks_individual', False),
                        help="Plot all-peaks as one figure per extension")
    parser.add_argument("--no-plot_all_peaks_individual", dest="plot_all_peaks_individual", action="store_false",
                        help="Disable individual all-peaks plots when enabled by JSON config")
    parser.add_argument("--plot_all_peaks_together", action="store_true",
                        default=_config_default(config, 'plot_all_peaks_together', False),
                        help="Plot all-peaks as a combined 2x2 subplot")
    parser.add_argument("--no-plot_all_peaks_together", dest="plot_all_peaks_together", action="store_false",
                        help="Disable combined all-peaks subplot")
    parser.add_argument("--plot_nonlinearity_individual", action="store_true",
                        default=_config_default(config, 'plot_nonlinearity_individual', False),
                        help="Plot nonlinearity as one figure per extension")
    parser.add_argument("--no-plot_nonlinearity_individual", dest="plot_nonlinearity_individual", action="store_false",
                        help="Disable individual nonlinearity plots when enabled by JSON config")
    parser.add_argument("--plot_nonlinearity_together", action="store_true",
                        default=_config_default(config, 'plot_nonlinearity_together', False),
                        help="Plot nonlinearity as a combined 2x2 subplot")
    parser.add_argument("--no-plot_nonlinearity_together", dest="plot_nonlinearity_together", action="store_false",
                        help="Disable combined nonlinearity subplot")
    parser.add_argument("--plot_zero_one_electrons", action="store_true",
                        default=_config_default(config, 'plot_zero_one_electrons', False),
                        help="Also produce zero/one peak plots converted to electrons (in addition to ADU)")
    parser.add_argument("--no-plot_zero_one_electrons", dest="plot_zero_one_electrons", action="store_false",
                        help="Disable the electron-units zero/one peak plots")
    parser.add_argument("--plot_zero_one_yscale", type=str,
                        default=_config_default(config, 'plot_zero_one_yscale', 'linear'),
                        help="Y-axis scale for zero/one peak plots. Options: 'linear', 'log'")
    # argparse runs `type` over string defaults, which would choke on the 'default'
    # / 'none' keywords, so default to None here and fall back to the config below.
    parser.add_argument("--plot_zero_one_xlim", nargs=2, type=float, default=None,
                        metavar=('LOW', 'HIGH'),
                        help="X-axis limits for the zero/one peak plots. "
                             "Omit (or set 'default') to use the fit range; 'none' to autoscale.")
    parser.add_argument("--plot_zero_one_ylim", nargs=2, type=float, default=None,
                        metavar=('LOW', 'HIGH'),
                        help="Y-axis limits for the ADU zero/one peak plots. "
                             "Omit (or set 'default') for the auto range; 'none' to autoscale.")
    parser.add_argument("--plot_zero_one_ylim_electrons", nargs=2, type=float, default=None,
                        metavar=('LOW', 'HIGH'),
                        help="Y-axis limits for the electron-units zero/one peak plots. Separate from "
                             "--plot_zero_one_ylim since the electron peaks are taller (smaller sigma). "
                             "Omit (or set 'default') for the auto range; 'none' to autoscale.")
    parser.add_argument("--electron_fit_mode", type=str, choices=['transform', 'refit'],
                        default=_config_default(config, 'electron_fit_mode', None),
                        help="How to obtain the electron-unit double-Gaussian curve on the zero/one "
                             "peak plots: 'transform' (default) analytically rescales the converged ADU "
                             "fit (exact; mu_0=0 / mu_1=1 by construction, no refit); 'refit' fits the "
                             "double Gaussian again directly to the electron-unit histogram, letting the "
                             "peaks/widths re-optimise in electron space.")
    parser.add_argument("--show_titles", action="store_true",
                        default=_config_default(config, 'show_titles', True),
                        help="Show suptitles and per-axis titles on all plots")
    parser.add_argument("--no-show_titles", dest="show_titles", action="store_false",
                        help="Hide all suptitles and per-axis titles")
    parser.add_argument("--plot_zero_one_sharex", action="store_true",
                        default=_config_default(config, 'plot_zero_one_sharex', True),
                        help="Share x-axis range across the 2x2 zero/one peak subplots")
    parser.add_argument("--no-plot_zero_one_sharex", dest="plot_zero_one_sharex", action="store_false",
                        help="Allow independent x-axis ranges per zero/one peak subplot")
    parser.add_argument("--plot_all_peaks_sharex", action="store_true",
                        default=_config_default(config, 'plot_all_peaks_sharex', True),
                        help="Share x-axis range across the 2x2 all-peaks subplots")
    parser.add_argument("--no-plot_all_peaks_sharex", dest="plot_all_peaks_sharex", action="store_false",
                        help="Allow independent x-axis ranges per all-peaks subplot")
    parser.add_argument("--plot_nonlinearity_sharex", action="store_true",
                        default=_config_default(config, 'plot_nonlinearity_sharex', True),
                        help="Share x-axis range across the 2x2 nonlinearity subplots")
    parser.add_argument("--no-plot_nonlinearity_sharex", dest="plot_nonlinearity_sharex", action="store_false",
                        help="Allow independent x-axis ranges per nonlinearity subplot")
    parser.add_argument("--plot_zero_one_sharey", action="store_true",
                        default=_config_default(config, 'plot_zero_one_sharey', True),
                        help="Share y-axis range across the 2x2 zero/one peak subplots")
    parser.add_argument("--no-plot_zero_one_sharey", dest="plot_zero_one_sharey", action="store_false",
                        help="Allow independent y-axis ranges per zero/one peak subplot")
    parser.add_argument("--plot_all_peaks_sharey", action="store_true",
                        default=_config_default(config, 'plot_all_peaks_sharey', True),
                        help="Share y-axis range across the 2x2 all-peaks subplots")
    parser.add_argument("--no-plot_all_peaks_sharey", dest="plot_all_peaks_sharey", action="store_false",
                        help="Allow independent y-axis ranges per all-peaks subplot")
    parser.add_argument("--peak_number_labels_individual", action="store_true",
                        default=_config_default(config, 'peak_number_labels_individual', True),
                        help="Annotate each detected peak with its index (0, 1, 2, ...) on the individual all-peaks plots")
    parser.add_argument("--no-peak_number_labels_individual", dest="peak_number_labels_individual", action="store_false",
                        help="Hide per-peak index annotations on the individual all-peaks plots")
    parser.add_argument("--peak_number_labels_together", action="store_true",
                        default=_config_default(config, 'peak_number_labels_together', True),
                        help="Annotate each detected peak with its index (0, 1, 2, ...) on the 2x2 together all-peaks subplot")
    parser.add_argument("--no-peak_number_labels_together", dest="peak_number_labels_together", action="store_false",
                        help="Hide per-peak index annotations on the 2x2 together all-peaks subplot")
    parser.add_argument("--peak_number_label_size", type=float,
                        default=_config_default(config, 'peak_number_label_size', 8),
                        help="Font size for the per-peak index annotations on the all-peaks plot")
    parser.add_argument("--plot_nonlinearity_sharey", action="store_true",
                        default=_config_default(config, 'plot_nonlinearity_sharey', True),
                        help="Share y-axis range across the 2x2 nonlinearity subplots")
    parser.add_argument("--no-plot_nonlinearity_sharey", dest="plot_nonlinearity_sharey", action="store_false",
                        help="Allow independent y-axis ranges per nonlinearity subplot")

    parsed_args = parser.parse_args(args)

    if parsed_args.file_string is None:
        parser.error('file_string is required unless it is provided in the JSON config')

    # Fall back to the config keyword (or 'default') when no CLI zero/one limits were
    # given. Kept out of argparse's `type` because the 'default'/'none' string keywords
    # would fail the float coercion applied to a string default.
    if parsed_args.plot_zero_one_xlim is None:
        parsed_args.plot_zero_one_xlim = _config_default(config, 'plot_zero_one_xlim', 'default')
    if parsed_args.plot_zero_one_ylim is None:
        parsed_args.plot_zero_one_ylim = _config_default(config, 'plot_zero_one_ylim', 'default')
    if parsed_args.plot_zero_one_ylim_electrons is None:
        parsed_args.plot_zero_one_ylim_electrons = _config_default(config, 'plot_zero_one_ylim_electrons', 'default')

    # argparse's `choices` validates CLI strings but not config-supplied defaults, so
    # re-check electron_fit_mode here to also reject bad values coming from the JSON.
    if parsed_args.electron_fit_mode not in (None, 'transform', 'refit'):
        parser.error(f"electron_fit_mode must be 'transform' or 'refit' (got {parsed_args.electron_fit_mode!r})")

    # argparse's `type` validates CLI strings but not config-supplied defaults, so
    # re-check the window scales here to also reject values < 1 coming from the JSON.
    for _scale_arg in ('zero_one_window_left_scale', 'zero_one_window_right_scale'):
        if getattr(parsed_args, _scale_arg) < 1.0:
            parser.error(f"{_scale_arg} must be >= 1.0 (got {getattr(parsed_args, _scale_arg)})")
    if int(parsed_args.zero_one_n_bins) < 10:
        parser.error(f"zero_one_n_bins must be an integer >= 10 (got {parsed_args.zero_one_n_bins})")
    parsed_args.zero_one_n_bins = int(parsed_args.zero_one_n_bins)
    if float(parsed_args.zero_one_peakfind_density) < 1.0:
        parser.error(f"zero_one_peakfind_density must be >= 1 (got {parsed_args.zero_one_peakfind_density})")
    parsed_args.zero_one_peakfind_density = float(parsed_args.zero_one_peakfind_density)

    # Group a flat CLI int list into column pairs (and validate the count) up front, so
    # the stored value is the canonical pair / per-extension form.
    try:
        parsed_args.fit_cols = _coerce_fit_cols(parsed_args.fit_cols)
    except ValueError as e:
        parser.error(str(e))

    # Collapse a single gain seed to a scalar (per-extension list otherwise) and check
    # every supplied value is positive -- the gain (ADU/e-) is always > 0. argparse's
    # `type` only validates CLI floats, so this per-entry check also covers config values.
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
