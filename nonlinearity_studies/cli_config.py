"""
cli_config — CLI config loading & argument-normalization helpers for run-nonlinearity-studies, split out of run_nonlinearity_studies.py.
"""
import argparse
import glob
import hashlib
import json
from pathlib import Path

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
    'plot_charge_per_column_together',
    'plot_charge_per_column_individual',
    'plot_charge_per_column_figsize',
    'plot_charge_per_column_individual_figsize',
    'plot_zero_one_ADU_individual',
    'plot_zero_one_ADU_together',
    'plot_zero_one_electrons_individual',
    'plot_zero_one_electrons_together',
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
    'save_output',
    'output_dir',
    'show_plots',
    'verbose',
    'nimages',
    'extra_plot_title',
    'peak_finder_widths',
    'peak_finder_buffers',
    'peak_finder_prominences',
    'nonlinearity_peakfinder_bin_factor',
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
    'fit_cols',
    'zero_one_gain_guess',
    'plot_zero_one_individual_figsize',
    'plot_zero_one_subplots_figsize',
    'plot_all_peaks_xlim',
    'plot_all_peaks_ylim',
    'plot_all_peaks_yscale',
    'plot_all_peaks_window_bin_factor',
    'plot_all_peaks_individual_figsize',
    'plot_all_peaks_subplots_figsize',
    'plot_nonlinearity_individual_figsize',
    'plot_nonlinearity_subplots_figsize',
    'plot_nonlinearity_xlim',
    'plot_nonlinearity_ylim',
    'plot_all_peaks_individual',
    'plot_all_peaks_together',
    'plot_nonlinearity_individual',
    'plot_nonlinearity_together',
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


_RUN_HASH_EXCLUDE = {
    'config', 'json', 'verbose', 'save_output', 'output_dir', 'show_plots',
    'pedsub_cache_dir', 'force_pedsub',
    'plot_zero_one_ADU_individual', 'plot_zero_one_ADU_together',
    'plot_zero_one_electrons_individual', 'plot_zero_one_electrons_together',
    'plot_all_peaks_individual', 'plot_all_peaks_together',
    'plot_nonlinearity_individual', 'plot_nonlinearity_together',
    'plot_charge_per_column_together', 'plot_charge_per_column_individual',
    'plot_resolution_individual', 'plot_resolution_together',
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
