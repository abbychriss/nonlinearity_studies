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
import json
from pathlib import Path
import sys
import re

# Handle imports for both direct execution and module import
if __name__ == "__main__":
    # When run as script, add parent directory to path
    sys.path.insert(0, str(Path(__file__).parent))
    from stitch_fits import stitch_fits
    from nonlinearity_studies import (
        get_fits,
        pedestal_subtract,
        get_zero_one_peaks_ext, 
        get_all_peaks_ext, 
        get_nonlinearity_ext, 
        get_nonlinearity_at_ext,
        plot_zero_one_peaks, 
        plot_all_peaks, 
        plot_nonlinearity
    )
else:
    # When imported as module, use relative imports
    from .stitch_fits import stitch_fits
    from .nonlinearity_studies import (
        get_fits, 
        pedestal_subtract,
        get_zero_one_peaks_ext, 
        get_all_peaks_ext, 
        get_nonlinearity_ext, 
        get_nonlinearity_at_ext,
        plot_zero_one_peaks, 
        plot_all_peaks, 
        plot_nonlinearity
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


CONFIG_KEYS = {
    'file_string',
    'stitch_fits',
    'plot_zero_one',
    'plot_all_peaks',
    'get_nonlinearity_at',
    'plot_nonlinearity',
    'save_plots',
    'output_dir',
    'verbose',
    'nimages',
    'extra_plot_title',
    'peak_finder_widths',
    'peak_finder_buffers',
    'fit_range_right',
    'max_one_peak_sigma_ratio',
    'do_pedestal_subtraction',
    'n_std_to_mask',
    'pedestal_subtraction_axis',
    'use_biweight_loc',
    'use_biweight_midvar'
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


def _normalize_scalar_or_list(value):
    if isinstance(value, list) and len(value) == 1:
        return value[0]
    return value


def main(args=None):
    """
    The main executable function of the script.
    
    Args:
        args: The Namespace object containing the parsed command-line arguments.
              If None, arguments are parsed from command line.
    """
    if args is None:
        args = init_argparse()
    
    file_path = Path(args.file_string)

    if not args.stitch_fits:
        if not file_path.is_absolute():
            # For relative paths, find them from current directory
            if not file_path.exists():
                # Try common search patterns
                for search_path in Path('.').rglob(file_path.name):
                    file_path = search_path
                    break

    # Get values from argparse arguments
    do_stitch_images = args.stitch_fits
    do_plot_zero_one_peaks = args.plot_zero_one
    do_plot_all_peaks = args.plot_all_peaks
    get_nonlinearity_at_charges = args.get_nonlinearity_at

    # Unpack single value from list for cleaner interface
    if get_nonlinearity_at_charges is not None:
        get_nonlinearity_at_charges = _normalize_scalar_or_list(get_nonlinearity_at_charges)
    
    # Get right bound for parabolic fit from argparse (can be single value or list of values for each extension)
    fit_range_right_ext = _normalize_scalar_or_list(args.fit_range_right)

    do_get_nonlinearity_at = get_nonlinearity_at_charges is not None
    do_plot_nonlinearity = args.plot_nonlinearity
    save_plots = args.save_plots
    output_dir = args.output_dir
    verbose = args.verbose
    max_one_peak_sigma_ratio = args.max_one_peak_sigma_ratio

    stitched_dir_name = 'stitched-fits'

    if do_stitch_images:
        input_dir, input_pattern = _derive_data_path(args.file_string)

        if stitched_dir_name in input_dir.parts:
            fits_file_path = next(input_dir.glob(input_pattern), None)
            if fits_file_path is None:
                print('\nError: no files found matching the specified stitched FITS pattern.')
                sys.exit(1)
        else:
            stitched_file = stitch_fits(
                input_dir.parent,
                directory=input_dir.name,
                image=input_pattern,
                out_path=Path(input_dir.name) / stitched_dir_name,
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
        default_fig_path = fits_file_path / 'plots'

    if output_dir is not None:
        fig_path = Path(output_dir)
    else:
        fig_path = default_fig_path

    if save_plots:
        fig_path.mkdir(parents=True, exist_ok=True)

    image_name = fits_file_path.name

    fits_path = str(fits_file_path)
    print(f'Analyzing image: {fits_path}\n')
    
    # Load data from FITS file
    data_ext = get_fits(fits_path)

    # If user specifies pedestal subtraction process, apply to all data before doing analysis
    # This will subtract the average baseline charge of specified axis (row/column) from that axis but keeps the data in ADU
    if args.do_pedestal_subtraction:
        pedestal_subtracted_data_ext = np.array([pedestal_subtract(data, n_std_to_mask=args.n_std_to_mask, axis=args.pedestal_subtraction_axis, use_biweight_loc=args.use_biweight_loc, use_biweight_midvar=args.use_biweight_midvar) for data in data_ext])
        data_ext = pedestal_subtracted_data_ext

    if not isinstance(fit_range_right_ext, list):
        fit_range_right_ext = [fit_range_right_ext] * len(data_ext)
    elif len(fit_range_right_ext) != len(data_ext):
        print(
            f'\nError: fit_range_right must be a single value or one value per extension '
            f'({len(data_ext)} values for this file).\n'
        )
        sys.exit(1)

    # Extract number of stitched images from filename; --nimages arg or config overrides
    nimages = None
    match = re.search(r'_(\d+)_stitched', fits_path)
    if match:
        nimages = int(match.group(1))
    if args.nimages is not None:
        nimages = args.nimages

    # Fit zeroth and first electron peaks to double gaussians
    zero_one_counts_ext, zero_one_edges_ext, pedestals, gains, double_gauss_popts, zero_one_ranges = get_zero_one_peaks_ext(
        data_ext,
        n=100,
        fit_bounds='default',
        max_one_peak_sigma_ratio=max_one_peak_sigma_ratio,
    )

    # Apply scipy peak finder to find location of every electron peak
    counts_ext, edges_ext, peaks_ext, centers_ext, hist_ranges = get_all_peaks_ext(data_ext,
                                                                                widths=args.peak_finder_widths,
                                                                                buffers=args.peak_finder_buffers,
                                                                                pedestals=pedestals, 
                                                                                double_gauss_popts=double_gauss_popts,
                                                                                gains=gains,
                                                                                bins='default',
                                                                                flatten=True,
                                                                                do_convert_to_electrons=True,
                                                                                range_left='default', 
                                                                                range_right=1500, 
                                                                                bin_factor=10,
                                                                                print_values=verbose)

    # Fit parabola to nonlinearity curve
    peak_charge_e_ext, charge_minus_npeak_ext, parabola_coeffs, parabola_pcovs, = get_nonlinearity_ext(peaks_ext,
                                                                                                        centers_ext, 
                                                                                                        pedestals, 
                                                                                                        gains, 
                                                                                                        fit_range_right_ext, 
                                                                                                        do_convert_to_electrons=False,
                                                                                                        fit_bounds_low=-100, 
                                                                                                        fit_bounds_high=100)

    # Get nonlinearity at specified charge value(s)
    if do_get_nonlinearity_at:
        get_nonlinearity_at_ext(get_nonlinearity_at_charges, parabola_coeffs, parabola_pcovs, fit_range_right_ext)

    # Fit a double gaussian to zero + 1 electron peak in each extension
    if do_plot_zero_one_peaks:
        plot_zero_one_peaks(data_ext, 
                            zero_one_counts_ext,
                            zero_one_edges_ext, 
                            pedestals, 
                            gains, 
                            double_gauss_popts, 
                            zero_one_ranges,
                            individual_figsize=(7,5), 
                            subplots_figsize=(10,8),
                            xlim='default',
                            ylim='default',
                            additional_title=f'{args.extra_plot_title}' if args.extra_plot_title else '',
                            suptitle='Double-Gaussian Fit to Zero and One Electron Peaks',
                            nimages=nimages,
                            yscale='linear',
                            fontsize=8,
                            n=100, 
                            do_convert_to_electrons=True,
                            plot_individual=False,
                            plot_together=True,
                            save_plots=save_plots,
                            fig_path=str(fig_path), 
                            file=image_name, 
                            dpi=350)

    if do_plot_all_peaks:
        plot_all_peaks_range_left= 1000#np.min(np.array(hist_ranges).flatten())
        plot_all_peaks_range_right= 1010#np.max(np.array(hist_ranges).flatten())

        plot_all_peaks(counts_ext, 
                    peaks_ext, 
                    centers_ext,  
                    xlim=(plot_all_peaks_range_left, plot_all_peaks_range_right),
                    ylim= (0,30),#'none', 
                    yscale='linear',
                    plot_individual=False, 
                    plot_together=True, 
                    draw_lines=True, 
                    linecolor='r', 
                    linestyle='--',
                    individual_figsize=(7,6), 
                    additional_title=args.extra_plot_title,
                    suptitle='Peaks in Pixel Charge Distribution',
                    nimages=nimages,
                    save_plots=save_plots,
                    fig_path=str(fig_path),
                    file=image_name, 
                    dpi=350,)

    if do_plot_nonlinearity:
        plot_nonlinearity(peaks_ext,
                        parabola_coeffs, 
                        peak_charge_e_ext, 
                        charge_minus_npeak_ext,
                        fit_range_right_ext,
                        xlim='default', 
                        ylim='default',
                        individual_figsize=(6,5), 
                        subplots_figsize=(9,7),
                        additional_title=args.extra_plot_title,
                        suptitle='Pixel Charge Nonlinearity Curve',
                        nimages=nimages,
                        line_color='r', 
                        scatter_color='b', 
                        s=2, 
                        alpha=1,
                        plot_individual=False, 
                        plot_together=True, 
                        save_plots=save_plots, 
                        fig_path=str(fig_path), 
                        file=image_name, 
                        dpi=350)
        

def init_argparse(args=None):
    """
    Initializes the ArgumentParser object and defines arguments.
    """
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument("-c", "--config", type=str, default=None,
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

    parser.add_argument("-c", "--config", type=str, default=pre_args.config,
                        help="Path to JSON file containing command-line arguments. Explicit CLI arguments override JSON values.")
    parser.add_argument('file_string', type=str, nargs='?',
                       default=_config_default(config, 'file_string', None),
                       help='Absolute or relative path to image file (.fz or .fits accepted) if not stitching (ex: data/avg_img.fz), or to directory with images if stitching (can include glob pattern like "*", ex: data/03-12-2026/avg*.fz or data/*). May also be set in JSON config.')
    parser.add_argument("-f", "--stitch_fits", action="store_true",
                       default=_config_default(config, 'stitch_fits', False),
                       help="Stitch FITS files by extension")
    parser.add_argument("--no-stitch_fits", dest="stitch_fits", action="store_false",
                       help="Disable FITS stitching when enabled by JSON config")
    parser.add_argument("-z", "--plot_zero_one", action="store_true",
                       default=_config_default(config, 'plot_zero_one', False),
                       help="Plot fits to zero/one electron peaks")
    parser.add_argument("--no-plot_zero_one", dest="plot_zero_one", action="store_false",
                       help="Disable zero/one peak plotting when enabled by JSON config")
    parser.add_argument("-a", "--plot_all_peaks", action="store_true",
                       default=_config_default(config, 'plot_all_peaks', False),
                       help="Plot entire charge distribution with line at each peak")
    parser.add_argument("--no-plot_all_peaks", dest="plot_all_peaks", action="store_false",
                       help="Disable all-peaks plotting when enabled by JSON config")
    parser.add_argument("-g", "--get_nonlinearity_at", nargs='+', type=float,
                       default=_config_default(config, 'get_nonlinearity_at', None),
                       help="Estimate nonlinearity at specified charge value(s) using parabolic fit")
    parser.add_argument("-n", "--plot_nonlinearity", action="store_true",
                       default=_config_default(config, 'plot_nonlinearity', False),
                       help="Plot nonlinearity curve with quadratic fit")
    parser.add_argument("--no-plot_nonlinearity", dest="plot_nonlinearity", action="store_false",
                       help="Disable nonlinearity plotting when enabled by JSON config")
    parser.add_argument("-s", "--save_plots", action="store_true",
                       default=_config_default(config, 'save_plots', False),
                       help="Save all plots as jpeg images")
    parser.add_argument("--no-save_plots", dest="save_plots", action="store_false",
                       help="Disable plot saving when enabled by JSON config")
    parser.add_argument("-o", "--output_dir", type=str,
                       default=_config_default(config, 'output_dir', None),
                       help="Directory to save all plots")
    parser.add_argument("-v", "--verbose", action="store_true",
                       default=_config_default(config, 'verbose', False),
                       help="Print verbose output")
    parser.add_argument("--no-verbose", dest="verbose", action="store_false",
                       help="Disable verbose output when enabled by JSON config")
    parser.add_argument("--extra_plot_title", type=str,
                       default=_config_default(config, 'extra_plot_title', ''),
                        help="Additional title added at beginning of all plot titles ('<Additional title>+<Default title>)")
    parser.add_argument("--nimages", type=int,
                       default=_config_default(config, 'nimages', None),
                        help="Number of images used in the stack. If not provided, extracted from filename for stitched files.")
    parser.add_argument("--peak_finder_widths", nargs='+', type=float,
                       default=_config_default(config, 'peak_finder_widths', [0.9, 0.9, 0.9, 0.9]),
                        help="Widths for peak finder. Represents maximum width around each peak in e-.")
    parser.add_argument("--peak_finder_buffers", nargs='+', type=int,
                       default=_config_default(config, 'peak_finder_buffers', [3, 3, 3, 3]),
                        help="Buffers for peak finder. Used to calculate minimum distance between neighboring peaks by d = nbins_per_electron - buffer. Units are in number of bins.")
    parser.add_argument("--fit_range_right", nargs='+', type=int,
                       default=_config_default(config, 'fit_range_right', 500),
                        help="Charge value in electrons to fit nonlinearity curve up to, can be list of 4 values (one per extension) or integer (used for all extensions)")
    parser.add_argument("--max_one_peak_sigma_ratio", type=float,
                        default=_config_default(config, 'max_one_peak_sigma_ratio', 1.5),
                        help="Maximum allowed one-electron width relative to the inferred zero-peak width")
    parser.add_argument("--do_pedestal_subtraction", action="store_true",
                        default=_config_default(config, 'do_pedestal_subtraction', True),
                        help="Whether to do pedestal subtraction along one or two axes as preprocess step")
    parser.add_argument("--n_std_to_mask", type=float,
                        default=_config_default(config, 'n_std_to_mask', 1.5),
                        help="Number of standard deviations from mean charge along axis to mask for pedestal subtraction")
    parser.add_argument("--pedestal_subtraction_axis", type=str,
                        default=_config_default(config, 'pedestal_subtraction_axis', 'row'),
                        help="Which axis to compute pedestals across. Options: 'row', 'col', 'row_then_col','col_then_row'")
    parser.add_argument("--use_biweight_loc", action="store_true",
                        default=_config_default(config, 'use_biweight_loc', True),
                        help="If true, uses Tukey biweight location (more robust - iteratively gets rid of outliers). If false, uses simple average.")
    parser.add_argument("--use_biweight_midvar", action="store_true",
                        default=_config_default(config, 'use_biweight_midvar', True),
                        help="If true, uses Tukey biweight midvariance. If false, uses standard deviation.")
    

    parsed_args = parser.parse_args(args)

    if parsed_args.file_string is None:
        parser.error('file_string is required unless it is provided in the JSON config')

    return parsed_args


if __name__ == '__main__':
    args = init_argparse()
    main(args)
