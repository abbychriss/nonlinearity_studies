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
from pathlib import Path
import sys

# Handle imports for both direct execution and module import
if __name__ == "__main__":
    # When run as script, add parent directory to path
    sys.path.insert(0, str(Path(__file__).parent))
    from stitch_fits import stitch_fits
    from nonlinearity_studies import (
        get_fits, 
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

    do_stitch_images = args.stitch_fits
    do_plot_zero_one_peaks = args.plot_zero_one_peaks
    do_plot_all_peaks = args.plot_all_peaks
    get_nonlinearity_at_charges = args.get_nonlinearity_at
    # Unpack single value from list for cleaner interface
    if get_nonlinearity_at_charges is not None and len(get_nonlinearity_at_charges) == 1:
        get_nonlinearity_at_charges = get_nonlinearity_at_charges[0]
    do_get_nonlinearity_at = get_nonlinearity_at_charges is not None
    do_plot_nonlinearity = args.plot_nonlinearity
    save_plots = args.save_plots

    if do_stitch_images:
        # Split input path to get directory and image pattern
        image_dir, image_pattern = _derive_data_path(args.file_string)
        
        # Check if path already contains 'combined-fits'
        if 'combined-fits' in image_dir.parts:
            # If combined-fits is already in the path, find its parent as data_path
            combined_fits_idx = image_dir.parts.index('combined-fits')
            data_path = Path(*image_dir.parts[:combined_fits_idx])
            # Don't need to stitch - files are already there
            image_name = next(Path(image_dir).glob(image_pattern)).name
        else:
            # Raw images case - stitch them
            data_path = image_dir.parent
            stitched_file = stitch_fits(data_path, directory=image_dir.name, image=image_pattern, 
                                         out_path='combined-fits/', print_header=False)
            image_name = Path(stitched_file).name
    else:
        file_path = Path(args.file_string)
        
        # Check if path contains 'combined-fits' to find the base data_path
        if 'combined-fits' in file_path.parts:
            combined_fits_idx = file_path.parts.index('combined-fits')
            data_path = Path(*file_path.parts[:combined_fits_idx])
        else:
            data_path = file_path.parent
        
        image_name = file_path.name

    # Derive fig_path from data_path
    fig_path = data_path.parent / 'plots'
    fig_path.mkdir(parents=True, exist_ok=True)
    
    print(f'Analyzing image: {image_name}')
    
    # Get data from fits file
    if do_stitch_images:
        # In stitch mode, files are in combined-fits subdirectory
        fits_path = str(data_path / 'combined-fits' / image_name)
    else:
        # In non-stitch mode, use the input path directly
        fits_path = args.file_string
    
    data_ext = get_fits(fits_path)

    # Fit zeroth and first electron peaks to double gaussians
    zero_one_counts_ext, zero_one_edges_ext, pedestals, gains, \
    double_gauss_popts, zero_one_ranges = get_zero_one_peaks_ext(data_ext, n=100, fit_bounds='default')

    # Apply scipy peak finder to find location of every electron peak
    counts_ext, edges_ext, peaks_ext, centers_ext, hist_ranges = get_all_peaks_ext(data_ext, 
                                                                                widths=[0.3,0.3,0.6,0.3], 
                                                                                buffers=[2.1,2.1,2.1,2.1], 
                                                                                pedestals=pedestals, 
                                                                                double_gauss_popts=double_gauss_popts,
                                                                                gains=gains,
                                                                                bins='default',
                                                                                flatten=True,
                                                                                do_convert_to_electrons=True, 
                                                                                range_left='default', 
                                                                                range_right=2500, 
                                                                                bin_factor=8)

    # Fit parabola to nonlinearity curve
    fit_range_right_ext= [600,800,500,1000]

    peak_charge_e_ext, charge_minus_npeak_ext, parabola_coeffs, parabola_pcovs, nonlinearity_at_500 = get_nonlinearity_ext(peaks_ext,
                                                                                                                           centers_ext, 
                                                                                                                           pedestals, 
                                                                                                                           gains, 
                                                                                                                           fit_range_right_ext, 
                                                                                                                           do_convert_to_electrons=False,
                                                                                                                           fit_bounds_low=-100, 
                                                                                                                           fit_bounds_high=100,
                                                                                                                           print_values=True)

    # Get nonlinearity at specified charge value(s)
    if do_get_nonlinearity_at:
        get_nonlinearity_at_ext(get_nonlinearity_at_charges, parabola_coeffs, parabola_pcovs, fit_range_right_ext, print_values=True)

    # Fit a double gaussian to zero + 1 electron peak in each extension
    if do_plot_zero_one_peaks:
        plot_zero_one_peaks(data_ext, 
                            zero_one_counts_ext,
                            zero_one_edges_ext, 
                            pedestals, 
                            gains, 
                            double_gauss_popts, 
                            zero_one_ranges,
                            individual_figsize=(6,5), 
                            subplots_figsize=(9,7),
                            xlim='default',
                            #ylim=(0.00001,2e5),
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
        plot_all_peaks_range_left=np.min(np.array(hist_ranges).flatten())
        plot_all_peaks_range_right=np.max(np.array(hist_ranges).flatten())

        plot_all_peaks(counts_ext, 
                    peaks_ext, 
                    centers_ext,  
                    xlim=(plot_all_peaks_range_left, plot_all_peaks_range_right),
                    ylim='none', 
                    yscale='linear',
                    plot_individual=False, 
                    plot_together=True, 
                    draw_lines=True, 
                    linecolor='r', 
                    linestyle='--',
                    individual_figsize=(7,6), 
                    subplots_figsize=(9,7),
                    suptitle='Peaks in Pixel Charge Distribution',
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
                        suptitle='Pixel Charge Nonlinearity Curve (Nimages = 10)',
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
        

def init_argparse():
    """
    Initializes the ArgumentParser object and defines arguments.
    """
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

    parser.add_argument('file_string', type=str, 
                       help='Absolute or relative path to image file (.fz or .fits accepted)')
    parser.add_argument("-f", "--stitch_fits", action="store_true", default=False, 
                       help="Stitch FITS files by extension")
    parser.add_argument("-z", "--plot_zero_one_peaks", action="store_true", default=False, 
                       help="Plot fits to zero/one electron peaks")
    parser.add_argument("-a", "--plot_all_peaks", action="store_true", default=False, 
                       help="Plot entire charge distribution with line at each peak")
    parser.add_argument("-g", "--get_nonlinearity_at", nargs='+', type=float, default=None,
                       help="Estimate nonlinearity at specified charge value(s) using parabolic fit")
    parser.add_argument("-n", "--plot_nonlinearity", action="store_true", default=False, 
                       help="Plot nonlinearity curve with quadratic fit")
    parser.add_argument("-s", "--save_plots", action="store_true", default=False, 
                       help="Save all plots as jpeg images")

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = init_argparse()
    main(args)