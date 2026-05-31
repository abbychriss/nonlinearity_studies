import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import biweight_location, biweight_midvariance
from scipy.optimize import curve_fit
from scipy.signal import find_peaks as scipy_find_peaks
import math

from pathlib import Path
from glob import glob

#---------------- ANALYSIS FUNCTIONS ----------------------------

#---------------- (0) Convert to electrons ----------------------
def convert_to_electrons(data, pedestal, gain, flatten=True):
    if flatten:
        data = np.array(data).flatten()
    data_electrons = (data - pedestal) / gain  # Subtract pedestal (mean ADU of zero electron peak) and divide by gain
    return data_electrons

#---------------- (1) Calculate noise/gain ----------------------
def _smooth_counts(counts, window=5):
    if len(counts) < 3 or window <= 1:
        return counts.astype(float)

    window = min(window, len(counts))
    if window % 2 == 0:
        window -= 1
    if window < 3:
        return counts.astype(float)

    kernel = np.ones(window) / window
    return np.convolve(counts, kernel, mode='same')


def _make_histogram(data, hist_range, n, max_bins=4000):
    left, right = hist_range
    if not np.isfinite(left) or not np.isfinite(right) or right <= left:
        raise ValueError(f'Invalid histogram range: {hist_range}')

    nbins = min(max_bins, max(50, int(n * (right - left))))
    data_window = data[(data > left) & (data < right)]
    counts, edges = np.histogram(data_window, bins=nbins, range=(left, right))
    centers = 0.5 * (edges[:-1] + edges[1:])

    if len(centers) < 2:
        raise ValueError(f'Could not make a useful histogram for range {hist_range}')

    return data_window, counts, edges, centers


def _estimate_peak_width(centers, counts, peak_index):
    smooth_counts = _smooth_counts(counts)
    bin_width = centers[1] - centers[0]
    peak_height = smooth_counts[peak_index]

    if peak_height <= 0:
        return bin_width

    half_max = 0.5 * peak_height

    left_index = peak_index
    while left_index > 0 and smooth_counts[left_index] > half_max:
        left_index -= 1

    right_index = peak_index
    while right_index < len(smooth_counts) - 1 and smooth_counts[right_index] > half_max:
        right_index += 1

    if right_index <= left_index:
        return bin_width

    return max((centers[right_index] - centers[left_index]) / 2.355, bin_width)


def _clip_to_bounds(values, bounds):
    low, high = [np.array(bound, dtype=float) for bound in bounds]
    values = np.array(values, dtype=float)
    eps = np.maximum(1e-12, 1e-9 * np.maximum(1, np.abs(high - low)))
    return np.minimum(np.maximum(values, low + eps), high - eps)


def _auto_zero_one_setup(data, zero_one_test_range, n, max_one_peak_sigma_ratio=1.5):
    if max_one_peak_sigma_ratio is not None and max_one_peak_sigma_ratio <= 0:
        raise ValueError('max_one_peak_sigma_ratio must be positive or None')

    use_auto_range = (
        zero_one_test_range is None
        or (isinstance(zero_one_test_range, str) and zero_one_test_range in ('auto', 'default'))
    )

    if use_auto_range:
        range_left, range_right = np.percentile(data, [0.1, 80])

        if not np.isfinite(range_left) or not np.isfinite(range_right) or range_right <= range_left:
            range_left, range_right = np.percentile(data, [0, 90])

        if not np.isfinite(range_left) or not np.isfinite(range_right) or range_right <= range_left:
            range_left = np.min(data)
            range_right = np.max(data)
    else:
        range_left, range_right = zero_one_test_range

    _, counts_test, edges_test, centers_test = _make_histogram(data, (range_left, range_right), n)

    if max(counts_test) == 0:
        raise ValueError(
            f'No data found in zero-one test range {(range_left, range_right)}. '
            'Use zero_one_test_range="auto" or choose a range around the pedestal.'
        )

    smooth_test = _smooth_counts(counts_test)
    zero_peak_index = np.argmax(smooth_test)
    zero_peak_charge = centers_test[zero_peak_index]
    zero_peak_width = _estimate_peak_width(centers_test, counts_test, zero_peak_index)
    bin_width = centers_test[1] - centers_test[0]
    zero_peak_width = max(zero_peak_width, bin_width)

    search_left = zero_peak_charge - 5 * zero_peak_width
    search_right = zero_peak_charge + max(20 * zero_peak_width, 2 * (range_right - zero_peak_charge))

    data_high = np.percentile(data, 99)
    if np.isfinite(data_high):
        search_right = min(search_right, data_high)

    if search_right <= zero_peak_charge:
        search_right = zero_peak_charge + 20 * zero_peak_width

    _, search_counts, search_edges, search_centers = _make_histogram(data, (search_left, search_right), n)
    smooth_search = _smooth_counts(search_counts, window=9)
    search_bin_width = search_centers[1] - search_centers[0]
    peak_distance = max(1, int(2 * zero_peak_width / search_bin_width))
    peak_prominence = max(5, 0.002 * max(smooth_search))
    peak_indices, _ = scipy_find_peaks(
        smooth_search,
        prominence=peak_prominence,
        distance=peak_distance,
    )

    one_peak_min_charge = zero_peak_charge + max(zero_peak_width, 0.2)
    right_peak_indices = [
        peak_index for peak_index in peak_indices
        if search_centers[peak_index] > one_peak_min_charge
    ]

    if right_peak_indices:
        found_one_peak = True
        one_peak_index = max(right_peak_indices, key=lambda peak_index: smooth_search[peak_index])
        one_peak_charge = search_centers[one_peak_index]
        one_peak_height = smooth_search[one_peak_index]
    else:
        found_one_peak = False
        one_peak_charge = zero_peak_charge + 4 * zero_peak_width
        one_peak_height = 0.1 * smooth_search[np.argmin(np.abs(search_centers - zero_peak_charge))]

    gain_guess = max(one_peak_charge - zero_peak_charge, 4 * zero_peak_width)

    zero_one_left = zero_peak_charge - max(4 * zero_peak_width, 0.5 * gain_guess)
    zero_one_right = zero_peak_charge + max(1.8 * gain_guess, 8 * zero_peak_width)
    zero_one_range = [zero_one_left, zero_one_right]

    _, zero_one_counts, zero_one_edges, zero_one_centers = _make_histogram(data, zero_one_range, n, max_bins=2000)
    max_zero_one_counts = max(zero_one_counts)

    if max_zero_one_counts == 0:
        raise ValueError(f'No data found in inferred zero-one range {zero_one_range}')

    fit_left, fit_right = zero_one_range
    m0_margin = max(2 * zero_peak_width, 0.3 * gain_guess)
    if found_one_peak:
        m1_margin = max(0.02 * zero_peak_width, 0.5 * (zero_one_centers[1] - zero_one_centers[0]))
        m1_low = max(
            zero_peak_charge + max(zero_peak_width, 0.2 * gain_guess),
            one_peak_charge - m1_margin,
        )
        m1_high = min(
            fit_right,
            zero_peak_charge + max(1.7 * gain_guess, 8 * zero_peak_width),
            one_peak_charge + m1_margin,
        )
    else:
        m1_low = zero_peak_charge + max(zero_peak_width, 0.2 * gain_guess)
        m1_high = min(fit_right, zero_peak_charge + max(1.7 * gain_guess, 8 * zero_peak_width))

    if m1_high <= m1_low:
        m1_high = fit_right

    fit_bounds_low = [
        max((zero_one_centers[1] - zero_one_centers[0]) / 10, 1e-8),
        max(fit_left, zero_peak_charge - m0_margin),
        max((zero_one_centers[1] - zero_one_centers[0]) / 10, 1e-8),
        m1_low,
        0,
        0,
    ]
    fit_bounds_high = [
        max(gain_guess, 4 * zero_peak_width),
        min(fit_right, zero_peak_charge + m0_margin),
        max(1.5 * gain_guess, 6 * zero_peak_width),
        m1_high,
        2 * max_zero_one_counts,
        2 * max_zero_one_counts,
    ]

    if max_one_peak_sigma_ratio is not None:
        max_one_peak_sigma = max_one_peak_sigma_ratio * zero_peak_width
        min_one_peak_sigma = min(0.5 * zero_peak_width, 0.9 * max_one_peak_sigma)
        fit_bounds_low[2] = max(fit_bounds_low[2], min_one_peak_sigma)
        fit_bounds_high[2] = min(fit_bounds_high[2], max_one_peak_sigma)

        if fit_bounds_high[2] <= fit_bounds_low[2]:
            fit_bounds_high[2] = fit_bounds_low[2] * 1.01

    fit_bounds = (fit_bounds_low, fit_bounds_high)

    one_peak_bin = np.argmin(np.abs(zero_one_centers - one_peak_charge))
    one_peak_height = max(one_peak_height, zero_one_counts[one_peak_bin], 0.05 * max_zero_one_counts)
    p0 = [
        zero_peak_width,
        zero_peak_charge,
        max(1.5 * zero_peak_width, gain_guess / 3),
        one_peak_charge,
        max_zero_one_counts,
        one_peak_height,
    ]
    p0 = _clip_to_bounds(p0, fit_bounds)

    return zero_one_counts, zero_one_edges, p0, fit_bounds, zero_one_range


# Function finds noise and gain from input pixel charge data
# zero_one_test_range can be 'auto' or a range of charge (in ADU) to search for the zero and one electron peaks
def calculate_noise_gain(data, zero_one_test_range='auto', n=200, fit_bounds='default',
                         max_one_peak_sigma_ratio=1.5):

    data = np.array(data).flatten()
    data = data[np.isfinite(data)]

    if data.size == 0:
        raise ValueError('Input data contains no finite values')

    zero_one_counts, zero_one_edges, p0, auto_fit_bounds, zero_one_range = _auto_zero_one_setup(
        data,
        zero_one_test_range,
        n,
        max_one_peak_sigma_ratio=max_one_peak_sigma_ratio,
    )
    zero_one_centers = 0.5 * (zero_one_edges[:-1] + zero_one_edges[1:])

    if isinstance(fit_bounds, str) and fit_bounds in ('default', 'auto'):
        fit_bounds = auto_fit_bounds
    else:
        p0 = _clip_to_bounds(p0, fit_bounds)

    popt, pcov = curve_fit(double_gauss, zero_one_centers, zero_one_counts, p0=p0,
                           maxfev=20000, bounds=fit_bounds)
    
    # Extract pedestal, noise, gain, and rest of double gaussian coefficients from curve fit
    pedestal=tuple(popt)[1] # Pedestal is mean of zero electron peak
    noise=tuple(popt)[0] # Noise is standard deviation of zero electron peak 
    gain=tuple(popt)[3]-tuple(popt)[1] # Gain is difference between mean of one and zero electron peaks
    return zero_one_counts, zero_one_edges, pedestal, noise, gain, popt, zero_one_range

# Function calculates pedestal row by row in the image
# Options for axis are 'row', 'column', 'row_then_col', 'col_then_row' 
def pedestal_subtract(data, n_std_to_mask, axis='row', use_biweight_loc=True, use_biweight_midvar=True):

    data = np.array(data, dtype=float)

    def subtract_along(arr, ax):
        # ax=1 subtracts a per-row pedestal; ax=0 subtracts a per-column pedestal
        if use_biweight_midvar:
            std = np.sqrt(biweight_midvariance(arr, axis=ax))  # biweight_midvariance returns variance, not std
        else:
            std = np.std(arr, axis=ax)
        if use_biweight_loc:
            avg_charge = biweight_location(arr, axis=ax)
        else:
            avg_charge = np.mean(arr, axis=ax)

        avg_charge_b = np.expand_dims(avg_charge, axis=ax)
        std_b = np.expand_dims(std, axis=ax)
        mask = np.abs(arr - avg_charge_b) <= n_std_to_mask * std_b
        masked = np.where(mask, arr, np.nan)

        if use_biweight_loc:
            avg_pedestal = biweight_location(masked, axis=ax, ignore_nan=True)
        else:
            avg_pedestal = np.nanmean(masked, axis=ax)

        return arr - np.expand_dims(avg_pedestal, axis=ax)

    if axis == 'row':
        return subtract_along(data, ax=1)
    elif axis in ('column', 'col'):
        return subtract_along(data, ax=0)
    elif axis == 'row_then_col':
        return subtract_along(subtract_along(data, ax=1), ax=0)
    elif axis == 'col_then_row':
        return subtract_along(subtract_along(data, ax=0), ax=1)

    return data


_PEDSUB_HEADER_KEYS = ('PEDSUB_A', 'PEDSUB_N', 'PEDSUB_L', 'PEDSUB_V')


def _pedsub_cache_path(source_path, cache_dir=None):
    source = Path(source_path)
    base = cache_dir if cache_dir is not None else source.parent
    base = Path(base)
    if cache_dir is not None:
        base.mkdir(parents=True, exist_ok=True)
    return base / f'{source.stem}.pedsub.fits'


def _pedsub_header_matches(header, axis, n_std_to_mask, use_biweight_loc, use_biweight_midvar):
    if not all(k in header for k in _PEDSUB_HEADER_KEYS):
        return False
    return (
        header['PEDSUB_A'] == axis
        and float(header['PEDSUB_N']) == float(n_std_to_mask)
        and bool(header['PEDSUB_L']) == bool(use_biweight_loc)
        and bool(header['PEDSUB_V']) == bool(use_biweight_midvar)
    )


def pedestal_subtract_ext_cached(data_ext, source_path, n_std_to_mask, axis='row',
                                 use_biweight_loc=True, use_biweight_midvar=True,
                                 cache_dir=None, force=False, verbose=True):
    """Pedestal-subtract each extension, caching the result to a FITS file next to the source.

    On rerun, if the cache exists and its header params match the requested params, the cached
    arrays are loaded instead of recomputing. Pass force=True to bypass the cache.
    """
    cache_path = _pedsub_cache_path(source_path, cache_dir)

    if not force and cache_path.exists():
        with fits.open(str(cache_path)) as hdul:
            if _pedsub_header_matches(hdul[0].header, axis, n_std_to_mask,
                                       use_biweight_loc, use_biweight_midvar):
                if verbose:
                    print(f'Loading cached pedestal-subtracted data from {cache_path}')
                return [hdul[i].data.copy() for i in range(1, len(hdul))]
            elif verbose:
                print(f'Cached params at {cache_path} differ from current; recomputing.')

    if verbose:
        print('Computing pedestal subtraction...')
    pedsub_data_ext = [
        pedestal_subtract(data, n_std_to_mask=n_std_to_mask, axis=axis,
                          use_biweight_loc=use_biweight_loc, use_biweight_midvar=use_biweight_midvar)
        for data in data_ext
    ]

    primary = fits.PrimaryHDU()
    primary.header['PEDSUB_A'] = (axis, 'pedestal subtraction axis')
    primary.header['PEDSUB_N'] = (float(n_std_to_mask), 'n_std_to_mask')
    primary.header['PEDSUB_L'] = (bool(use_biweight_loc), 'use biweight location')
    primary.header['PEDSUB_V'] = (bool(use_biweight_midvar), 'use biweight midvariance')
    primary.header['SRC_FITS'] = (str(source_path)[-68:], 'source FITS file (truncated)')
    hdul_out = fits.HDUList([primary] + [fits.ImageHDU(data=arr) for arr in pedsub_data_ext])
    hdul_out.writeto(str(cache_path), overwrite=True)
    if verbose:
        print(f'Saved pedestal-subtracted cache to {cache_path}')

    return pedsub_data_ext


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

#---------------- (4) Get nonlinearity at charge value ----------------------------
# use coefficients of parabola fit to nonlinearity curve from fit_nonlinearity to estimate the nonlinearity
# at a single value or list of values
# required arguments: q (charge value that parabola equation is evaluated at, can be int, float or list), parabola_coeff (the curve_fit optimal parabola coefficients found in fit_nonlinearity)
# optional arguments: parabola_pcov, fit_range_right (used to ascertain how confident we should be in the nonlinearity value); print_values (if you want to print what the function outputs)
def get_nonlinearity_at(q, parabola_coeff, parabola_pcov=None, fit_range_right=None, print_values=False):
    a = parabola_coeff[0]
    b = parabola_coeff[1]
    c = parabola_coeff[2]

    print(f'Parabola best fit: Nonlin(q) = {np.round(a,5)}*q^2 + {np.round(b,5)}*q + {np.round(c,5)}\n')

    if isinstance(q, (float, int, np.floating, np.integer)):
        nonlinearity_at_q = np.round(a*q**2 + b*q + c, 3)
    elif isinstance(q, (list, np.ndarray)):
        nonlinearity_at_q = [np.round((a*q_i**2 + b*q_i + c), 3) for q_i in q]
    else:
        raise TypeError(f'q must be a number or list, got {type(q)}')

    if isinstance(q, (float, int, np.floating, np.integer)):
        print(f'Interpolated nonlinearity at {q} e- = {nonlinearity_at_q} e-\n')
    else:
        for i in range(len(q)):
            print(f'Interpolated nonlinearity at {q[i]} e- = {nonlinearity_at_q[i]} e-\n')
    
    if print_values:
        if fit_range_right is not None:
            print(f'Parabola fit was performed up to {fit_range_right} e-\n')
        if parabola_pcov is not None:
            print(f'Covariance of fit = {parabola_pcov}\n')

    return nonlinearity_at_q

#---------------- PLOTTING FUNCTIONS ----------------------------

#---------------- Plot zero-one peaks  -----------------------
# Usage: for plotting zero-one electron peaks from each extension on same subplot or individually by extension.
# Input data is list of 2D pixel charge arrays from all 4 extensions.
# xlim can be 'default', 'none', or tuple(left, right)
# ylim can be 'default', 'none' or tuple(bottom, top)
def plot_zero_one_peaks(data_ext, 
                        zero_one_counts_ext,
                        zero_one_edges_ext,
                        pedestals, 
                        gains, 
                        double_gauss_popts, 
                        zero_one_ranges,
                        individual_figsize=(7,5), 
                        subplots_figsize=(9,7),
                        xlim='default',
                        ylim='default',
                        suptitle='Double-Gaussian Fit to Zero and One Electron Peaks',
                        additional_title='',
                        nimages=10,
                        fontsize=7.5,
                        yscale='linear',
                        n=200, 
                        do_convert_to_electrons=False,
                        electron_fit_bounds=([0.001, -0.9, 0.1, 0.05, 0, 0], [1.0,  1,  1.0,  1.1, 1e10, 1e10]),
                        plot_individual=False,
                        plot_together=True,
                        do_plot_adu=True,
                        sharex=True,
                        sharey=True,
                        show_titles=True,
                        save_plots=False,
                        fig_path='./', file='zero_one_peaks',
                        dpi=350):

    fig_path = Path(fig_path)
    if file != 'zero_one_peaks':
        base_name = file[:-5] + '_zero_one_peaks'
    else:
        base_name = file
    fig_name = fig_path / base_name

    if plot_individual:
        for ext, data in (enumerate(data_ext) if do_plot_adu else []):
            data = np.array(data).flatten()

            zero_one_counts=zero_one_counts_ext[ext]
            zero_one_edges=zero_one_edges_ext[ext]
            pedestal=pedestals[ext]
            gain=gains[ext]
            double_gauss_popt=double_gauss_popts[ext]
            zero_one_range=zero_one_ranges[ext]

            fig, ax = plt.subplots(1, 1, figsize=individual_figsize, constrained_layout=True)
            if show_titles:
                fig.suptitle(f'{additional_title}{suptitle} in ADU (Nimages = {nimages}): EXT {ext + 1}')
            ax.set_xlabel('Charge (ADU)')
            ax.set_ylabel('N')

            double_gauss_coeff = tuple(double_gauss_popt)+(gain,)
            data_window = data[(data > zero_one_range[0]) & (data < zero_one_range[1])]
            nbins=int(n*(zero_one_range[1] - zero_one_range[0]))

            bin_width = zero_one_edges[1] - zero_one_edges[0]
            zero_one_centers = 0.5 * (zero_one_edges[:-1] + zero_one_edges[1:])

            if yscale=='log':
                zero_one_counts = np.maximum(zero_one_counts, 1) #need in order to prevent empty bars in histogram if there are any bins that have 0 counts (log(0) is infinite)
                ax.set_yscale('log')
            elif yscale!='linear':
                ax.set_yscale(yscale)
            ax.bar(zero_one_edges[:-1], zero_one_counts, edgecolor='none', linewidth=0, align='edge', width=np.diff(zero_one_edges))
            
            ax.plot(zero_one_centers, double_gauss(zero_one_centers, *double_gauss_popt), 'r',
                label=r'$\sigma_0$ = %5.3f, $\mu_0$ = %5.3f, $\sigma_1$ = %5.3f, $\mu_1$ = %5.3f,'%double_gauss_coeff[0:4]
                +'\n'+'$N_0$ = %5.3f, $N_1$ = %5.3f, gain = %5.3f ADU/$e^{–}$'%double_gauss_coeff[4:])
            ax.legend(loc="upper right", fontsize=fontsize)
            
            if xlim=='default':
                ax.set_xlim(zero_one_range[0],zero_one_range[1])
            elif xlim!='none':
                ax.set_xlim(xlim)
            
            if ylim=='default':
                if yscale=='log':
                    ax.set_ylim(0, np.max(zero_one_counts) * 1e4)
                elif yscale=='linear':
                    ax.set_ylim(0, np.max(zero_one_counts) + 2.5e4)
            elif ylim!='none':
                ax.set_ylim(ylim)

            if save_plots:
                output_path = fig_name.with_stem(fig_name.stem + f'_EXT{ext+1}').with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plots to {output_path}')
            plt.show()

        if do_convert_to_electrons:
            for ext, data in enumerate(data_ext):
                data=np.array(data).flatten()
                zero_one_range = zero_one_ranges[ext]
                pedestal = pedestals[ext]
                gain = gains[ext]

                data_window=data[(data > zero_one_range[0]) & (data < zero_one_range[1])]

                data_window_e = convert_to_electrons(data_window, pedestal, gain)
                zero_one_range_e = convert_to_electrons(zero_one_range, pedestal, gain)
                nbins = int(n * (zero_one_range_e[1] - zero_one_range_e[0]))

                zero_one_counts_e, zero_one_edges_e = np.histogram(data_window_e, bins=nbins, range=zero_one_range_e)
                zero_one_centers_e = 0.5 * (zero_one_edges_e[:-1] + zero_one_edges_e[1:])
                bin_width_e = zero_one_edges_e[1] - zero_one_edges_e[0]

                double_gauss_popt_e, _ = curve_fit(double_gauss, zero_one_centers_e, zero_one_counts_e, maxfev=2000,
                                                                     bounds=electron_fit_bounds)

                fig, ax = plt.subplots(1, 1, figsize=individual_figsize, constrained_layout=True)
                if show_titles:
                    fig.suptitle(rf'{additional_title}{suptitle} in $e^-$ (Nimages = {nimages}): EXT {ext + 1}')

                if yscale=='log':
                    zero_one_counts_e = np.maximum(zero_one_counts_e, 1)
                    ax.set_yscale('log')
                elif yscale!='linear':
                    ax.set_yscale(yscale)

                ax.bar(zero_one_edges_e[:-1], zero_one_counts_e, align='edge', edgecolor='none', linewidth=0, width=np.diff(zero_one_edges_e))
                ax.set_xlabel(r'Charge ($e^–$)')
                ax.set_ylabel('N')
                ax.plot(zero_one_centers_e, double_gauss(zero_one_centers_e, *double_gauss_popt_e), 'r',
                    label=r'$\sigma_0$ = %5.3f $e^{–}$, $\mu_0$ = %5.3f $e^{–}$, $\sigma_1$ = %5.3f $e^{–}$, $\mu_1$ = %5.3f $e^{–}$'%tuple(double_gauss_popt_e)[0:4])
                ax.legend(loc="upper right", fontsize=fontsize)

                if xlim=='default':
                    ax.set_xlim(zero_one_range_e[0], zero_one_range_e[1])
                elif xlim!='none':
                    ax.set_xlim(xlim)

                if ylim=='default':
                    if yscale=='log':
                        ax.set_ylim(0, np.max(zero_one_counts_e) * 1e4)
                    elif yscale=='linear':
                        ax.set_ylim(0, np.max(zero_one_counts_e) + 2.5e4)
                elif ylim!='none':
                    ax.set_ylim(ylim)

                if save_plots:
                    output_path = fig_name.with_stem(fig_name.stem + f'_electrons_EXT{ext+1}').with_suffix('.jpeg')
                    plt.savefig(str(output_path), dpi=dpi)
                    print(f'Saved plot to {output_path}')
                plt.show()

    if plot_together:

        if do_plot_adu:
            fig, axs = plt.subplots(2, 2, figsize=subplots_figsize, constrained_layout=True, sharex=sharex, sharey=sharey)
            if show_titles:
                fig.suptitle(f'{additional_title}{suptitle} in ADU (Nimages = {nimages})')
            axs = axs.flatten()

            for ext, data in enumerate(data_ext):
                data = np.array(data).flatten()
                zero_one_counts=zero_one_counts_ext[ext]
                zero_one_edges=zero_one_edges_ext[ext]
                pedestal=pedestals[ext]
                gain=gains[ext]
                double_gauss_popt=double_gauss_popts[ext]
                zero_one_range=zero_one_ranges[ext]

                ax = axs[ext]
                double_gauss_coeff = tuple(double_gauss_popt)+(gain,)
                data_window = data[(data > zero_one_range[0]) & (data < zero_one_range[1])]
                nbins=int(n*(zero_one_range[1]-zero_one_range[0]))

                zero_one_centers = 0.5 * (zero_one_edges[:-1] + zero_one_edges[1:])
                bin_width = zero_one_edges[1] - zero_one_edges[0]

                if yscale=='log':
                    zero_one_counts = np.maximum(zero_one_counts, 1) #need in order to prevent empty bars in histogram if there are any bins that have 0 counts
                    ax.set_yscale('log')
                elif yscale!='linear':
                    ax.set_yscale(yscale)

                ax.bar(zero_one_centers, zero_one_counts, align='center', edgecolor='none', linewidth=0, width=bin_width)

                ax.plot(zero_one_centers, double_gauss(zero_one_centers, *double_gauss_popt), 'r',
                    label=r'$\sigma_0$ = %5.3f, $\mu_0$ = %5.3f, $\sigma_1$ = %5.3f, $\mu_1$ = %5.3f,'%double_gauss_coeff[0:4]
                    +'\n'+'$N_0$ = %5.3f, $N_1$ = %5.3f, gain = %5.3f ADU/$e^{–}$'%double_gauss_coeff[4:])

                ax.set_xlabel('Charge (ADU)')
                ax.set_ylabel('N')
                if show_titles:
                    ax.set_title(f'EXT {ext + 1}')
                ax.legend(loc="upper right", fontsize=fontsize)

                if xlim=='default':
                    ax.set_xlim(zero_one_range[0],zero_one_range[1])
                elif xlim!='none':
                    ax.set_xlim(xlim)

                if ylim=='default':
                    if yscale=='log':
                        ax.set_ylim(0, np.max(zero_one_counts) * 1e4)
                    elif yscale=='linear':
                        ax.set_ylim(0, np.max(zero_one_counts) + 2.5e4)
                elif ylim!='none':
                    ax.set_ylim(ylim)

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
            plt.show()

        if do_convert_to_electrons:
            fig, axs = plt.subplots(2, 2, figsize=subplots_figsize, constrained_layout=True, sharex=sharex, sharey=sharey)
            if show_titles:
                fig.suptitle(rf'{additional_title}{suptitle} in $e^-$ (Nimages = {nimages})')
            axs = axs.flatten()

            for ext, data in enumerate(data_ext):
                ax = axs[ext]

                data = np.array(data).flatten()
                zero_one_range = zero_one_ranges[ext]
                pedestal = pedestals[ext]
                gain = gains[ext]

                data_window = data[(data > zero_one_range[0]) & (data < zero_one_range[1])]

                data_window_e = convert_to_electrons(data_window, pedestal, gain)
                zero_one_range_e = convert_to_electrons(zero_one_range, pedestal, gain) 
                nbins = int(n * (zero_one_range_e[1] - zero_one_range_e[0]))

                zero_one_counts_e, zero_one_edges_e = np.histogram(data_window_e, bins=nbins, range=zero_one_range_e)
                zero_one_centers_e = 0.5 * (zero_one_edges_e[:-1] + zero_one_edges_e[1:])
                bin_width_e = zero_one_edges_e[1] - zero_one_edges_e[0]

                zero_one_counts=zero_one_counts_ext[ext]
                zero_one_edges=zero_one_edges_ext[ext]
                
                double_gauss_popt_e, _ = curve_fit(double_gauss, zero_one_centers_e, zero_one_counts_e, maxfev=2000, bounds=electron_fit_bounds)
                
                if yscale=='log':
                    zero_one_counts_e = np.maximum(zero_one_counts_e, 1) #need in order to prevent empty bars in histogram if there are any bins that have 0 counts
                    ax.set_yscale('log')

                elif yscale!='linear':
                    ax.set_yscale(yscale)

                ax.bar(zero_one_centers_e, zero_one_counts_e, align='center', edgecolor='none', linewidth=0, width=bin_width_e)

                if show_titles:
                    ax.set_title(f'EXT {ext + 1}')
                ax.plot(zero_one_centers_e, double_gauss(zero_one_centers_e, *double_gauss_popt_e), 'r',
                    label=r'$\sigma_0$ = %5.3f $e^{–}$, $\mu_0$ = %5.3f $e^{–}$, $\sigma_1$ = %5.3f $e^{–}$, $\mu_1$ = %5.3f $e^{–}$'%tuple(double_gauss_popt_e)[0:4])
                ax.legend(loc="upper right", fontsize=fontsize)
                ax.set_xlabel(r'Charge ($e^–$)')
                ax.set_ylabel('N')

                if xlim=='default':
                    ax.set_xlim(zero_one_range_e[0], zero_one_range_e[1])
                elif xlim!='none':
                    ax.set_xlim(xlim)

                if ylim=='default':
                    if yscale=='log':
                        ax.set_ylim(0, np.max(zero_one_counts_e) * 1e4)
                    elif yscale=='linear':
                        ax.set_ylim(0, np.max(zero_one_counts_e) + 2.5e4)
                elif ylim!='none':
                    ax.set_ylim(ylim)

            for i in (0, 1):
                axs[i].set_xlabel('')
                axs[i].tick_params(labelbottom=True)
            for i in (1, 3):
                axs[i].set_ylabel('')
                axs[i].tick_params(labelleft=True)

            if save_plots:
                output_path = fig_name.with_stem(fig_name.stem + '_electrons').with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            plt.show()


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
                   individual_figsize=(6,5), subplots_figsize=(9,7),
                   additional_title='',
                   suptitle='Peaks in Pixel Charge Distribution',
                   nimages=10,
                   sharex=True,
                   sharey=True,
                   show_titles=True,
                   save_plots=False,
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
                fig.suptitle(f'{additional_title}{suptitle}: EXT {ext + 1}')
            ax.bar(centers, counts, align='center', edgecolor='none', linewidth=0, width=bin_width)
            ax.set_xlabel(r'Charge ($e^-$)')
            ax.set_ylabel('N')
            if yscale!='linear':
                ax.set_yscale(yscale)
            ax.set_xlim(xlim)
            if ylim!='none':
                ax.set_ylim(ylim)

            # draw vertical lines and labels at each peak
            if draw_lines:
                for i,p in enumerate(peaks):
                    peak_x = centers[p]
                    peak_y = counts[p]

                    ax.axvline(peak_x, linestyle=linestyle, color=linecolor)
                    ax.text(peak_x,
                        peak_y,
                        f"{i}",
                        verticalalignment='bottom',
                        horizontalalignment='center',
                        color=linecolor,
                        fontsize=10)
                    
            if save_plots:
                output_path = fig_name.with_stem(fig_name.stem + f'_EXT{ext+1}').with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            plt.show()

    if plot_together:
        fig, axs = plt.subplots(2,2,figsize=subplots_figsize,constrained_layout=True,sharex=sharex,sharey=sharey)
        axs=axs.flatten()
        if show_titles:
            fig.suptitle(f'{additional_title}{suptitle}')

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
            ax.set_xlim(xlim)
            if ylim!='none':
                ax.set_ylim(ylim)
            if show_titles:
                ax.set_title(f'EXT {ext + 1}')

            # draw vertical lines and labels at each peak
            if draw_lines:
                for i,p in enumerate(peaks):
                    peak_x = centers[p]
                    peak_y = counts[p]

                    ax.axvline(peak_x, linestyle=linestyle, color=linecolor)
                    ax.text(peak_x,
                        peak_y,
                        f"{i}",
                        verticalalignment='bottom',
                        horizontalalignment='center',
                        color=linecolor,
                        fontsize=10)

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
        plt.show()


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
                      suptitle='Pixel Charge Nonlinearity Curve',
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
                fig.suptitle(f'{additional_title}{suptitle} (Nimages = {nimages}): EXT {ext + 1}')
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
            ax.legend(loc="upper right", fontsize=8)
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
            plt.show()

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
            if show_titles:
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
        plt.show()


#---------------- UTILITY FUNCTIONS ----------------------------

#---------------- Get Fits ----------------------------
def get_fits(file_input):
    """
    Load FITS extensions from a file.

    Parameters
    ----------
    file_input : str or Path
        Path to the FITS file (absolute or relative to current working directory)

    Returns
    -------
    ext_charge : list
        List of data arrays from extensions 1–4
    """
    file_path = Path(file_input).resolve()
    
    # Check if file exists
    if not file_path.exists():
        raise FileNotFoundError(f"FITS file not found: {file_path}")

    # Load FITS file
    with fits.open(str(file_path)) as hdu_list:
        ext_charge = [hdu_list[i].data for i in range(1, 5)]

    return ext_charge


#---------------- Return data for each extensions in a list from pixel charge data for all extensions
def _value_for_extension(value, ext, n_ext):
    if isinstance(value, (list, tuple)) and len(value) == n_ext:
        if all(isinstance(v, (list, tuple, np.ndarray)) and len(v) == 2 for v in value):
            return value[ext]
    return value


def get_zero_one_peaks_ext(data_ext,
                           n=200, fit_bounds='default', zero_one_test_range='auto',
                           max_one_peak_sigma_ratio=1.5):
    zero_one_counts_ext = []
    zero_one_edges_ext = []
    pedestals = []
    gains = []
    double_gauss_popts = []
    zero_one_ranges = []
    for ext, data in enumerate(data_ext):

        zero_one_test_range_ext = _value_for_extension(zero_one_test_range, ext, len(data_ext))

        zero_one_counts, zero_one_edges, pedestal, noise, gain, double_gauss_popt, zero_one_range = calculate_noise_gain(
            data,
            zero_one_test_range=zero_one_test_range_ext,
            n=n,
            fit_bounds=fit_bounds,
            max_one_peak_sigma_ratio=max_one_peak_sigma_ratio,
        )
        zero_one_counts_ext.append(zero_one_counts)
        zero_one_edges_ext.append(zero_one_edges)
        pedestals.append(pedestal)
        gains.append(gain)
        double_gauss_popts.append(double_gauss_popt)
        zero_one_ranges.append(zero_one_range)

    return zero_one_counts_ext, zero_one_edges_ext, pedestals, gains, double_gauss_popts, zero_one_ranges
        

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

def get_nonlinearity_at_ext(q, parabola_coeffs, parabola_pcovs, fit_range_right_ext):
    
    nonlinearity_at_q_ext = []
    for ext, parabola_coeff in enumerate(parabola_coeffs):
        print('\n***********************************************************')
        print(f'\nEXT {ext+1}:\n')
        parabola_pcov = parabola_pcovs[ext]
        fit_range_right = fit_range_right_ext[ext]
        nonlinearity_at_q = get_nonlinearity_at(q, parabola_coeff, parabola_pcov, fit_range_right)
        nonlinearity_at_q_ext.append(nonlinearity_at_q)
    return nonlinearity_at_q_ext

#---------------- Curves ----------------------------
def double_gauss(x, s0, m0, s1, m1, N0, N1):
    return N0 * np.exp(-(x-m0)**2/(2*s0**2)) + N1 * np.exp(-(x-m1)**2/(2*s1**2))

def parabola(x, a, b, c):
    return a*x**2 + b*x + c
