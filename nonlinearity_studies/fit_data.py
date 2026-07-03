"""
fit_data — split out of nonlinearity_studies.py.
"""
import math
import numpy as np
from pedestal_subtract.core import convert_to_electrons
from scipy.optimize import curve_fit
from scipy.signal import find_peaks as scipy_find_peaks

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
                   bin_factor=10,
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


def parabola(x, a, b, c):
    return a*x**2 + b*x + c
