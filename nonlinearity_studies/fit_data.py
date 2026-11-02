"""
fit_data — split out of nonlinearity_studies.py.
"""
import math
import numpy as np
from pedestal_subtract.core import convert_to_electrons
from scipy.optimize import curve_fit
from scipy.signal import find_peaks as scipy_find_peaks

# peak_finder_buffers are quoted normalized to this binning, i.e. the number in
# the configs is the buffer you would use at 10 bins per electron.
BUFFER_NORM_BIN_FACTOR = 10


def scale_buffers(buffers, bin_factor, n_ext):
    """Convert normalized peak_finder_buffers to bins at the binning in use.

    find_all_peaks subtracts the buffer from the bin factor to get find_peaks'
    minimum neighbor-peak distance (d = bin_factor - buffer, in bins), so a
    fixed buffer only means a fixed physical spacing at a fixed binning.
    Quoting it normalized to BUFFER_NORM_BIN_FACTOR bins per electron and
    scaling it here means retuning the bin factor leaves the peak separation
    the buffer asks for unchanged. ``bin_factor`` is bins per ELECTRON, which
    is what both pipelines key their coarse peak finding on -- the
    gain-independent one histograms ADU, so it converts its bins-per-electron
    factor to bins per ADU for the histogram itself, but the buffer is scaled
    here from the same bins-per-electron number.

    Accepts a scalar, a length-1 list, or one value per extension, and always
    returns one value per extension.
    """
    if not isinstance(buffers, list):
        buffers = [buffers] * n_ext
    elif len(buffers) == 1:
        buffers = buffers * n_ext
    return [b * bin_factor / BUFFER_NORM_BIN_FACTOR for b in buffers]


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
                   nonlinearity_peakfinder_bin_factor=10,
                   prominence=None):

    if flatten:
        data=np.array(data).flatten()

    if do_convert_to_electrons:
        data = convert_to_electrons(data, pedestal, gain, flatten=False)

    if range_left=='default':
        if do_convert_to_electrons:
            range_left = -0.5
        else:
            # If finding electron peaks in ADU, start at the left side of the zero
            # electron peak. 3 sigma out is far enough down the zero peak's flank to
            # capture it whole; going further left (the old -1 ADU of extra margin,
            # inherited from the electron-space default) only admits low-count bumps
            # in the pedestal tail, which the peak finder then numbers as peak 0.
            range_left=pedestal-3*noise

    hist_range = (range_left, range_right)

    if bins=='default':
       bins=math.floor((hist_range[1]-hist_range[0])*nonlinearity_peakfinder_bin_factor)

    counts, edges = np.histogram(data, bins=bins,range=hist_range)
    centers = 0.5 * (edges[1:] + edges[:-1])
    peaks, properties = scipy_find_peaks(counts, height=0, width=width*nonlinearity_peakfinder_bin_factor, distance=nonlinearity_peakfinder_bin_factor-buffer, prominence=prominence)

    # Also return the (electron-converted, if applicable) flattened data, so callers
    # that need to re-histogram a subset of it (e.g. plot_all_peaks zooming into a
    # narrow charge window at a finer nonlinearity_peakfinder_bin_factor than was used here) don't have to
    # redo the flattening/ADU->electron conversion themselves.
    return counts, edges, peaks, centers, properties, hist_range, data


def resolve_known_points(known_points, n_ext):
    """Per-extension list of known (charge, nonlinearity) points. A flat / [[q, nl], ...] value applies to every
    extension; a list with one entry per extension ([[[q, nl], ...], None, ...]) applies per extension."""
    per_ext = isinstance(known_points, (list, tuple)) and any(
        k is None or (isinstance(k, (list, tuple)) and len(k) and isinstance(k[0], (list, tuple))) for k in known_points)
    if not per_ext:
        return [known_points] * n_ext
    if len(known_points) != n_ext:
        raise ValueError(f'Per-extension known_points needs {n_ext} entries, got {len(known_points)}')
    return list(known_points)


def fit_nonlinearity(peaks, centers, pedestal, gain, fit_range_right, do_convert_to_electrons=False, fit_bounds_low=-100, fit_bounds_high=100, fix_intercept=True, known_points=None):
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
    x = peak_charge_e[:fit_idx]
    y = charge_minus_npeak[:fit_idx]

    # Known (charge, nonlinearity) pairs are fit alongside the peaks below fit_range_right.
    if known_points:
        known = np.asarray(known_points, dtype=float).reshape(-1, 2)
        x = np.concatenate([x, known[:, 0]])
        y = np.concatenate([y, known[:, 1]])

    if fix_intercept:
        # Nonlin(0) = 0 identically: peak number 0 is the zero-electron peak, so
        # its charge is the pedestal and charge_minus_npeak[0] = 0 by construction.
        # Left free, c absorbs the misfit of a parabola against a curve that stops
        # being quadratic well before fit_range_right -- with ~1500 peaks in the
        # fit, most of them at high charge, least squares buys residual there by
        # lifting the whole curve, which shows up as the parabola sitting ~0.5 e-
        # above a scatter that is pinned to the origin. Dropping c from the fit
        # removes that freedom; the returned coefficients keep their (a, b, c)
        # shape, with c = 0, so every caller of ``parabola`` is unaffected.
        coeff, pcov = curve_fit(_parabola_no_offset, x, y, maxfev=2000,
                                bounds=(fit_bounds_low, fit_bounds_high))
        parabola_coeff = np.array([coeff[0], coeff[1], 0.0])
        parabola_pcov = np.zeros((3, 3))
        parabola_pcov[:2, :2] = pcov
    else:
        parabola_coeff, parabola_pcov = curve_fit(parabola, x, y,
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


def _parabola_no_offset(x, a, b):
    """``parabola`` with the constant term held at zero; see fit_nonlinearity."""
    return a*x**2 + b*x
