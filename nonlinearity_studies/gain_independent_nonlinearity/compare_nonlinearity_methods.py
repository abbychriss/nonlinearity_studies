"""
Computes variance in nonlinearity with respect to the gain
Sweeps gain from [1.0, 1.4] in steps of 0.01 and samples nonlin(500e-)
Computes standard deviation across the gains

Computes variance in alpha and alpha' with respect to peak number
Samples alpha and alpha' at several peak numbers and computes standard deviation
Want to understand how much the spread in alpha/alpha' (fit parameters) is due to the spread in gain, also compute the maximum change in g_fit corresponding to the different alphas

Plots alpha and/or alpha' in same plot against index with g_fit(i) (from the alpha calculation) as secondary axis
"""

import numpy as np
import matplotlib.pyplot as plt

from ..nonlinearity import get_nonlinearity_ext, get_nonlinearity_at_ext
from .spectrum_fit import fit_charge_spectrum_running_ext
from .gain_scan import scan_gain_nonlinearity_ext
from .deltas import fit_alpha_from_S_ext

def nonlin_at_q_vs_gain(peaks_ext, centers_ext, pedestals, q=500, gain_range=np.linspace(1.0, 1.4, 400), fit_range_right_ext=1000.0, fit_bounds_low=-100, fit_bounds_high=100, fix_intercept=True):
    """Returns a 2d list of shape (len(gain_range), 4). The first four elements of each list are the nonlinearity at q (in electrons) in each extension. Also returns ``gain_range`` and ``q``, the charge in e- at which the nonlinearity is computed."""

    nonlinearity_at_gains_ext = []

    # get_nonlinearity_ext promotes a scalar fit range to one per extension
    # itself, but get_nonlinearity_at_ext indexes it as a list, so it is
    # promoted here for both.
    if not isinstance(fit_range_right_ext, list):
        fit_range_right_ext = [fit_range_right_ext] * len(peaks_ext)

    for gain in gain_range:

        gains = [gain for i in range(4)]

        peak_charge_e_ext, charge_minus_npeak_ext, parabola_coeffs, parabola_pcovs = get_nonlinearity_ext(peaks_ext, centers_ext, pedestals, gains, fit_range_right_ext, do_convert_to_electrons=True, fit_bounds_low=fit_bounds_low, fit_bounds_high=fit_bounds_high, fix_intercept=fix_intercept)

        nonlin_at_q_ext = get_nonlinearity_at_ext(q, parabola_coeffs, parabola_pcovs, fit_range_right_ext, save_path=None)

        nonlinearity_at_gains_ext.append(nonlin_at_q_ext)

    # Transpose list so that it's the standard shape (one list per extension rather than one list per gain tested)
    nonlinearity_at_gains_ext = np.array(nonlinearity_at_gains_ext).T

    return nonlinearity_at_gains_ext, gain_range, q

def _position(run, i):
        """Array position of peak number ``i`` in a running-fit dict; None if not fit."""
        if i is None:
            return None
        pos = int(np.searchsorted(run["i"], i))
        return pos if pos < run["i"].size and run["i"][pos] == i else None

def alpha_vs_index_ext(fit_ADU_spectrum_coeffs, gain_optimization_out, test_i=[]):
    """returns list of dicts (one per extension) of ``alpha`` and corresponding index at the maximum peak number ``alpha_at_max_i`, ``max_i``, optimal (lowest |nonlin|) index ``alpha_at_best_i``, ``best_i``, and other indices supplied by user ``alpha_at_test_i`, ``test_i`, the standard deviation of the alphas sampled ``std_alphas``, the maximum difference of the gains sampled ``max_gain_diff``. Also returns a list of dicts (one per extension) ``index_alpha_gfit_sorted_ext`` containing all the indices, alphas, and gains sampled, sorted by index."""

    if test_i=='all':
        min_max_i = min([fit_ADU_spectrum_coeffs[ext]["i"][-1] for ext in range(len(fit_ADU_spectrum_coeffs))])
        max_min_i = max([fit_ADU_spectrum_coeffs[ext]["i"][0] for ext in range(len(fit_ADU_spectrum_coeffs))])
        
        test_i = list(range(max_min_i,min_max_i+1))

    alpha_dict_ext = [{},{},{},{}]
    index_alpha_gfit_sorted_ext = [{},{},{},{}]

    for ext in range(4):
        
        max_i = fit_ADU_spectrum_coeffs[ext]["i"][-1]
        alpha_at_max_i = fit_ADU_spectrum_coeffs[ext]["alpha"][-1]
        alpha_dict_ext[ext]["max_i"] = max_i
        alpha_dict_ext[ext]["alpha_at_max_i"] = alpha_at_max_i

        # best_i is a peak number, not a position in the running arrays, so it
        # goes through _position before it can index alpha or g_fit. It is None
        # when every candidate gain gave a NaN nonlinearity; it is recorded as
        # NaN rather than None so it survives the vstack/argsort below.
        best_i = gain_optimization_out[ext]["best_i"]
        best_i_position = _position(fit_ADU_spectrum_coeffs[ext], best_i)
        alpha_at_best_i = fit_ADU_spectrum_coeffs[ext]["alpha"][best_i_position] if best_i_position is not None else np.nan
        g_fit_at_best_i = fit_ADU_spectrum_coeffs[ext]["g_fit"][best_i_position] if best_i_position is not None else np.nan
        alpha_dict_ext[ext]["best_i"] = best_i if best_i is not None else np.nan
        alpha_dict_ext[ext]["alpha_at_best_i"] = alpha_at_best_i

        g_fits_sampled = [fit_ADU_spectrum_coeffs[ext]["g_fit"][-1], g_fit_at_best_i]

        # Only the test indices inside this extension's fit domain can be
        # sampled, so the rest are dropped and what is recorded is what was
        # actually sampled rather than what was asked for -- otherwise the
        # index / alpha / g_fit rows come out different lengths below.
        test_pos = [(i, _position(fit_ADU_spectrum_coeffs[ext], i)) for i in test_i]
        test_pos = [(i, p) for i, p in test_pos if p is not None]
        test_i_used = [i for i, _ in test_pos]
        alpha_dict_ext[ext]["test_i"] = test_i_used

        alpha_at_test_i = [fit_ADU_spectrum_coeffs[ext]["alpha"][p] for _, p in test_pos]
        alpha_dict_ext[ext]["alpha_at_test_i"] = alpha_at_test_i

        alphas_sampled = [alpha_at_max_i, alpha_at_best_i, *alpha_at_test_i]
        indices_sampled = [max_i, alpha_dict_ext[ext]["best_i"], *test_i_used]
        g_fits_sampled.extend([fit_ADU_spectrum_coeffs[ext]["g_fit"][p] for _, p in test_pos])

        all_values_sampled = np.vstack((indices_sampled, alphas_sampled, g_fits_sampled))
        indices = np.argsort(all_values_sampled[0, :])
        all_values_sampled_sorted = all_values_sampled[:, indices]

        index_alpha_gfit_sorted_ext[ext]["i"] = all_values_sampled_sorted[0]
        index_alpha_gfit_sorted_ext[ext]["alpha"] = all_values_sampled_sorted[1]
        index_alpha_gfit_sorted_ext[ext]["g_fit"] = all_values_sampled_sorted[2]

        # A missing best_i leaves a NaN in both lists; the builtin max/min and
        # np.std would propagate it into the summary instead of reporting the
        # spread of the samples that do exist.
        max_gain_diff = np.nanmax(g_fits_sampled) - np.nanmin(g_fits_sampled)
        std_alphas = np.nanstd(alphas_sampled)

        alpha_dict_ext[ext]["std_alphas"] = std_alphas
        alpha_dict_ext[ext]["max_gain_diff"] = max_gain_diff

    return alpha_dict_ext, index_alpha_gfit_sorted_ext

    
def alpha_prime_vs_index_ext(S_ext, alpha_dict_ext, test_i=[]):
    """ Returns a list of dicts (one per extension) containing the indices sampled ``i`` and the corresponding ``alpha'``
    """

    all_indices_sampled_ext =  [[],[],[],[]]
    index_alpha_prime_ext = [{},{},{},{}]

    for ext, alpha_dict in enumerate(alpha_dict_ext):
        max_i = alpha_dict["max_i"]
        best_i = alpha_dict["best_i"]
        # A missing best_i is carried as NaN, and n_max=NaN would raise inside
        # peak_number_mask's int() while n_max=None would silently fit every
        # peak instead of skipping, so non-finite indices are dropped here.
        all_indices_sampled_ext[ext].extend(
            [i for i in [max_i, best_i, *test_i] if i is not None and np.isfinite(i)])

    for ext in range(4):
        index_alpha_prime_ext[ext]["i"] = all_indices_sampled_ext[ext]
        index_alpha_prime_ext[ext]["alpha_prime"] = [fit_alpha_from_S_ext(S_ext, i)[ext]["slope"] for i in all_indices_sampled_ext[ext]]

    return index_alpha_prime_ext

