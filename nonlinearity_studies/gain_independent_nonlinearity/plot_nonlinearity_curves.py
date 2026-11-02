"""
Computes nonlinearity from electron-converted charge spectra
using several different methods for computing the gain
1) Standard g from pedestal subtract zero/one peaks
2) Standard g(1) from sub-bin localized peaks
3) g_fit from entire fitted ADU spectrum
4) g_fit(i) for several i in [0,1500] where i is peak number including optimal
5) g_dir from ADU spectrum
"""

from pathlib import Path

import matplotlib.pyplot as plt

from pedestal_subtract.core import _finish_fig
from pedestal_subtract.plotting import format_title

from ..nonlinearity import get_nonlinearity_ext
from ..fit_data import parabola
from .compare_nonlinearity_methods import _position

"""Function that computes the residual nonlinearity by using the direct gain at each peak number"""

"""Compute the nonlinearity with different gains, then plot them all together on the same plot overlaid
The tricky thing is that different gains result in different conversions to e- units, but I think this is okay for a rough comparison.
"""
def plot_nonlinearity_curves(zero_one_gains, sub_bin_gains, fit_ADU_spectrum_coeffs, gain_optimization_out, test_i, g_dirs, peaks_ext, centers_ext, pedestals, fit_range_right_ext,
                             fix_intercept=True,
                             plot_individual=False, plot_together=False,
                             individual_figsize=(9, 6.5), subplots_figsize=(16, 9),
                             sharex=True, sharey=False, show_titles=True,
                             additional_title='', nimages=10, suptitle='Nonlinearity by Gain Method',
                             save_plots=False, show_plots=True,
                             fig_path='./', file='gain_indep_nonlinearity', dpi=350):
    """Given peaks, centers, pedestals in ADU and gains (one per extension) from zero/one double gaussian fit, zero/one sub-bin method, g_fit at several i including the best i, and direct/avg gain, compute nonlinearity curves for all and plot on same axes.
    Generates ``g_fits_at_i`` as an intermediate step given ``test_i`` (list of indices to sample from) and ``fit_ADU_spectrum_coeffs``, the output of ``fit_charge_spectrum_running_ext``. It is a list containing one list per index ``i`` of g_fit(i) for each extension.
    ``g_fits_at_best_i`` is generated from ``gain_optimization_out``. It is a list of four lists, the first element of each being the index ``i`` of each extension corresponding to the optimal nonlinearity, and the second being the ``g_fit(i)`` at the corresponding ``i``.

    ``fix_intercept`` is passed through to every parabola fit drawn here; see
    ``fit_data.fit_nonlinearity``. With it off, a fit_range_right that reaches
    past the quadratic regime leaves the curves offset from the scatter at low
    peak numbers, since the free constant term absorbs the high-charge misfit.

    The figure-size / shared-axis / title / saving arguments are the same set the
    plots in ``plots.py`` take, so ``run`` can hand this the one ``plot_kwargs``
    dict it builds for all of them.
    """

    fig_name = Path(fig_path) / f'{Path(file).stem}_nonlinearity_curves'

    if plot_individual or plot_together:
        zero_one_peak_charge_e_ext, zero_one_charge_minus_npeak_ext, zero_one_parabola_coeffs, _ = get_nonlinearity_ext(peaks_ext, centers_ext, pedestals, zero_one_gains, fit_range_right_ext, do_convert_to_electrons=True, fit_bounds_low=-100, fit_bounds_high=100, fix_intercept=fix_intercept)

        sub_bin_peak_charge_e_ext, sub_bin_charge_minus_npeak_ext, sub_bin_parabola_coeffs, _ = get_nonlinearity_ext(peaks_ext, centers_ext, pedestals, sub_bin_gains, fit_range_right_ext, do_convert_to_electrons=True, fit_bounds_low=-100, fit_bounds_high=100, fix_intercept=fix_intercept)

        # An i outside an extension's fit domain has no g_fit(i), and _position
        # returns None there. Indexing an ndarray with None inserts an axis
        # rather than raising, which would put a whole array where a scalar gain
        # belongs, so only the indices every extension fit are kept. The labels
        # follow test_i_used, not test_i, so they stay aligned with the curves.
        test_i_used = []
        g_fits_at_i = []
        for i in test_i:
            positions = [_position(fit_ADU_spectrum_coeffs[ext], i) for ext in range(4)]
            if any(p is None for p in positions):
                continue
            test_i_used.append(i)
            g_fits_at_i.append([fit_ADU_spectrum_coeffs[ext]["g_fit"][p]
                                for ext, p in enumerate(positions)])

        g_fits_i_peak_charge_e_ext = []
        g_fits_i_charge_minus_npeak_ext = []
        g_fits_i_parabola_coeffs = []
        for g_fits in g_fits_at_i:
            g_fit_peak_charge_e_ext, g_fit_charge_minus_npeak_ext, g_fit_parabola_coeffs, _ = get_nonlinearity_ext(peaks_ext, centers_ext, pedestals, g_fits, fit_range_right_ext, do_convert_to_electrons=True, fit_bounds_low=-100, fit_bounds_high=100, fix_intercept=fix_intercept)
            g_fits_i_peak_charge_e_ext.append(g_fit_peak_charge_e_ext)
            g_fits_i_charge_minus_npeak_ext.append(g_fit_charge_minus_npeak_ext)
            g_fits_i_parabola_coeffs.append(g_fit_parabola_coeffs)

        # best_gain is the gain the nonlinearity is recomputed with; best_i is
        # the peak number it came from, and is only used to label the curve.
        g_fit_at_best_i = [gain_optimization_out[ext]["best_gain"] for ext in range(4)]
        best_i_at_best = [gain_optimization_out[ext]["best_i"] for ext in range(4)]
        g_fit_best_i_peak_charge_e_ext, g_fit_best_i_charge_minus_npeak_ext, g_fit_best_i_parabola_coeffs, _ = get_nonlinearity_ext(peaks_ext, centers_ext, pedestals, g_fit_at_best_i, fit_range_right_ext, do_convert_to_electrons=True, fit_bounds_low=-100, fit_bounds_high=100, fix_intercept=fix_intercept)

        g_dir_peak_charge_e_ext, g_dir_charge_minus_npeak_ext, g_dir_parabola_coeffs, _ = get_nonlinearity_ext(peaks_ext, centers_ext, pedestals, g_dirs, fit_range_right_ext, do_convert_to_electrons=True, fit_bounds_low=-100, fit_bounds_high=100, fix_intercept=fix_intercept)

    if plot_individual:
        for ext, peaks in enumerate(peaks_ext):
            fig, ax = plt.subplots(1, 1, figsize=individual_figsize, constrained_layout=True)

            ax.grid()

            ax.scatter(zero_one_peak_charge_e_ext[ext], zero_one_charge_minus_npeak_ext[ext], color='b', alpha=0.7, s=4)
            ax.plot(zero_one_peak_charge_e_ext[ext], parabola(zero_one_peak_charge_e_ext[ext], *zero_one_parabola_coeffs[ext]),
                                    label=f'Zero-One Gain: g = {zero_one_gains[ext]:.4g} ADU/e-')

            ax.scatter(sub_bin_peak_charge_e_ext[ext], sub_bin_charge_minus_npeak_ext[ext], color='b', alpha=0.7, s=4)
            ax.plot(sub_bin_peak_charge_e_ext[ext], parabola(sub_bin_peak_charge_e_ext[ext], *sub_bin_parabola_coeffs[ext]),
                                                label=f'Sub-Bin Gain: g(1) = {sub_bin_gains[ext]:.4g} ADU/e-')

            for n, index in enumerate(test_i_used):
                ax.scatter(g_fits_i_peak_charge_e_ext[n][ext], g_fits_i_charge_minus_npeak_ext[n][ext], color='b', alpha=0.7, s=4)
                ax.plot(g_fits_i_peak_charge_e_ext[n][ext], parabola(g_fits_i_peak_charge_e_ext[n][ext], *g_fits_i_parabola_coeffs[n][ext]),  label=f'g_fit(i={index}) = {g_fits_at_i[n][ext]:.4g} ADU/e-', )

            ax.scatter(g_fit_best_i_peak_charge_e_ext[ext], g_fit_best_i_charge_minus_npeak_ext[ext], color='b', alpha=0.7, s=4)
            ax.plot(g_fit_best_i_peak_charge_e_ext[ext], parabola(g_fit_best_i_peak_charge_e_ext[ext], *g_fit_best_i_parabola_coeffs[ext]), label=f'Optimal g_fit(i={best_i_at_best[ext]}) = {g_fit_at_best_i[ext]:.4g} ADU/e-')

            ax.scatter(g_dir_peak_charge_e_ext[ext], g_dir_charge_minus_npeak_ext[ext], color='b', alpha=0.7, s=4)
            ax.plot(g_dir_peak_charge_e_ext[ext], parabola(g_dir_peak_charge_e_ext[ext], *g_dir_parabola_coeffs[ext]),
            label=f'Direct Gain: g_dir(500e-) = {g_dirs[ext]:.4g} ADU/e-')

            ax.legend(loc="upper right")
            if show_titles:
                ax.set_title(format_title(additional_title, suptitle, nimages, f'EXT {ext + 1}'),
                             pad=10)
            ax.set_xlabel(r'Measured Pixel Charge ($e^-$)')
            ax.set_ylabel(r'Measured Pixel Charge - Peak n. ($e^-$) ')

            if save_plots:
                output_path = fig_name.with_stem(fig_name.stem + f'_EXT{ext+1}').with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            _finish_fig(show_plots)

    if plot_together:
        fig, axs = plt.subplots(2, 2, figsize=subplots_figsize, constrained_layout=True,
                                sharex=sharex, sharey=sharey)
        axs=axs.flatten()
        if show_titles:
            fig.suptitle(format_title(additional_title, suptitle, nimages))

        for ext in range(len(zero_one_peak_charge_e_ext)):
            ax = axs[ext]
            ax.grid()

            ax.scatter(zero_one_peak_charge_e_ext[ext], zero_one_charge_minus_npeak_ext[ext], color='b', alpha=0.7, s=4)
            ax.plot(zero_one_peak_charge_e_ext[ext], parabola(zero_one_peak_charge_e_ext[ext], *zero_one_parabola_coeffs[ext]),
                                    label=f'Zero-One Gain: g = {zero_one_gains[ext]:.4g} ADU/e-')

            ax.scatter(sub_bin_peak_charge_e_ext[ext], sub_bin_charge_minus_npeak_ext[ext], color='b', alpha=0.7, s=4)
            ax.plot(sub_bin_peak_charge_e_ext[ext], parabola(sub_bin_peak_charge_e_ext[ext], *sub_bin_parabola_coeffs[ext]),
                                                label=f'Sub-Bin Gain: g = {sub_bin_gains[ext]:.4g} ADU/e-')

            for n, index in enumerate(test_i_used):
                ax.scatter(g_fits_i_peak_charge_e_ext[n][ext], g_fits_i_charge_minus_npeak_ext[n][ext], color='b', alpha=0.7, s=4)
                ax.plot(g_fits_i_peak_charge_e_ext[n][ext], parabola(g_fits_i_peak_charge_e_ext[n][ext], *g_fits_i_parabola_coeffs[n][ext]), label=f'g_fit(i={index}) = {g_fits_at_i[n][ext]:.4g} ADU/e-', )

            ax.scatter(g_fit_best_i_peak_charge_e_ext[ext], g_fit_best_i_charge_minus_npeak_ext[ext], color='b', alpha=0.7, s=4)
            ax.plot(g_fit_best_i_peak_charge_e_ext[ext], parabola(g_fit_best_i_peak_charge_e_ext[ext], *g_fit_best_i_parabola_coeffs[ext]), label=f'Optimal g_fit(i={best_i_at_best[ext]}) = {g_fit_at_best_i[ext]:.4g} ADU/e-')

            ax.scatter(g_dir_peak_charge_e_ext[ext], g_dir_charge_minus_npeak_ext[ext], color='b', alpha=0.7, s=4)
            ax.plot(g_dir_peak_charge_e_ext[ext], parabola(g_dir_peak_charge_e_ext[ext], *g_dir_parabola_coeffs[ext]),
            label=f'Direct Gain: g_dir(500e-) = {g_dirs[ext]:.4g} ADU/e-')

            ax.legend(loc="upper right")
            ax.set_title(f'EXT {ext + 1}')
            ax.set_xlabel(r'Measured Pixel Charge ($e^-$)')
            ax.set_ylabel(r'Measured Pixel Charge - Peak n. ($e^-$) ')

        if save_plots:
            output_path = fig_name.with_suffix('.jpeg')
            plt.savefig(str(output_path), dpi=dpi)
            print(f'Saved plot to {output_path}')
        _finish_fig(show_plots)



    

    
