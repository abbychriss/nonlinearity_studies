"""
nonlinearity — split out of nonlinearity_studies.py.
"""
import numpy as np
from pathlib import Path

from .fit_data import find_all_peaks, fit_nonlinearity

def get_nonlinearity_at(q, parabola_coeff, parabola_pcov=None, fit_range_right=None, print_values=False, out_lines=None):
    a = parabola_coeff[0]
    b = parabola_coeff[1]
    c = parabola_coeff[2]

    def _emit(s):
        print(s)
        if out_lines is not None:
            out_lines.append(s)

    _emit(f'Parabola best fit: Nonlin(q) = {np.round(a,5)}*q^2 + {np.round(b,5)}*q + {np.round(c,5)}\n')

    if isinstance(q, (float, int, np.floating, np.integer)):
        nonlinearity_at_q = np.round(a*q**2 + b*q + c, 3)
    elif isinstance(q, (list, np.ndarray)):
        nonlinearity_at_q = [np.round((a*q_i**2 + b*q_i + c), 3) for q_i in q]
    else:
        raise TypeError(f'q must be a number or list, got {type(q)}')

    if isinstance(q, (float, int, np.floating, np.integer)):
        _emit(f'Interpolated nonlinearity at {q} e- = {nonlinearity_at_q} e-\n')
    else:
        for i in range(len(q)):
            _emit(f'Interpolated nonlinearity at {q[i]} e- = {nonlinearity_at_q[i]} e-\n')

    if print_values:
        if fit_range_right is not None:
            _emit(f'Parabola fit was performed up to {fit_range_right} e-\n')
        if parabola_pcov is not None:
            _emit(f'Covariance of fit = {parabola_pcov}\n')

    return nonlinearity_at_q


def get_all_peaks_ext(data_ext, widths, buffers, pedestals, double_gauss_popts, gains, bins='default', flatten=True, do_convert_to_electrons=True, range_left='default', range_right=2000, nonlinearity_peakfinder_bin_factor=10, prominences=None, print_values=False):
    counts_ext = []
    edges_ext = []
    peaks_ext = []
    centers_ext = []
    hist_ranges = []
    data_electrons_ext = []

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

        counts, edges, peaks, centers, properties, hist_range, data_electrons = find_all_peaks(data,
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
                                                                               nonlinearity_peakfinder_bin_factor=nonlinearity_peakfinder_bin_factor,
                                                                               prominence=prominence)

        counts_ext.append(counts)
        edges_ext.append(edges)
        peaks_ext.append(peaks)
        centers_ext.append(centers)
        hist_ranges.append(hist_range)
        data_electrons_ext.append(data_electrons)
    if print_values:
        print('***********************************************************')

    return counts_ext, edges_ext, peaks_ext, centers_ext, hist_ranges, data_electrons_ext


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


def get_nonlinearity_at_ext(q, parabola_coeffs, parabola_pcovs, fit_range_right_ext, save_path=None):
    out_lines = [] if save_path is not None else None

    nonlinearity_at_q_ext = []
    for ext, parabola_coeff in enumerate(parabola_coeffs):
        print('\n***********************************************************')
        header = f'EXT {ext+1}:'
        print(f'\n{header}\n')
        if out_lines is not None:
            out_lines.append('')
            out_lines.append(header)
            out_lines.append('')
        parabola_pcov = parabola_pcovs[ext]
        fit_range_right = fit_range_right_ext[ext]
        nonlinearity_at_q = get_nonlinearity_at(q, parabola_coeff, parabola_pcov, fit_range_right, out_lines=out_lines)
        nonlinearity_at_q_ext.append(nonlinearity_at_q)

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        with save_path.open('w', encoding='utf-8') as f:
            f.write('\n'.join(out_lines))
        print(f'Saved nonlinearity-at-charge text output to {save_path}')

    return nonlinearity_at_q_ext
