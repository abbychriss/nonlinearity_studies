"""
plotting — split out of nonlinearity_studies.py.
"""
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from pathlib import Path
from pedestal_subtract.core import _finish_fig
from pedestal_subtract.plotting import format_title

from .fit_data import _gauss_single, _make_comb_model, parabola, resolve_known_points
from .resolution import classify_resolution
from .summary import _normalize_charges

_RESOLUTION_VERDICT_COLORS = {
    'well resolved': 'tab:green', 'marginal': 'tab:orange', 'unresolved': 'tab:red',
}


def set_plot_style(font_size, legend_font_size, use_tex=False, font_family=None):
    """Apply the text style of the font_size / legend_font_size / use_tex / font_family config options."""
    # Uniform text size for everything except legends, which use legend_font_size.
    plt.rcParams.update({'font.size': font_size, 'legend.fontsize': legend_font_size, 'legend.title_fontsize': legend_font_size,
                         'axes.titlesize': font_size, 'axes.labelsize': font_size,
                         'xtick.labelsize': font_size, 'ytick.labelsize': font_size,
                         'text.usetex': use_tex})
    # Same font handling as meta_studies: use_tex always means Computer Modern Roman and
    # overrides font_family.
    if use_tex:
        plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Computer Modern Roman']})
    elif font_family:
        plt.rcParams.update({'font.family': font_family, 'axes.unicode_minus': False})


def _sci(value, sig=3):
    r"""A number as mathtext scientific notation: -9.07e-06 -> '-9.073 \cdot 10^{-6}'.

    '%e' would put a literal 'e-06' in the label, which reads as a typo next to
    the LaTeX around it. The x10 factor is dropped when the exponent is zero, so
    ordinary-sized numbers stay plain. Must be used inside a mathtext ($...$)
    span.
    """
    value = float(value)
    if not np.isfinite(value):
        return r'\mathrm{nan}'
    mantissa, exponent = f'{value:.{sig}e}'.split('e')
    exponent = int(exponent)
    if exponent == 0:
        return mantissa
    return r'%s \cdot 10^{%d}' % (mantissa, exponent)


def _term(value, sig=3, sci=True):
    r"""A term to append to an equation, its sign carried as the operator.

    Returns e.g. r'- 9.073 \cdot 10^{-6}', never '+ (-9.073e-06)': a negative
    coefficient reads as a subtraction from what came before instead of an added
    negative. ``sci=False`` keeps the magnitude in plain decimal, for
    coefficients that are not small enough to want an exponent.
    """
    value = float(value)
    body = _sci(abs(value), sig) if sci else f'{abs(value):.{sig}f}'
    return r'%s %s' % ('-' if value < 0 else '+', body)


def _parabola_label(coeff):
    r"""Legend label for a parabola fit, with negative coefficients subtracted rather than
    rendered as '+ -a' (e.g. '-1.234 \cdot 10^{-5} x^2 - 0.045 x + 1.230')."""
    a, b, c = coeff
    return r'$%s x^2 %s x %s$' % (_sci(a), _term(b, sci=False), _term(c, sci=False))


def _xlim_window(xlim_e):
    """Resolve xlim_e to a (left, right) tuple, or None if it's not a concrete window (e.g. 'none')."""
    is_window = (
        isinstance(xlim_e, (list, tuple)) and len(xlim_e) == 2
        and all(isinstance(v, (int, float, np.integer, np.floating)) for v in xlim_e)
    )
    if not is_window:
        return None
    return min(xlim_e), max(xlim_e)


def _auto_ylim_in_window(centers, counts, xlim_e, yscale='linear', pad=0.05):
    """Tight y-limits framing only the histogram bars inside the x-window.

    Used when ylim is 'none' for plot_all_peaks: rather than letting matplotlib
    autoscale to the full distribution (dominated by the tall low-charge peaks), zoom
    to the tallest bar whose center falls within ``xlim_e`` so the peaks in the
    displayed x-range are actually visible. ``xlim_e`` may be a (left, right) pair; any
    other value (e.g. the sentinel 'none') means "use every bar". Returns (bottom, top),
    or None when the window contains no bars (so the caller leaves autoscale untouched).
    """
    centers = np.asarray(centers, dtype=float)
    counts = np.asarray(counts, dtype=float)

    window = _xlim_window(xlim_e)
    if window is not None:
        left, right = window
        in_window = (centers >= left) & (centers <= right)
    else:
        in_window = np.ones(centers.shape, dtype=bool)

    selected = counts[in_window]
    selected = selected[np.isfinite(selected)]
    if selected.size == 0:
        return None

    top = float(np.max(selected))
    if top <= 0:
        return None

    if yscale == 'log':
        positive = selected[selected > 0]
        bottom = float(np.min(positive)) * 0.5 if positive.size else 0.5
        top *= (1.0 + pad) * 1.5
    else:
        bottom = 0.0
        top *= (1.0 + pad)
    return (bottom, top)


def _windowed_hist(data, window, bin_factor):
    """Re-histogram ``data`` over just ``window`` at ``bin_factor`` bins per x-axis unit.

    The unit is whatever the data is in -- electrons for the standard analysis, ADU for
    the gain-independent one. Used by plot_all_peaks when ``data_ext``/``window_bin_factor``
    are given, so a zoomed-in view can use a bin width independent of (typically finer than)
    whatever peak-finder bin factor was used over the full spectrum in get_all_peaks_ext.
    Only the samples inside ``window`` are binned, which keeps this fast even though
    peak-finding itself still ran over the full charge range.
    """
    bins = math.floor((window[1] - window[0]) * bin_factor)
    data = np.asarray(data)
    data_window = data[(data > window[0]) & (data < window[1])]
    counts, edges = np.histogram(data_window, bins=bins, range=window)
    centers = 0.5 * (edges[1:] + edges[:-1])
    return counts, centers


def _draw_hist(ax, centers, counts, bin_width):
    """Filled step histogram as ONE artist (ax.bar makes a Rectangle per bin: ~50 s for 4 x 16k bins)."""
    ax.stairs(counts, np.append(centers - bin_width / 2, centers[-1] + bin_width / 2), fill=True)


BINS_PER_E = r'bins/$e^-$'


def _bin_factor_text(bin_factor, unit=BINS_PER_E, gain=None):
    """Render one bin factor with its unit, e.g. ``40 bins/ADU = 102 bins/e-``.

    ``unit`` is whatever the x axis is binned in: bins per electron for the
    standard analysis, bins per ADU for the gain-independent one. When the axis
    is in ADU the same number means a much narrower bin, so ``gain`` (that
    extension's ADU per electron, used for display only) adds the electron
    equivalent alongside it, which is what makes an ADU plot's binning directly
    comparable to a bins-per-electron plot's. Returns None when there is no
    factor to report.
    """
    if bin_factor is None:
        return None
    text = f'{bin_factor:.4g} {unit}'
    if unit != BINS_PER_E and gain is not None and np.isfinite(gain):
        text += f' = {bin_factor * gain:.0f} {BINS_PER_E}'
    return text


SUB_BIN_FIT_COLOR = 'purple'
GAUSS_FIT_COLOR = 'green'


def _draw_sub_bin_fits(ax, fits, bin_width, alpha=0.6, n_pts=25):
    """Overlay every successful sub-bin fit on the drawn histogram.

    The log-parabola's vertex form IS a Gaussian, so both methods' curves come
    from the same (amp, mu, sigma) triple and differ only in which peaks they
    cover and how wide a window they were fit over. Each curve is drawn over
    exactly the window that produced it (fit_lo..fit_hi), and all curves of one
    method are concatenated into a single NaN-separated polyline -- one Line2D
    and one legend entry for ~1500 peaks instead of one per peak.

    ``amp`` is a height in counts of the histogram the fit was made on, which is
    binned at sub_bin_peakfinder_bin_factor; ``bin_width`` is the bin width of
    the histogram actually being drawn. Counts scale with bin width, so the
    curves are rescaled by their ratio to sit on the drawn bars.
    """
    if fits is None or fits['ok'].size == 0:
        return
    amp_scale = abs(bin_width) / abs(fits['bin_width'])
    n_total = fits['ok'].size
    t = np.linspace(0.0, 1.0, n_pts)

    for mask, color, name in (
            (fits['ok'] & ~fits['refit'], SUB_BIN_FIT_COLOR, 'Log parabola fit'),
            (fits['ok'] & fits['refit'], GAUSS_FIT_COLOR, 'Gaussian fits')):
        idx = np.flatnonzero(mask)
        if idx.size == 0:
            continue
        lo = fits['fit_lo'][idx][:, None]
        hi = fits['fit_hi'][idx][:, None]
        x = lo + (hi - lo) * t[None, :]
        y = _gauss_single(x, (fits['amp'][idx] * amp_scale)[:, None],
                          fits['mu'][idx][:, None], fits['sigma'][idx][:, None])
        # a NaN between consecutive curves lifts the pen, so the single line
        # never connects one peak's window to the next
        nan_col = np.full((idx.size, 1), np.nan)
        ax.plot(np.concatenate([x, nan_col], axis=1).ravel(),
                np.concatenate([y, nan_col], axis=1).ravel(),
                color=color, alpha=alpha, linewidth=1,
                label=f'{name} ({idx.size} / {n_total} peaks found)')


def _peak_line_legend(ax, raw, sub_bin, linecolor, linestyle,
                      sub_bin_linecolor, sub_bin_linestyle,
                      raw_bin_factor=None, sub_bin_bin_factor=None,
                      display_bin_factor=None, unit=BINS_PER_E, gain=None,
                      fontsize=None):
    """Name the two peak-location methods the vertical lines come from, plus the
    binning the histogram itself is drawn at.

    The lines themselves are drawn once per peak, so the legend is built from
    empty proxy lines instead -- one entry per method, not one per peak. Each
    method carries the bin factor its peaks were found at, which the two need not
    share, so it is appended per entry rather than quoted once for the plot. The
    bars can be drawn at a third, independent factor (plot_all_peaks' window
    rebinning), so ``display_bin_factor`` gets its own entry whenever the drawn
    binning is not simply the peak finder's.
    """
    def _label(name, bin_factor, kind=None):
        text = _bin_factor_text(bin_factor, unit, gain)
        if text is None:
            return name
        return f'{name} ({text})' if kind is None else f'{name} ({kind}: {text})'

    # The two peak entries name no kind: the method is already in the entry's own
    # name, so only the histogram's factor needs "display:" to say it is the
    # binning the bars are drawn at rather than one peaks were found at.
    # TEMPORARILY disabled: scipy peaks legend entry (the lines themselves are still drawn).
    # if raw:
    #     ax.plot([], [], linestyle=linestyle, color=linecolor,
    #             label=_label('scipy peaks', raw_bin_factor))
    if sub_bin:
        ax.plot([], [], linestyle=sub_bin_linestyle, color=sub_bin_linecolor,
                label=_label('sub-bin peaks', sub_bin_bin_factor))
    if display_bin_factor is not None:
        ax.plot([], [], marker='s', linestyle='none', color='C0',
                label=_label('histogram', display_bin_factor, 'display'))
    # Built from whatever is labeled on the axes, so the sub-bin fit curves
    # drawn by _draw_sub_bin_fits get their entries without being passed here.
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=fontsize)


def plot_all_peaks(counts_ext,
                   peaks_ext,
                   centers_ext,
                   xlim=(1000,1010),
                   ylim='none',
                   yscale='linear',
                   data_ext=None,
                   window_bin_factor=None,
                   plot_individual=False, plot_together=False,
                   raw_peak_location_vertical_line=True, linecolor='r', linestyle='--',
                   raw_peak_number=True, peak_number_labels_size=8,
                   sub_bin_peaks_ext=None, sub_bin_linecolor='tab:orange', sub_bin_linestyle='--',
                   sub_bin_peak_location_vertical_line=True, sub_bin_peak_number=True,
                   sub_bin_fits_ext=None, sub_bin_fit_alpha=0.6, legend_fontsize=None,
                   peakfinder_bin_factor=None, sub_bin_bin_factor=None,
                   bin_factor_unit=BINS_PER_E, bin_factor_gains=None,
                   individual_figsize=(7,6), subplots_figsize=(13,9),
                   additional_title='',
                   suptitle='Pixel Charge Peak-Finder',
                   nimages=10,
                   sharex=True,
                   sharey=True,
                   show_titles=True,
                   save_plots=False,
                   show_plots=True,
                   fig_path='./', file='peak_finder',
                   dpi=350,
                   xlabel=r'Charge ($e^-$)'):

    fig_path = Path(fig_path)
    if file != 'peak_finder':
        base_name = file[:-5] + '_peak_finder'
    else:
        base_name = file
    fig_name = fig_path / base_name

    # The bars carry their own binning only when the window rebinning below runs;
    # otherwise they are the peak finder's histogram, already named in the legend.
    display_bin_factor = window_bin_factor if data_ext is not None else None

    if plot_individual:
        for ext, peaks in enumerate(peaks_ext):
            centers=centers_ext[ext]
            xlim_e = xlim[ext] if isinstance(xlim, list) else xlim
            ylim_e = ylim[ext] if isinstance(ylim, list) else ylim

            if data_ext is not None and window_bin_factor is not None:
                window = _xlim_window(xlim_e)
                if window is None:
                    raise ValueError(
                        f'plot_all_peaks: data_ext/window_bin_factor rebinning requires a '
                        f'concrete (left, right) xlim window, got xlim={xlim_e!r} for EXT {ext + 1}.'
                    )
                counts, plot_centers = _windowed_hist(data_ext[ext], window, window_bin_factor)
            else:
                counts = counts_ext[ext]
                plot_centers = centers
            bin_width = plot_centers[1] - plot_centers[0]

            fig, ax = plt.subplots(1, 1, figsize=individual_figsize, constrained_layout=True)
            if show_titles:
                ax.set_title(format_title(additional_title, suptitle, nimages, f'EXT {ext + 1}'), pad=10)
            _draw_hist(ax, plot_centers, counts, bin_width)
            if sub_bin_fits_ext is not None:
                _draw_sub_bin_fits(ax, sub_bin_fits_ext[ext], bin_width, alpha=sub_bin_fit_alpha)
            ax.set_xlabel(xlabel)
            ax.set_ylabel('N')
            if yscale!='linear':
                ax.set_yscale(yscale)
            if xlim_e!='none':
                ax.set_xlim(xlim_e)
            if ylim_e != 'none':
                ax.set_ylim(ylim_e)
            else:
                auto_ylim = _auto_ylim_in_window(plot_centers, counts, xlim_e, yscale)
                if auto_ylim is not None:
                    ax.set_ylim(auto_ylim)

            sub_bin_peaks = sub_bin_peaks_ext[ext] if sub_bin_peaks_ext is not None else None

            # draw vertical lines and labels only at peaks within the displayed x-window
            if raw_peak_location_vertical_line or raw_peak_number or sub_bin_peak_location_vertical_line or sub_bin_peak_number:
                window = _xlim_window(xlim_e)
                for i,p in enumerate(peaks):
                    peak_x = centers[p]
                    if window is not None and not (window[0] <= peak_x <= window[1]):
                        continue

                    if raw_peak_location_vertical_line:
                        ax.axvline(peak_x, linestyle=linestyle, color=linecolor)
                    if raw_peak_number:
                        label_y = 0.93 if (i % 2 == 0) else 0.83
                        ax.text(peak_x,
                            label_y,
                            f"{i}",
                            transform=ax.get_xaxis_transform(),
                            verticalalignment='top',
                            horizontalalignment='center',
                            color=linecolor,
                            fontsize=peak_number_labels_size)

                    if sub_bin_peaks is not None:
                        sub_x = sub_bin_peaks[i]
                        if sub_bin_peak_location_vertical_line:
                            ax.axvline(sub_x, linestyle=sub_bin_linestyle, color=sub_bin_linecolor)
                        if sub_bin_peak_number:
                            sub_label_y = 0.73 if (i % 2 == 0) else 0.63
                            ax.text(sub_x,
                                sub_label_y,
                                f"{i}",
                                transform=ax.get_xaxis_transform(),
                                verticalalignment='top',
                                horizontalalignment='center',
                                color=sub_bin_linecolor,
                                fontsize=peak_number_labels_size)
            _peak_line_legend(ax,
                              raw_peak_location_vertical_line,
                              sub_bin_peaks is not None and sub_bin_peak_location_vertical_line,
                              linecolor, linestyle, sub_bin_linecolor, sub_bin_linestyle,
                              raw_bin_factor=peakfinder_bin_factor,
                              sub_bin_bin_factor=sub_bin_bin_factor
                              if sub_bin_bin_factor is not None else peakfinder_bin_factor,
                              display_bin_factor=display_bin_factor,
                              unit=bin_factor_unit,
                              gain=bin_factor_gains[ext]
                              if bin_factor_gains is not None else None,
                              fontsize=legend_fontsize)

            if save_plots:
                output_path = fig_name.with_stem(fig_name.stem + f'_EXT{ext+1}').with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            _finish_fig(show_plots)

    if plot_together:
        fig, axs = plt.subplots(2,2,figsize=subplots_figsize,constrained_layout=True,sharex=sharex,sharey=sharey)
        axs=axs.flatten()
        if show_titles:
            fig.suptitle(format_title(additional_title, suptitle, nimages))

        _auto_ylims = []
        for ext, peaks in enumerate(peaks_ext):
            centers=centers_ext[ext]
            ax = axs[ext]
            xlim_e = xlim[ext] if isinstance(xlim, list) else xlim
            ylim_e = ylim[ext] if isinstance(ylim, list) else ylim

            if data_ext is not None and window_bin_factor is not None:
                window = _xlim_window(xlim_e)
                if window is None:
                    raise ValueError(
                        f'plot_all_peaks: data_ext/window_bin_factor rebinning requires a '
                        f'concrete (left, right) xlim window, got xlim={xlim_e!r} for EXT {ext + 1}.'
                    )
                counts, plot_centers = _windowed_hist(data_ext[ext], window, window_bin_factor)
            else:
                counts = counts_ext[ext]
                plot_centers = centers
            bin_width = plot_centers[1] - plot_centers[0]

            _draw_hist(ax, plot_centers, counts, bin_width)
            if sub_bin_fits_ext is not None:
                _draw_sub_bin_fits(ax, sub_bin_fits_ext[ext], bin_width, alpha=sub_bin_fit_alpha)
            ax.set_xlabel(xlabel)
            ax.set_ylabel('N')
            if yscale!='linear':
                ax.set_yscale(yscale)
            if xlim_e!='none':
                ax.set_xlim(xlim_e)
            if ylim_e != 'none':
                ax.set_ylim(ylim_e)
            else:
                auto_ylim = _auto_ylim_in_window(plot_centers, counts, xlim_e, yscale)
                if auto_ylim is not None:
                    ax.set_ylim(auto_ylim)
                    _auto_ylims.append(auto_ylim)
            ax.set_title(f'EXT {ext + 1}')

            sub_bin_peaks = sub_bin_peaks_ext[ext] if sub_bin_peaks_ext is not None else None

            # draw vertical lines and labels only at peaks within the displayed x-window
            if raw_peak_location_vertical_line or raw_peak_number or sub_bin_peak_location_vertical_line or sub_bin_peak_number:
                window = _xlim_window(xlim_e)
                for i,p in enumerate(peaks):
                    peak_x = centers[p]
                    if window is not None and not (window[0] <= peak_x <= window[1]):
                        continue

                    if raw_peak_location_vertical_line:
                        ax.axvline(peak_x, linestyle=linestyle, color=linecolor)
                    if raw_peak_number:
                        label_y = 0.93 if (i % 2 == 0) else 0.83
                        ax.text(peak_x,
                            label_y,
                            f"{i}",
                            transform=ax.get_xaxis_transform(),
                            verticalalignment='top',
                            horizontalalignment='center',
                            color=linecolor,
                            fontsize=peak_number_labels_size)

                    if sub_bin_peaks is not None:
                        sub_x = sub_bin_peaks[i]
                        if sub_bin_peak_location_vertical_line:
                            ax.axvline(sub_x, linestyle=sub_bin_linestyle, color=sub_bin_linecolor)
                        if sub_bin_peak_number:
                            sub_label_y = 0.73 if (i % 2 == 0) else 0.63
                            ax.text(sub_x,
                                sub_label_y,
                                f"{i}",
                                transform=ax.get_xaxis_transform(),
                                verticalalignment='top',
                                horizontalalignment='center',
                                color=sub_bin_linecolor,
                                fontsize=peak_number_labels_size)
            _peak_line_legend(ax,
                              raw_peak_location_vertical_line,
                              sub_bin_peaks is not None and sub_bin_peak_location_vertical_line,
                              linecolor, linestyle, sub_bin_linecolor, sub_bin_linestyle,
                              raw_bin_factor=peakfinder_bin_factor,
                              sub_bin_bin_factor=sub_bin_bin_factor
                              if sub_bin_bin_factor is not None else peakfinder_bin_factor,
                              display_bin_factor=display_bin_factor,
                              unit=bin_factor_unit,
                              gain=bin_factor_gains[ext]
                              if bin_factor_gains is not None else None,
                              fontsize=legend_fontsize)

        # With sharey=True the per-axis set_ylim above lets the last extension dictate
        # the shared range; when auto-zooming to 'none', widen to span every extension's
        # window so no extension's peaks get clipped.
        if sharey and _auto_ylims:
            axs[0].set_ylim(min(b for b, _ in _auto_ylims), max(t for _, t in _auto_ylims))

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
        _finish_fig(show_plots)


def _resolve_nonlinearity_xbounds(xlim_e, peak_charge_e, fit_range_right):
    """Resolve an extension's x-limits to concrete (left, right) bounds.

    Mirrors how plot_nonlinearity resolves xlim for display, so the auto ylim is
    measured over exactly the x-range that will be shown.
    """
    if xlim_e == 'default':
        return -100.0, float(fit_range_right) + 500.0
    if xlim_e == 'none':
        x = np.asarray(peak_charge_e, dtype=float)
        return float(np.min(x)), float(np.max(x))
    left, right = xlim_e
    return float(left), float(right)


def _auto_nonlinearity_ylim(peak_charge_e, charge_minus_npeak, xbounds, pad=10.0):
    """[min - pad, max + pad] of the y values whose x falls within xbounds."""
    x = np.asarray(peak_charge_e, dtype=float)
    y = np.asarray(charge_minus_npeak, dtype=float)
    left, right = xbounds
    mask = (x >= left) & (x <= right)
    yv = y[mask] if np.any(mask) else y
    return float(np.min(yv)) - pad, float(np.max(yv)) + pad


def plot_nonlinearity(peaks_ext,
                      parabola_coeffs,
                      peak_charge_e_ext,
                      charge_minus_npeak_ext,
                      fit_range_right_ext,
                      xlim='default', ylim='default',
                      individual_figsize=(6,5), subplots_figsize=(13,9),
                      additional_title='',
                      suptitle='Nonlinearity Curve',
                      nimages=10,
                      line_color='r', 
                      scatter_color='b', 
                      s=6, 
                      alpha=0.5,
                      known_points=None,
                      parabola_coeffs_no_known=None,
                      known_colors=('orange', 'green', 'purple', 'brown', 'magenta', 'cyan'),
                      overlay=None, overlay_names=('Original', 'NL-corrected'),
                      overlay_scatter_color='purple', overlay_line_color='orange',
                      plot_individual=False,
                      plot_together=False,
                      sharex=True,
                      sharey=True,
                      show_titles=True,
                      save_plots=False,
                      show_plots=True,
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

    known_points_ext = resolve_known_points(known_points, len(peaks_ext))

    def _known_xy(ext):
        kp = known_points_ext[ext]
        return np.asarray([] if kp is None else kp, dtype=float).reshape(-1, 2).T

    def _plot_no_known(ax, ext):
        if parabola_coeffs_no_known is not None:
            xs = peak_charge_e_ext[ext]
            ax.plot(xs, parabola(xs, *parabola_coeffs_no_known[ext]), 'k--',
                    label='no known points: ' + _parabola_label(parabola_coeffs_no_known[ext]))

    def _plot_known(ax, ext):
        for i, (x, y) in enumerate(zip(*_known_xy(ext))):
            ax.scatter(x, y, c=known_colors[i % len(known_colors)], s=30, zorder=5, label=f'known NL({x:g}e-) = {y:.2f}e-')

    def _plot_data(ax, ext):
        """The nonlinearity points and fitted curve; with ``overlay`` (a second dataset as
        (peak_charge_e_ext, charge_minus_npeak_ext, parabola_coeffs)) both are drawn on
        the same axes, each point set and curve named in the legend."""
        datasets = [(peak_charge_e_ext, charge_minus_npeak_ext, parabola_coeffs, scatter_color, line_color,
                     overlay_names[0] if overlay is not None else None)]
        if overlay is not None:
            datasets.append((*overlay, overlay_scatter_color, overlay_line_color, overlay_names[1]))
        for xs, ys, coeffs, sc, lc, name in datasets:
            fit_label = _parabola_label(coeffs[ext])
            ax.plot(xs[ext], parabola(xs[ext], *coeffs[ext]), color=lc,
                    label=fit_label if name is None else f'{name} fit: {fit_label}')
            ax.scatter(xs[ext], ys[ext], c=sc, s=s, alpha=alpha, label=None if name is None else f'{name} data')

    # Precompute auto ylims: per extension, [min-10, max+10] of the nonlinearity
    # points within the displayed xlim. With sharey, every extension uses the
    # global min/max (which need not come from the same extension).
    auto_ylims = None
    if ylim == 'auto':
        auto_ylims = []
        for ext in range(len(peak_charge_e_ext)):
            xlim_e = xlim[ext] if isinstance(xlim, list) else xlim
            xbounds = _resolve_nonlinearity_xbounds(
                xlim_e, peak_charge_e_ext[ext], fit_range_right_ext[ext])
            kx, ky = _known_xy(ext)
            auto_ylims.append(_auto_nonlinearity_ylim(
                np.concatenate([peak_charge_e_ext[ext], kx]), np.concatenate([charge_minus_npeak_ext[ext], ky]), xbounds))
        if sharey:
            g_lo = min(lo for lo, _ in auto_ylims)
            g_hi = max(hi for _, hi in auto_ylims)
            auto_ylims = [(g_lo, g_hi)] * len(auto_ylims)

    if plot_individual:
        for ext, peaks in enumerate(peaks_ext):
            fig, ax = plt.subplots(1, 1, figsize=individual_figsize, constrained_layout=True)
            if show_titles:
                ax.set_title(format_title(additional_title, suptitle, nimages, f'EXT {ext + 1}'), pad=10)
            ax.grid()

            peak_charge_e=peak_charge_e_ext[ext]
            charge_minus_npeak=charge_minus_npeak_ext[ext]
            fit_range_right=fit_range_right_ext[ext]
            xlim_e = xlim[ext] if isinstance(xlim, list) else xlim
            ylim_e = ylim[ext] if isinstance(ylim, list) else ylim

            _plot_data(ax, ext)
            _plot_no_known(ax, ext)
            _plot_known(ax, ext)
            ax.legend(loc="upper right")
            ax.set_xlabel(r'Measured Pixel Charge ($e^-$)')
            ax.set_ylabel(r'Measured Pixel Charge - Peak n. ($e^-$) ')

            if ylim_e=='auto':
                ax.set_ylim(auto_ylims[ext])
            elif ylim_e=='default':
                fit_idx = int(np.searchsorted(peak_charge_e, fit_range_right))
                ky = _known_xy(ext)[1]
                ax.set_ylim(min([*charge_minus_npeak, *ky])-10,
                            max([*charge_minus_npeak[:fit_idx], *ky])+5)
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
            _finish_fig(show_plots)

    if plot_together:
        fig, axs = plt.subplots(2, 2, figsize=subplots_figsize, constrained_layout=True, sharex=sharex, sharey=sharey)
        axs=axs.flatten()
        if show_titles:
            fig.suptitle(format_title(additional_title, suptitle, nimages))
        for ext, peak_charge_e in enumerate(peak_charge_e_ext):
            ax = axs[ext]
            ax.grid()

            charge_minus_npeak=charge_minus_npeak_ext[ext]
            fit_range_right=fit_range_right_ext[ext]
            xlim_e = xlim[ext] if isinstance(xlim, list) else xlim
            ylim_e = ylim[ext] if isinstance(ylim, list) else ylim

            _plot_data(ax, ext)
            _plot_no_known(ax, ext)
            _plot_known(ax, ext)
            ax.legend(loc="upper right")
            ax.set_title(f'EXT {ext + 1}')
            ax.set_xlabel(r'Measured Pixel Charge ($e^-$)')
            ax.set_ylabel(r'Measured Pixel Charge - Peak n. ($e^-$) ')

            if ylim_e=='auto':
                ax.set_ylim(auto_ylims[ext])
            elif ylim_e=='default':
                fit_idx = int(np.searchsorted(peak_charge_e, fit_range_right))
                ky = _known_xy(ext)[1]
                ax.set_ylim(min([*charge_minus_npeak, *ky])-10,
                            max([*charge_minus_npeak[:fit_idx], *ky])+5)
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
        _finish_fig(show_plots)


def plot_images(panels, panel_labels=None, stack='horizontal',
                vmin=None, vmax=None, norm=LogNorm, norm_kwargs=None, cmap='gray', cbar_label='Pixel charge (ADU)',
                individual_figsize=(10, 5),
                additional_title='',
                suptitle='Image',
                nimages=10,
                show_titles=True,
                save_plots=False,
                show_plots=True,
                fig_path='./', file='image', suffix='image',
                dpi=350):
    """One figure per extension showing 2-D pixel arrays on a log greyscale with a colorbar.

    ``panels`` is a list of images (each a per-extension list of 2-D arrays). Every
    panel of an extension is drawn side by side (``stack='horizontal'``) or one above the
    other (``'vertical'``), sharing one norm and one colorbar so they are directly comparable.

    ``norm`` is a matplotlib normalization class (``LogNorm``, the default, or e.g.
    ``SymLogNorm`` -- which also shows values <= 0, with ``norm_kwargs=dict(linthresh=1)``
    -- or ``PowerNorm`` with ``norm_kwargs=dict(gamma=0.3)``). Unless ``vmin``/``vmax``
    are given, it spans the 1st-99.9th percentile of the extension's pixels across its
    panels (the positive ones only for ``LogNorm``, which cannot show anything <= 0).
    """
    labels = panel_labels if panel_labels is not None else [None] * len(panels)
    for ext in range(len(panels[0])):
        images = [np.asarray(panel[ext]) for panel in panels]
        pixels = np.concatenate([im[im > 0] if norm is LogNorm else im.ravel() for im in images])
        panel_norm = norm(vmin=vmin if vmin is not None else np.percentile(pixels, 1),
                          vmax=vmax if vmax is not None else np.percentile(pixels, 99.9),
                          **(norm_kwargs or {}))

        shape = (1, len(images)) if stack == 'horizontal' else (len(images), 1)
        fig, axs = plt.subplots(*shape, figsize=individual_figsize, constrained_layout=True,
                                sharex=True, sharey=True, squeeze=False)
        axs = axs.ravel()
        for ax, image, label in zip(axs, images, labels):
            # matplotlib's default antialiasing: 'nearest' would sample only ~1 column in 8 of a
            # 3500-column image, dropping thin features such as tracks
            mappable = ax.imshow(image, cmap=cmap, norm=panel_norm, aspect='auto')
            ax.set_xlabel('Column')
            ax.set_ylabel('Row')
            if label:
                ax.set_title(label)
        fig.colorbar(mappable, ax=axs.tolist(), label=cbar_label)
        if show_titles:
            title = format_title(additional_title, suptitle, nimages, f'EXT {ext + 1}')
            if len(images) > 1:
                fig.suptitle(title)
            else:
                axs[0].set_title(title, pad=10)

        if save_plots:
            output_path = Path(fig_path) / f'{Path(file).stem}_{suffix}_EXT{ext + 1}.jpeg'
            plt.savefig(str(output_path), dpi=dpi)
            print(f'Saved plot to {output_path}')
        _finish_fig(show_plots)


def _draw_resolution_window(ax, res, sigma_well=0.25, sigma_limit=0.5):
    """Draw one resolution window (data bars, comb components, comb sum, null,
    and a verdict annotation) onto an existing Axes."""
    xw, yw = res['xw'], res['yw']
    bw = res['bin_width']
    comb, _ = _make_comb_model(res['means'])
    xf = np.linspace(res['lo'], res['hi'], 2000)

    ax.bar(xw, yw, width=bw, align='center', color='0.8', edgecolor='none',
           label='data', zorder=1)

    sigma = res['comb_popt'][0]
    amps = res['comb_popt'][1:]
    for i, (amp, mu) in enumerate(zip(amps, res['means'])):
        ax.plot(xf, _gauss_single(xf, amp, mu, sigma), color='tab:blue', lw=0.8,
                alpha=0.7, zorder=2, label='components' if i == 0 else None)
    ax.plot(xf, comb(xf, *res['comb_popt']), color='tab:red', lw=2.0, zorder=4,
            label='n-Gaussian fit')
    if res['null_popt'] is not None:
        ax.plot(xf, _gauss_single(xf, *res['null_popt']), color='k', lw=1.3,
                ls='--', zorder=3, label='single-Gaussian null')

    verdict = res.get('verdict') or classify_resolution(res, sigma_well, sigma_limit)
    color = _RESOLUTION_VERDICT_COLORS.get(verdict, 'k')
    ax.set_xlabel(r'Charge ($e^-$)')
    ax.set_ylabel('N')
    txt = (rf'$\sigma$ = {res["sigma_e"]:.3f} $\pm$ {res["sigma_e_err"]:.3f} $e^-$'+'\n'
           f'peaks fit: {res["n_components"]}')
    ax.text(0.02, 0.97, txt, transform=ax.transAxes, va='top', ha='left',
            fontsize=plt.rcParams['legend.fontsize'], zorder=11,
            bbox=dict(boxstyle='square', fc='white', ec=color, lw=1.5, alpha=0.85))
    # Draw the legend last and force it above all data bars, fit curves, and the
    # annotation box so it is never occluded.
    leg = ax.legend(loc='upper right')
    leg.set_zorder(12)
    return verdict


def plot_resolution(results_ext, charges,
                    individual_figsize=(8, 6.5), subplots_figsize=(13, 9),
                    additional_title='',
                    suptitle='Single-Electron Resolution',
                    nimages=10,
                    sigma_well=0.25, sigma_limit=0.5,
                    plot_individual=False, plot_together=True,
                    sharex=False, sharey=False,
                    show_titles=True,
                    save_plots=False, show_plots=True,
                    fig_path='./', file='resolution', dpi=350):
    """Plot the resolution windows from resolution_at_charge_ext.

    ``plot_individual`` makes one figure per (extension, charge). ``plot_together``
    makes one 2x2 figure per charge (one subplot per extension).
    """
    charges = _normalize_charges(charges)
    fig_path = Path(fig_path)
    if file != 'resolution':
        base_name = file[:-5] + '_resolution'
    else:
        base_name = file
    fig_name = fig_path / base_name
    n_ext = len(results_ext)

    if plot_individual:
        for rows in results_ext:
            for res in rows:
                if res is None:
                    continue
                fig, ax = plt.subplots(figsize=individual_figsize,
                                       constrained_layout=True)
                _draw_resolution_window(ax, res, sigma_well, sigma_limit)
                if show_titles:
                    ax.set_title(
                        format_title(additional_title, suptitle, nimages, f'EXT {res["ext"]}'),
                        pad=10)
                if save_plots:
                    output_path = fig_name.with_stem(
                        fig_name.stem + f'_EXT{res["ext"]}_q{res["q"]:.0f}'
                    ).with_suffix('.jpeg')
                    plt.savefig(str(output_path), dpi=dpi)
                    print(f'Saved plot to {output_path}')
                _finish_fig(show_plots)

    if plot_together:
        for j, q in enumerate(charges):
            present = [results_ext[e][j] for e in range(n_ext)
                       if j < len(results_ext[e]) and results_ext[e][j] is not None]
            if not present:
                continue
            fig, axs = plt.subplots(2, 2, figsize=subplots_figsize,
                                    constrained_layout=True, sharex=sharex,
                                    sharey=sharey)
            axs = axs.flatten()
            if show_titles:
                fig.suptitle(format_title(additional_title, suptitle, nimages))
            for e in range(len(axs)):
                res = results_ext[e][j] if e < n_ext and j < len(results_ext[e]) else None
                if res is None:
                    axs[e].set_visible(False)
                    continue
                _draw_resolution_window(axs[e], res, sigma_well, sigma_limit)
                if show_titles:
                    axs[e].set_title(f'EXT {res["ext"]}')
            if save_plots:
                output_path = fig_name.with_stem(
                    fig_name.stem + f'_q{q:.0f}'
                ).with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            _finish_fig(show_plots)
