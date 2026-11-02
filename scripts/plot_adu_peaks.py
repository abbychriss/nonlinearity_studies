#!/usr/bin/env python3
"""ADU-spectrum peak plot, in the same format as nonlinearity_studies.plot_all_peaks.

Companion to gain_independent_nonlinearity.py. Draws the pedestal-subtracted
charge spectrum in ADU with the same decorations as the production
``plot_all_peaks``: bar histogram, dashed vertical lines at each detected peak,
alternating peak-number labels, per-extension individual figures and/or a 2x2
"together" grid, the same auto-ylim-in-window behaviour, the same optional
re-binning of a zoomed window, and the same .jpeg saving convention.

Three deliberate differences from the production version, all of them the point
of the exercise rather than stylistic drift:

  1. x-axis is ADU, not electrons. No gain is applied anywhere.

  2. Vertical lines are drawn at the SUB-BIN peak positions (mu), not at the
     bin center ``centers[p]``. On this data the bin width is ~0.05 ADU while
     the effect under study is ~0.001 ADU/e-, so bin-center quantization would
     be larger than the signal.

  3. Labels show the ASSIGNED ELECTRON NUMBER, not the ordinal index. The
     production version labels with ``f"{i}"``, the position in the peak list,
     which silently mislabels every peak after a missed one. Passing ``n_ext``
     from the gain-independent analysis labels with the real electron number
     and lets gaps be drawn distinctly (see ``gap_ext``).

Typical use, from a notebook or the CLI below::

    from plot_adu_peaks import plot_all_peaks_adu
    plot_all_peaks_adu(counts_ext, mu_ext, n_ext, centers_ext,
                       xlim=(0, 40), xlim_units='electrons', gain_ext=gains,
                       plot_together=True, show_plots=True)
"""
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:                                    # same helper the production plots use
    from pedestal_subtract.core import _finish_fig
except ImportError:                     # keep this module importable standalone
    def _finish_fig(show_plots):
        if show_plots:
            plt.show()
        plt.close()


def _xlim_window(xlim_e):
    """Resolve xlim to a (left, right) tuple, or None if it is not a concrete window."""
    is_window = (
        isinstance(xlim_e, (list, tuple)) and len(xlim_e) == 2
        and all(isinstance(v, (int, float, np.integer, np.floating)) for v in xlim_e)
    )
    if not is_window:
        return None
    return min(xlim_e), max(xlim_e)


def _auto_ylim_in_window(centers, counts, xlim_e, yscale='linear', pad=0.05):
    """Tight y-limits framing only the histogram bars inside the x-window.

    Same behaviour as the production helper: without this, autoscale is
    dominated by the enormous pedestal peak and every other peak is invisible.
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


def _windowed_hist_adu(data, window, bin_factor, gain):
    """Re-histogram ``data`` over ``window`` (ADU) at ``bin_factor`` bins per ELECTRON.

    ``bin_factor`` keeps the production meaning -- bins per electron, not bins
    per ADU -- so the same number gives the same visual resolution as it does
    in the electron-space plot. Converting requires the gain, which is used
    here for display binning ONLY and never enters any measured quantity.
    """
    if gain is None or not np.isfinite(gain) or gain <= 0:
        raise ValueError('plot_all_peaks_adu: window_bin_factor rebinning needs a '
                         'finite positive gain_ext for the extension being plotted.')
    span_e = (window[1] - window[0]) / float(gain)
    bins = max(1, math.floor(span_e * bin_factor))
    data = np.asarray(data)
    data_window = data[(data > window[0]) & (data < window[1])]
    counts, edges = np.histogram(data_window, bins=bins, range=window)
    centers = 0.5 * (edges[1:] + edges[:-1])
    return counts, centers


def _resolve_window(xlim_e, xlim_units, gain, pedestal=0.0):
    """Return the (left, right) window in ADU, converting from electrons if asked.

    In ADU the peaks sit ~1 gain apart (~1.3 ADU here), so a window that is
    natural to think about in electrons ("show me peaks 500 to 510") is an
    awkward ADU number. ``xlim_units='electrons'`` lets the caller keep
    thinking in electrons while everything plotted and measured stays in ADU.
    """
    window = _xlim_window(xlim_e)
    if window is None or xlim_units == 'ADU':
        return window
    if xlim_units != 'electrons':
        raise ValueError(f"xlim_units must be 'ADU' or 'electrons', got {xlim_units!r}")
    if gain is None or not np.isfinite(gain) or gain <= 0:
        raise ValueError("plot_all_peaks_adu: xlim_units='electrons' needs a finite "
                         "positive gain_ext for the extension being plotted.")
    return (pedestal + window[0] * gain, pedestal + window[1] * gain)


def _draw_peak_decorations(ax, mu, n_lab, gap, window, linecolor, linestyle,
                           gap_linecolor, gap_linestyle, show_labels, label_size):
    """Vertical lines + alternating peak-number labels, for peaks inside the window.

    Label parity keys off the electron number rather than the loop index, so a
    peak sits at the same height whatever window it is displayed in.
    """
    for j in range(len(mu)):
        peak_x = mu[j]
        if window is not None and not (window[0] <= peak_x <= window[1]):
            continue

        is_gap = bool(gap[j]) if gap is not None else False
        color = gap_linecolor if is_gap else linecolor
        ax.axvline(peak_x, linestyle=(gap_linestyle if is_gap else linestyle),
                   color=color)

        if show_labels:
            lab = int(n_lab[j])
            label_y = 0.93 if (lab % 2 == 0) else 0.83
            ax.text(peak_x, label_y, f'{lab}',
                    transform=ax.get_xaxis_transform(),
                    verticalalignment='top', horizontalalignment='center',
                    color=color, fontsize=label_size)


def _rebin(counts, centers, k):
    """Sum groups of k adjacent bins. Exact -- no re-reading of the raw data.

    The ANALYSIS histogram has to be fine (~4 bins per readout-noise sigma) or
    sub-bin localization has nothing to work with. That is far finer than is
    useful to look at: at 25 bins/e- the bars are dominated by single-count
    Poisson spikes and the peak shape is hard to see. Rebinning for DISPLAY
    only keeps the two requirements separate.
    """
    k = int(k)
    if k <= 1:
        return np.asarray(counts, dtype=float), np.asarray(centers, dtype=float)
    counts = np.asarray(counts, dtype=float)
    centers = np.asarray(centers, dtype=float)
    n_keep = (counts.size // k) * k
    c = counts[:n_keep].reshape(-1, k).sum(axis=1)
    x = centers[:n_keep].reshape(-1, k).mean(axis=1)
    return c, x


def plot_all_peaks_adu(counts_ext,
                       mu_ext,
                       n_ext,
                       centers_ext,
                       xlim=(1000, 1010),
                       xlim_units='electrons',
                       rebin=1,
                       ylim='none',
                       yscale='linear',
                       data_ext=None,
                       window_bin_factor=None,
                       gain_ext=None,
                       pedestal_ext=None,
                       gap_ext=None,
                       plot_individual=False, plot_together=False,
                       draw_lines=True, linecolor='r', linestyle='--',
                       gap_linecolor='tab:orange', gap_linestyle='-',
                       peak_number_labels_individual=True,
                       peak_number_labels_together=True,
                       peak_number_label_size=8,
                       individual_figsize=(7, 6), subplots_figsize=(13, 9),
                       additional_title='',
                       suptitle='Peaks in Pixel Charge Distribution (ADU)',
                       nimages=10,
                       sharex=True,
                       sharey=True,
                       show_titles=True,
                       save_plots=False,
                       show_plots=True,
                       fig_path='./', file='peak_finder_adu',
                       dpi=350):
    """Plot the pedestal-subtracted ADU spectrum with peak lines and numbers.

    Parameters mirror ``nonlinearity_studies.plotting.plot_all_peaks`` except
    where noted.

    counts_ext, centers_ext : per-extension histogram counts and bin centers,
        in ADU (as returned by the gain-independent analysis).
    mu_ext : per-extension SUB-BIN peak positions in ADU. This replaces the
        production ``peaks_ext`` (indices into ``centers``); pass positions,
        not indices.
    n_ext : per-extension assigned electron numbers, used for the labels.
    gap_ext : optional per-extension bool array, True where a peak follows a
        detected gap. Those lines are drawn in ``gap_linecolor`` with
        ``gap_linestyle`` so a missed peak is visible at a glance.
    xlim_units : 'electrons' (default) or 'ADU'. Electrons are converted with
        ``gain_ext`` for DISPLAY ONLY.
    gain_ext, pedestal_ext : per-extension scalars, needed only for
        ``xlim_units='electrons'`` and/or ``window_bin_factor`` rebinning.
    """
    fig_path = Path(fig_path)
    if file != 'peak_finder_adu':
        base_name = file[:-5] + '_peak_finder_adu' if file.endswith('.fits') \
            else file + '_peak_finder_adu'
    else:
        base_name = file
    fig_name = fig_path / base_name

    n_plots = len(mu_ext)

    def _per_ext(seq, ext, default=None):
        if seq is None:
            return default
        return seq[ext] if isinstance(seq, (list, tuple, np.ndarray)) else seq

    def _prepare(ext):
        """Resolve everything this extension needs to be drawn."""
        centers = np.asarray(centers_ext[ext], dtype=float)
        mu = np.asarray(mu_ext[ext], dtype=float)
        n_lab = np.asarray(n_ext[ext])
        gap = None if gap_ext is None else np.asarray(gap_ext[ext], dtype=bool)

        gain = _per_ext(gain_ext, ext)
        pedestal = _per_ext(pedestal_ext, ext, 0.0) or 0.0

        xlim_e = xlim[ext] if isinstance(xlim, list) else xlim
        ylim_e = ylim[ext] if isinstance(ylim, list) else ylim
        window = _resolve_window(xlim_e, xlim_units, gain, pedestal)

        if data_ext is not None and window_bin_factor is not None:
            if window is None:
                raise ValueError(
                    f'plot_all_peaks_adu: data_ext/window_bin_factor rebinning '
                    f'requires a concrete (left, right) xlim window, got '
                    f'xlim={xlim_e!r} for EXT {ext + 1}.')
            counts, plot_centers = _windowed_hist_adu(
                data_ext[ext], window, window_bin_factor, gain)
        else:
            counts = np.asarray(counts_ext[ext], dtype=float)
            plot_centers = centers
            k = rebin[ext] if isinstance(rebin, (list, tuple, np.ndarray)) else rebin
            counts, plot_centers = _rebin(counts, plot_centers, k)
        return counts, plot_centers, mu, n_lab, gap, window, ylim_e

    def _style(ax, counts, plot_centers, window, ylim_e):
        bin_width = plot_centers[1] - plot_centers[0]
        ax.bar(plot_centers, counts, align='center', edgecolor='none',
               linewidth=0, width=bin_width)
        ax.set_xlabel('Charge (ADU)')
        ax.set_ylabel('N')
        if yscale != 'linear':
            ax.set_yscale(yscale)
        if window is not None:
            ax.set_xlim(window)
        auto = None
        if ylim_e != 'none':
            ax.set_ylim(ylim_e)
        else:
            auto = _auto_ylim_in_window(plot_centers, counts, window, yscale)
            if auto is not None:
                ax.set_ylim(auto)
        return auto

    if plot_individual:
        for ext in range(n_plots):
            counts, plot_centers, mu, n_lab, gap, window, ylim_e = _prepare(ext)
            fig, ax = plt.subplots(1, 1, figsize=individual_figsize,
                                   constrained_layout=True)
            if show_titles:
                ax.set_title(
                    f'{additional_title}{suptitle} (Nimages = {nimages}): EXT {ext + 1}',
                    fontsize=12, pad=10)
            _style(ax, counts, plot_centers, window, ylim_e)
            if draw_lines:
                _draw_peak_decorations(
                    ax, mu, n_lab, gap, window, linecolor, linestyle,
                    gap_linecolor, gap_linestyle,
                    peak_number_labels_individual, peak_number_label_size)

            if save_plots:
                output_path = (fig_name.with_stem(fig_name.stem + f'_EXT{ext+1}')
                               .with_suffix('.jpeg'))
                output_path.parent.mkdir(parents=True, exist_ok=True)
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            _finish_fig(show_plots)

    if plot_together:
        nrow = 2 if n_plots > 2 else 1
        ncol = 2 if n_plots > 1 else 1

        # Resolve every extension up front so we can decide whether sharing the
        # x-axis is even legitimate.
        prepared = [_prepare(ext) for ext in range(n_plots)]
        windows = [p[5] for p in prepared]

        # CRITICAL DIFFERENCE FROM THE ELECTRON-SPACE PLOT: each extension has
        # its own gain, so a window specified in electrons maps to a DIFFERENT
        # ADU range per extension. Sharing the x-axis then forces every panel
        # into the last extension's range, and the other extensions' peaks are
        # drawn off-screen -- silently, since the vertical lines are simply not
        # visible. The production plot_all_peaks is immune because its x-axis is
        # already in electrons. Share only when the resolved windows agree.
        share_x = sharex
        if share_x and all(w is not None for w in windows) and len(windows) > 1:
            spans = [w[1] - w[0] for w in windows]
            tol = 0.01 * max(spans)
            if (max(w[0] for w in windows) - min(w[0] for w in windows) > tol or
                    max(w[1] for w in windows) - min(w[1] for w in windows) > tol):
                share_x = False
                lo, hi = min(w[0] for w in windows), max(w[1] for w in windows)
                print('plot_all_peaks_adu: per-extension ADU windows differ '
                      f'({lo:.4g} to {hi:.4g} ADU across extensions) because each '
                      'extension has its own gain; disabling sharex so every '
                      'panel shows its own peaks.')

        fig, axs = plt.subplots(nrow, ncol, figsize=subplots_figsize,
                                constrained_layout=True, sharex=share_x, sharey=sharey)
        axs = np.atleast_1d(axs).flatten()
        if show_titles:
            fig.suptitle(f'{additional_title}{suptitle} (Nimages = {nimages})')

        _auto_ylims = []
        for ext in range(n_plots):
            counts, plot_centers, mu, n_lab, gap, window, ylim_e = prepared[ext]
            ax = axs[ext]
            auto = _style(ax, counts, plot_centers, window, ylim_e)
            if auto is not None:
                _auto_ylims.append(auto)
            ax.set_title(f'EXT {ext + 1}')
            if draw_lines:
                _draw_peak_decorations(
                    ax, mu, n_lab, gap, window, linecolor, linestyle,
                    gap_linecolor, gap_linestyle,
                    peak_number_labels_together, peak_number_label_size)

        # With sharey=True the per-axis set_ylim lets the last extension dictate
        # the shared range; widen to span every extension's window so no
        # extension's peaks get clipped.
        if sharey and _auto_ylims:
            axs[0].set_ylim(min(b for b, _ in _auto_ylims),
                            max(t for _, t in _auto_ylims))

        for i in range(len(axs)):
            if i < len(axs) - ncol:
                axs[i].set_xlabel('')
            axs[i].tick_params(labelbottom=True)
            if i % ncol:
                axs[i].set_ylabel('')
            axs[i].tick_params(labelleft=True)
        for i in range(n_plots, len(axs)):
            axs[i].set_visible(False)

        if save_plots:
            output_path = fig_name.with_suffix('.jpeg')
            output_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(str(output_path), dpi=dpi)
            print(f'Saved plot to {output_path}')
        _finish_fig(show_plots)


def plot_from_results(results, xlim=(1000, 1010), xlim_units='electrons',
                      display_bin_factor=8.0, **kwargs):
    """Convenience wrapper: plot straight from gain_independent_nonlinearity results.

    ``results`` is the list of per-extension dicts returned by the analysis
    (each carrying its ``_arrays`` payload), so a notebook can go from analysis
    to plot without unpacking anything by hand.

    ``display_bin_factor`` is the bins-per-electron actually drawn. The analysis
    histogram is much finer than this by necessity (sub-bin localization needs
    ~4 bins per readout-noise sigma, which here is ~25 bins/e-); the extra
    resolution just shows up as single-count Poisson spikes when plotted. Bins
    are summed in exact integer groups to reach roughly this density, so set it
    to None to draw the analysis binning untouched.
    """
    counts_ext = [r['_arrays']['counts'] for r in results]
    centers_ext = [r['_arrays']['centers'] for r in results]
    n_ext = [r['_arrays']['n'] for r in results]
    gap_ext = [r['_arrays']['gap_before'] for r in results]
    # _arrays['Q'] is referenced to the zero peak; put it back on the raw ADU
    # axis the histogram uses by adding the zero-peak position back on.
    mu_ext = []
    for r in results:
        Q = r['_arrays']['Q']
        n = r['_arrays']['n']
        Q0_axis = r['_arrays'].get('Q0_axis')
        if Q0_axis is None:                     # fall back: locate n == 0
            Q0_axis = 0.0
        mu_ext.append(Q + Q0_axis)
    gain_ext = [r['q1_ADU'] for r in results]
    pedestal_ext = [r['_arrays'].get('Q0_axis', 0.0) for r in results]

    if display_bin_factor is None:
        rebin = 1
    else:
        rebin = [max(1, int(round((r['q1_ADU'] / r['bin_width_ADU'])
                                  / display_bin_factor))) for r in results]

    return plot_all_peaks_adu(counts_ext, mu_ext, n_ext, centers_ext,
                              xlim=xlim, xlim_units=xlim_units, rebin=rebin,
                              gain_ext=gain_ext, pedestal_ext=pedestal_ext,
                              gap_ext=gap_ext, **kwargs)


def _main(argv=None):
    """Plot directly from the CSVs written by gain_independent_nonlinearity.py."""
    import argparse
    p = argparse.ArgumentParser(
        description='Plot the ADU charge spectrum with peak lines and electron-number '
                    'labels, in the plot_all_peaks format.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('fits_file', help='The FITS image the analysis was run on.')
    p.add_argument('--results_dir', required=True,
                   help='Directory containing the *_extN_peaks.csv files.')
    p.add_argument('--ext', type=int, nargs='+', default=None,
                   help='1-based extensions to plot. Default: all found.')
    p.add_argument('--xlim', type=float, nargs=2, default=[0, 40])
    p.add_argument('--xlim_units', choices=['electrons', 'ADU'], default='electrons')
    p.add_argument('--yscale', default='linear', choices=['linear', 'log'])
    p.add_argument('--window_bin_factor', type=float, default=None,
                   help='Re-bin the displayed window at this many bins per electron.')
    p.add_argument('--pedsub_cache_dir', default=None)
    p.add_argument('--pedestal_subtraction_axis', default='row')
    p.add_argument('--n_std_to_mask', type=float, default=1.5)
    p.add_argument('--together', action='store_true')
    p.add_argument('--individual', action='store_true')
    p.add_argument('--save', action='store_true')
    p.add_argument('--no_show', action='store_true')
    p.add_argument('--fig_path', default='./')
    args = p.parse_args(argv)

    from pedestal_subtract.core import get_fits, pedestal_subtract_ext_cached

    fits_path = Path(args.fits_file).expanduser().resolve()
    stem = fits_path.name.split('.')[0]
    rdir = Path(args.results_dir)

    found = sorted(rdir.glob(f'{stem}_ext*_peaks.csv'))
    if not found:
        raise SystemExit(f'No {stem}_ext*_peaks.csv files in {rdir}')
    exts = ([e for e in args.ext] if args.ext
            else [int(f.stem.split('_ext')[1].split('_')[0]) for f in found])

    data_ext = get_fits(str(fits_path))
    data_ext = pedestal_subtract_ext_cached(
        data_ext, source_path=fits_path, n_std_to_mask=args.n_std_to_mask,
        axis=args.pedestal_subtraction_axis, cache_dir=args.pedsub_cache_dir,
        verbose=True)

    counts_ext, centers_ext, mu_ext, n_ext, gap_ext, gain_ext, dsel = \
        [], [], [], [], [], [], []
    for e in exts:
        d = np.genfromtxt(rdir / f'{stem}_ext{e}_peaks.csv', delimiter=',', names=True)
        n = d['n']
        Q = d['Q_ADU']
        q1 = float(Q[n == 1][0])
        flat = np.asarray(data_ext[e - 1]).ravel()
        flat = flat[np.isfinite(flat)]
        # Q is referenced to the zero peak; the histogram axis is raw ADU, so
        # recover the zero-peak position as the median of the pedestal core.
        Q0 = float(np.median(flat[np.abs(flat - np.median(flat)) < 0.5 * q1]))
        lo, hi = (args.xlim if args.xlim_units == 'ADU'
                  else (Q0 + args.xlim[0] * q1, Q0 + args.xlim[1] * q1))
        bw = q1 / (args.window_bin_factor or 25.0)
        c, edges = np.histogram(flat, bins=max(10, int((hi - lo) / bw)), range=(lo, hi))
        counts_ext.append(c)
        centers_ext.append(0.5 * (edges[1:] + edges[:-1]))
        mu_ext.append(Q + Q0)
        n_ext.append(n)
        gap_ext.append(d['gap_before'].astype(bool))
        gain_ext.append(q1)
        dsel.append(flat)

    plot_all_peaks_adu(
        counts_ext, mu_ext, n_ext, centers_ext,
        xlim=tuple(args.xlim), xlim_units=args.xlim_units, yscale=args.yscale,
        gain_ext=gain_ext,
        pedestal_ext=[float(np.median(d[np.abs(d - np.median(d)) < 0.5 * g]))
                      for d, g in zip(dsel, gain_ext)],
        gap_ext=gap_ext,
        plot_together=args.together or not args.individual,
        plot_individual=args.individual,
        save_plots=args.save, show_plots=not args.no_show,
        fig_path=args.fig_path, file=fits_path.name)


if __name__ == '__main__':
    _main()
