"""
plotting — split out of nonlinearity_studies.py.
"""
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from pedestal_subtract.core import _finish_fig

from .fit_data import _gauss_single, _make_comb_model, parabola
from .resolution import classify_resolution
from .summary import _normalize_charges

_RESOLUTION_VERDICT_COLORS = {
    'well resolved': 'tab:green', 'marginal': 'tab:orange', 'unresolved': 'tab:red',
}


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


def plot_all_peaks(counts_ext,
                   peaks_ext,
                   centers_ext,
                   xlim=(500,510),
                   ylim='none',
                   yscale='linear', 
                   plot_individual=False, plot_together=False,
                   draw_lines=True, linecolor='r', linestyle='--',
                   peak_number_labels_individual=True, peak_number_labels_together=True, peak_number_label_size=8,
                   individual_figsize=(7,6), subplots_figsize=(13,9),
                   additional_title='',
                   suptitle='Peaks in Pixel Charge Distribution',
                   nimages=10,
                   sharex=True,
                   sharey=True,
                   show_titles=True,
                   save_plots=False,
                   show_plots=True,
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
                ax.set_title(f'{additional_title}{suptitle} (Nimages = {nimages}): EXT {ext + 1}', fontsize=12, pad=10)
            ax.bar(centers, counts, align='center', edgecolor='none', linewidth=0, width=bin_width)
            ax.set_xlabel(r'Charge ($e^-$)')
            ax.set_ylabel('N')
            if yscale!='linear':
                ax.set_yscale(yscale)
            xlim_e = xlim[ext] if isinstance(xlim, list) else xlim
            ylim_e = ylim[ext] if isinstance(ylim, list) else ylim
            if xlim_e!='none':
                ax.set_xlim(xlim_e)
            if ylim_e != 'none':
                ax.set_ylim(ylim_e)
            else:
                auto_ylim = _auto_ylim_in_window(centers, counts, xlim_e, yscale)
                if auto_ylim is not None:
                    ax.set_ylim(auto_ylim)

            # draw vertical lines and labels only at peaks within the displayed x-window
            if draw_lines:
                window = _xlim_window(xlim_e)
                for i,p in enumerate(peaks):
                    peak_x = centers[p]
                    if window is not None and not (window[0] <= peak_x <= window[1]):
                        continue

                    ax.axvline(peak_x, linestyle=linestyle, color=linecolor)
                    if peak_number_labels_individual:
                        label_y = 0.93 if (i % 2 == 0) else 0.83
                        ax.text(peak_x,
                            label_y,
                            f"{i}",
                            transform=ax.get_xaxis_transform(),
                            verticalalignment='top',
                            horizontalalignment='center',
                            color=linecolor,
                            fontsize=peak_number_label_size)

            if save_plots:
                output_path = fig_name.with_stem(fig_name.stem + f'_EXT{ext+1}').with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            _finish_fig(show_plots)

    if plot_together:
        fig, axs = plt.subplots(2,2,figsize=subplots_figsize,constrained_layout=True,sharex=sharex,sharey=sharey)
        axs=axs.flatten()
        if show_titles:
            fig.suptitle(f'{additional_title}{suptitle} (Nimages = {nimages})')

        _auto_ylims = []
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
            xlim_e = xlim[ext] if isinstance(xlim, list) else xlim
            ylim_e = ylim[ext] if isinstance(ylim, list) else ylim
            if xlim_e!='none':
                ax.set_xlim(xlim_e)
            if ylim_e != 'none':
                ax.set_ylim(ylim_e)
            else:
                auto_ylim = _auto_ylim_in_window(centers, counts, xlim_e, yscale)
                if auto_ylim is not None:
                    ax.set_ylim(auto_ylim)
                    _auto_ylims.append(auto_ylim)
            ax.set_title(f'EXT {ext + 1}')

            # draw vertical lines and labels only at peaks within the displayed x-window
            if draw_lines:
                window = _xlim_window(xlim_e)
                for i,p in enumerate(peaks):
                    peak_x = centers[p]
                    if window is not None and not (window[0] <= peak_x <= window[1]):
                        continue

                    ax.axvline(peak_x, linestyle=linestyle, color=linecolor)
                    if peak_number_labels_together:
                        label_y = 0.93 if (i % 2 == 0) else 0.83
                        ax.text(peak_x,
                            label_y,
                            f"{i}",
                            transform=ax.get_xaxis_transform(),
                            verticalalignment='top',
                            horizontalalignment='center',
                            color=linecolor,
                            fontsize=peak_number_label_size)

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
                      suptitle='Pixel Charge Nonlinearity',
                      nimages=10,
                      line_color='r', 
                      scatter_color='b', 
                      s=2, 
                      alpha=0.5,
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
            auto_ylims.append(_auto_nonlinearity_ylim(
                peak_charge_e_ext[ext], charge_minus_npeak_ext[ext], xbounds))
        if sharey:
            g_lo = min(lo for lo, _ in auto_ylims)
            g_hi = max(hi for _, hi in auto_ylims)
            auto_ylims = [(g_lo, g_hi)] * len(auto_ylims)

    if plot_individual:
        for ext, peaks in enumerate(peaks_ext):
            fig, ax = plt.subplots(1, 1, figsize=individual_figsize, constrained_layout=True)
            if show_titles:
                ax.set_title(f'{additional_title}{suptitle} (Nimages = {nimages}): EXT {ext + 1}', fontsize=12, pad=10)
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
            ax.legend(loc="upper right", fontsize=10)
            ax.set_xlabel(r'Measured Pixel Charge ($e^-$)')
            ax.set_ylabel(r'Measured Pixel Charge - Peak n. ($e^-$) ')

            if ylim_e=='auto':
                ax.set_ylim(auto_ylims[ext])
            elif ylim_e=='default':
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
            _finish_fig(show_plots)

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
            ax.set_title(f'EXT {ext + 1}')
            ax.set_xlabel(r'Measured Pixel Charge ($e^-$)')
            ax.set_ylabel(r'Measured Pixel Charge - Peak n. ($e^-$) ')

            if ylim_e=='auto':
                ax.set_ylim(auto_ylims[ext])
            elif ylim_e=='default':
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
    txt = (rf'$\sigma$ = {res["sigma_e"]:.3f} $\pm$ {res["sigma_e_err"]:.3f} $e^-$ ($\sigma_0$ = {res["sigma_seed_e"]:.3f})'+'\n'
           f'peaks fit: {res["n_components"]} of {res["expected_peaks"]} expected\n'
           rf'reduced $\chi^2$ = {res["reduced_chi2"]:.2f}'
           f'   $\\Delta$AIC = {res["delta_aic"]:.0f}\n'
           f'verdict: {verdict.upper()}'
           rf'  ($\sigma$<{sigma_well} resolved, <{sigma_limit} marginal)')
    ax.text(0.02, 0.97, txt, transform=ax.transAxes, va='top', ha='left',
            fontsize=9, zorder=11,
            bbox=dict(boxstyle='round', fc='white', ec=color, lw=1.5, alpha=0.85))
    # Draw the legend last and force it above all data bars, fit curves, and the
    # annotation box so it is never occluded.
    leg = ax.legend(loc='upper right', fontsize=8)
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
                        f'{additional_title}{suptitle} (Nimages = {nimages}): '
                        f'EXT {res["ext"]}, q = {res["q"]:.0f} $e^-$  '
                        f'(window {res["lo"]:.0f}-{res["hi"]:.0f}, '
                        f'{res["n_components"]}/{res["expected_peaks"]} peaks)',
                        fontsize=11, pad=10)
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
                fig.suptitle(f'{additional_title}{suptitle} at q = {q:g} $e^-$ '
                             f'(Nimages = {nimages})')
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
