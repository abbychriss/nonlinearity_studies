"""
Figures for the gain-independent analysis.

Every plot here is the same picture -- a per-extension scatter against peak
number, optionally with a fitted curve, a reference line, or both -- so the
individual / 2x2-together subplot (shared axes, suptitle, hidden inner
labels, residual panels under the fits, saving) lives once in
``_plot_scatter_ext`` and each public function is a thin wrapper that says what
to draw. The layout conventions match
``plotting.plot_all_peaks`` / ``plot_nonlinearity``.

Styling is otherwise left at matplotlib's defaults on purpose: every panel
carries a major grid under the data, colors are only chosen where a line has to
be told apart from the scatter, and the marker size/opacity are set so a
thousand-peak scatter stays readable.
"""

import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from pedestal_subtract.core import _finish_fig
from pedestal_subtract.plotting import format_title

from ..plotting import _sci, _term
from .spectrum_fit import _polyval_asc


def _plot_scatter_ext(x_ext, y_ext, xlabel, ylabel, suptitle, base_suffix,
                      overlay_ext=None, hline_ext=None, vline_ext=None,
                      top_axis_ext=None,
                      resid_ext=None, resid_ylabel=None, resid_overlay_ext=None,
                      overlay_kwargs=None, resid_overlay_kwargs=None,
                      scatter_label=None, legend_text_ext=None, yscale='linear',
                      plot_individual=False, plot_together=False,
                      individual_figsize=(7, 6), subplots_figsize=(13, 9),
                      additional_title='', nimages=10, sharex=True, sharey=True,
                      show_titles=True, save_plots=False, show_plots=True,
                      fig_path='./', file='gain_indep_nonlinearity', dpi=350,
                      s=10, alpha=0.5):
    """Scatter y_ext[ext] vs x_ext[ext], one figure per extension and/or a 2x2 grid.

    overlay_ext / hline_ext / vline_ext are per-extension and optional:
        overlay_ext[ext] = (x, y, label)   a curve drawn over the scatter, with
                                           ``overlay_kwargs`` passed to plot
        hline_ext[ext]   = (value, label)  a horizontal reference line, or a
                                           list of them; an optional third
                                           element is a dict of line kwargs,
                                           where an 'x0' entry starts the line
                                           at that x instead of spanning the axis
        vline_ext[ext]   = (value, label)  a vertical reference line, same
                                           list / kwargs-dict forms as hline_ext
    Any of them may be None for an individual extension.

    ``top_axis_ext[ext] = (positions, labels, axis_label)`` adds a second x axis
    along the top of the scatter, with fixed ticks at ``positions`` (in the main
    axis' own data coordinates) reading ``labels`` -- for a companion quantity
    that is paired to x point by point rather than by a function of it. Build
    the ticks with ``_companion_axis_ticks``.

    ``scatter_label`` labels the scatter itself in the legend; ``legend_text_ext``
    is a per-extension string (or list of strings) that joins the legend with no
    line or marker of its own, for a number worth quoting next to the fits.

    Wherever a fitted curve is drawn, a short panel under the scatter carries
    the residuals about it, sharing the x axis. Extensions whose fit failed
    (overlay None) keep a single full-height panel, so a 2x2 grid can mix the
    two without the columns going ragged. ``resid_ext`` puts a residual of the
    caller's own choosing in that panel instead -- one array per extension,
    aligned with x_ext, labelled ``resid_ylabel`` -- so a plot can carry a
    residual that is not taken about the overlay, or with no overlay at all.
    ``resid_overlay_ext[ext] = (x, y, label)`` draws a fitted curve in that
    panel instead of over the scatter, in red unless ``resid_overlay_kwargs``
    says otherwise, and with a small legend so it fits in the short panel.
    """
    if not (plot_individual or plot_together):
        return

    fig_name = Path(fig_path) / f'{Path(file).stem}_{base_suffix}'

    if resid_ext is None and overlay_ext is not None:
        resid_ext = [_overlay_residual(x, y, o)
                     for x, y, o in zip(x_ext, y_ext, overlay_ext)]
    if resid_ext is not None and all(r is None for r in resid_ext):
        resid_ext = None
    resid_ylabel = resid_ylabel or _residual_ylabel(ylabel)

    def _has_resid(ext):
        return resid_ext is not None and resid_ext[ext] is not None

    def _draw(ax, ext, ax_resid=None):
        # Grid first, and pinned below the data, so every scatter, fit and
        # reference line in both panels renders on top of it.
        for panel in (ax, ax_resid):
            if panel is not None:
                panel.grid(True, which='major')
                panel.set_axisbelow(True)
        ax.scatter(x_ext[ext], y_ext[ext], s=s, alpha=alpha, label=scatter_label)
        labelled = scatter_label is not None
        if overlay_ext is not None and overlay_ext[ext] is not None:
            ox, oy, olabel = overlay_ext[ext]
            ax.plot(ox, oy, label=olabel, **(overlay_kwargs or {}))
            labelled = labelled or olabel is not None
        if hline_ext is not None and hline_ext[ext] is not None:
            for entry in _as_line_list(hline_ext[ext]):
                value, hlabel = entry[0], entry[1]
                hkwargs = dict(entry[2]) if len(entry) > 2 else {}
                # An 'x0' in the kwargs is a start in data coordinates, so the
                # line runs from there to the last point of the scatter instead
                # of across the whole axis and out into the left margin.
                x0 = hkwargs.pop('x0', None)
                if x0 is None or not len(x_ext[ext]):
                    ax.axhline(value, label=hlabel, **hkwargs)
                else:
                    ax.plot([x0, np.max(x_ext[ext])], [value] * 2,
                            label=hlabel, **hkwargs)
                labelled = labelled or hlabel is not None
        if vline_ext is not None and vline_ext[ext] is not None:
            for entry in _as_line_list(vline_ext[ext]):
                value, vlabel = entry[0], entry[1]
                ax.axvline(value, label=vlabel, **(entry[2] if len(entry) > 2 else {}))
                labelled = labelled or vlabel is not None
        if legend_text_ext is not None and legend_text_ext[ext] is not None:
            texts = legend_text_ext[ext]
            for text in ([texts] if isinstance(texts, str) else texts):
                ax.plot([], [], ' ', label=text)
                labelled = True
        if labelled:
            ax.legend()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_yscale(yscale)

        if top_axis_ext is not None and top_axis_ext[ext] is not None:
            positions, tick_labels, top_label = top_axis_ext[ext]
            # A twinned axis rather than ax.secondary_xaxis: the companion values
            # are paired to x sample by sample, not by an invertible function of
            # it, so the ticks are pinned by hand and only the limits have to
            # follow the scatter (which a shared x axis can still change later).
            ax_top = ax.twiny()
            ax_top.set_xlim(ax.get_xlim())
            ax.callbacks.connect(
                'xlim_changed', lambda a, top=ax_top: top.set_xlim(a.get_xlim()))
            ax_top.set_xticks(positions, labels=tick_labels)
            ax_top.set_xlabel(top_label)

        if ax_resid is not None:
            ax_resid.scatter(x_ext[ext], resid_ext[ext], s=s, alpha=alpha)
            ax_resid.axhline(0.0, color='k', ls='--', alpha=0.8)
            if resid_overlay_ext is not None and resid_overlay_ext[ext] is not None:
                rx, ry, rlabel = resid_overlay_ext[ext]
                ax_resid.plot(rx, ry, label=rlabel,
                              **(resid_overlay_kwargs or dict(color='red')))
                if rlabel is not None:
                    ax_resid.legend()
            ax_resid.set_ylabel(resid_ylabel)
            # The residual panel is the bottom of the pair, so it takes over the
            # x label and the tick labels from the scatter above it.
            ax_resid.set_xlabel(xlabel)
            ax.set_xlabel('')
            ax.tick_params(labelbottom=False)

    def _bottom(main, resid):
        """The axis of a panel that carries the x label -- residual if there is one."""
        return main if resid is None else resid

    if plot_individual:
        for ext in range(len(y_ext)):
            if _has_resid(ext):
                fig, (ax, ax_resid) = plt.subplots(
                    2, 1, figsize=individual_figsize, constrained_layout=True,
                    sharex=True, gridspec_kw={'height_ratios': [3, 1]})
            else:
                fig, ax = plt.subplots(1, 1, figsize=individual_figsize,
                                       constrained_layout=True)
                ax_resid = None
            if show_titles:
                ax.set_title(format_title(additional_title, suptitle, nimages, f'EXT {ext + 1}'),
                             pad=10)
            _draw(ax, ext, ax_resid)
            if save_plots:
                output_path = fig_name.with_stem(fig_name.stem + f'_EXT{ext+1}').with_suffix('.jpeg')
                plt.savefig(str(output_path), dpi=dpi)
                print(f'Saved plot to {output_path}')
            _finish_fig(show_plots)

    if plot_together:
        # One gridspec row per panel-half rather than plt.subplots(2, 2), so a
        # scatter and its residuals can sit in the same grid cell. A panel with
        # no residuals spans both of its rows and looks exactly as it did before.
        fig = plt.figure(figsize=subplots_figsize, constrained_layout=True)
        gs = fig.add_gridspec(4, 2, height_ratios=[3, 1, 3, 1])
        if show_titles:
            fig.suptitle(format_title(additional_title, suptitle, nimages))

        mains, resids = [], []
        for ext in range(len(y_ext)):
            row, col = 2 * (ext // 2), ext % 2
            share = dict(sharex=mains[0] if (mains and sharex) else None,
                         sharey=mains[0] if (mains and sharey) else None)
            if _has_resid(ext):
                ax = fig.add_subplot(gs[row, col], **share)
                first_resid = next((r for r in resids if r is not None), None)
                ax_resid = fig.add_subplot(gs[row + 1, col], sharex=ax,
                                           sharey=first_resid if sharey else None)
            else:
                ax = fig.add_subplot(gs[row:row + 2, col], **share)
                ax_resid = None
            mains.append(ax)
            resids.append(ax_resid)
            _draw(ax, ext, ax_resid)
            ax.set_title(f'EXT {ext + 1}')

        # Same shared-axis convention as the other 2x2 "together" plots: hide the
        # x-label on the top row and the y-label on the right column, since the
        # shared axis already carries that information.
        for i in (0, 1):
            bottom = _bottom(mains[i], resids[i])
            bottom.set_xlabel('')
            bottom.tick_params(labelbottom=True)
        for i in (1, 3):
            mains[i].set_ylabel('')
            mains[i].set_yscale(yscale)
            mains[i].tick_params(labelleft=True)
            if resids[i] is not None:
                resids[i].set_ylabel('')
                resids[i].tick_params(labelleft=True)

        if save_plots:
            output_path = fig_name.with_suffix('.jpeg')
            plt.savefig(str(output_path), dpi=dpi)
            print(f'Saved plot to {output_path}')
        _finish_fig(show_plots)


def _companion_axis_ticks(x, companion, n_ticks=8, fmt='%g'):
    """Ticks along x labelled with the companion value paired to each sample.

    The companion (the peak number i behind each candidate gain) is paired to x
    point by point but is not a function of it: g_fit(i) wiggles, so nearly the
    same gain can come back at a later i and there is no exact gain -> i inverse
    to hand ``secondary_xaxis``. Ticks are therefore placed at evenly spaced x
    targets and each is labelled with the companion of the nearest sample, which
    is unambiguous as long as the wiggle stays small next to the tick spacing.
    Returns (positions, labels), or None if x has nothing finite to place them on.
    """
    x = np.asarray(x, dtype=float)
    companion = np.asarray(companion)
    good = np.isfinite(x)
    if good.sum() < 2:
        return None
    x, companion = x[good], companion[good]
    targets = np.linspace(x.min(), x.max(), n_ticks)
    # unique() also drops the duplicates left where several targets fall in the
    # same gap between sampled gains, so those ticks are simply not drawn.
    idx = np.unique(np.abs(x[None, :] - targets[:, None]).argmin(axis=1))
    return x[idx], [fmt % c for c in companion[idx]]


def _as_line_list(entry):
    """One (value, label[, kwargs]) reference line, or a list of them, as a list."""
    return entry if isinstance(entry[0], (tuple, list)) else [entry]


def _overlay_residual(x, y, overlay):
    """Scatter minus the curve drawn over it, sampled at the scatter's own x.

    The curve is normally evaluated on exactly the scatter's x grid, where this
    is a plain subtraction; np.interp covers the caller whose curve lives on a
    different grid (g_fit(i) drawn over the direct gain), and NaNs the points
    that fall outside it rather than extrapolating.
    """
    if overlay is None:
        return None
    ox, oy, _ = overlay
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ox = np.asarray(ox, dtype=float)
    oy = np.asarray(oy, dtype=float)
    if ox.shape == x.shape and np.array_equal(ox, x):
        return y - oy
    return y - np.interp(x, ox, oy, left=np.nan, right=np.nan)


def _g_fit_legend_text(fits_ext):
    r"""Per-extension 'g_fit = ...' legend entry: the fitted linear coefficient
    is not attached to any drawn line, so it joins the legend as text only."""
    return [r'$g_\mathrm{fit} = %s$ ADU' % _sci(fit['g_fit'])
            if np.isfinite(fit['g_fit']) else None
            for fit in fits_ext]


def _residual_ylabel(ylabel):
    """'residual' carrying the unit the main y label ends with, if it has one."""
    unit = re.search(r'\(([^()]*)\)\s*$', ylabel)
    return f'residual\n({unit.group(1)})' if unit else 'residual'


def plot_gain_indep_nonlinearity_ext(delta_g_ext, line_fits_ext=None, **kwargs):
    """delta_g(n) = q(n) - q(n-1) - g vs n, with the fitted line through it."""
    x_ext = [np.arange(1, len(d) + 1) for d in delta_g_ext]
    overlay_ext = None
    if line_fits_ext is not None:
        overlay_ext = [
            (x, fit['slope'] * x + fit['intercept'],
             r'$(%s)\,n %s$' % (_sci(fit['slope']), _term(fit['intercept'], sci=False)))
            if np.isfinite(fit['slope']) else None
            for x, fit in zip(x_ext, line_fits_ext)
        ]
    kwargs.setdefault('suptitle', 'Differential of the Gain')
    kwargs.setdefault('overlay_kwargs', dict(color='red'))
    _plot_scatter_ext(x_ext, delta_g_ext,
                      xlabel=r'peak number $n$',
                      ylabel=r'$q(n) - q(n-1) - g$ (ADU)',
                      base_suffix='gain_indep_nonlin',
                      overlay_ext=overlay_ext, **kwargs)


def plot_charge_spectrum_fit_ext(q_ext, fits_ext, cubic_fits_ext=None, **kwargs):
    """The charge spectrum's departure from a straight line, with the fitted
    nonlinear term on top.

    Plotting q(n) itself would just show a line of slope ~g with the curvature
    invisible, so the linear part of the fit is subtracted from both: the
    scatter is q(n) - (a0 + g_fit*n) and the curve is the quadratic fit's
    alpha*n^2. When the cubic fit is supplied, what it adds on top of that
    quadratic goes in the residual panel, where the residual about alpha*n^2 is
    the only place it is large enough to see.
    """
    x_ext = [np.arange(len(q)) for q in q_ext]
    y_ext = []
    overlay_ext = []
    resid_overlay_ext = []
    for ext, (q, fit) in enumerate(zip(q_ext, fits_ext)):
        n = x_ext[ext]
        linear = fit['a0'] + fit['g_fit'] * n
        y_ext.append(np.asarray(q, dtype=float) - linear)
        if not np.isfinite(fit['alpha']):
            overlay_ext.append(None)
            resid_overlay_ext.append(None)
            continue
        quad = fit['alpha'] * n ** 2
        overlay_ext.append(
            (n, quad, r'$\alpha n^2$, $\alpha = %s$ ADU' % _sci(fit['alpha'])))
        if cubic_fits_ext is not None and np.isfinite(cubic_fits_ext[ext]['beta']):
            cubic = cubic_fits_ext[ext]
            resid_overlay_ext.append(
                (n, _polyval_asc(cubic['coeffs'], n) - linear - quad,
                 r'cubic: $\alpha n^2 + \beta n^3$, $\alpha = %s$ ADU, $\beta = %s$ ADU'
                 % (_sci(cubic['alpha']), _sci(cubic['beta']))))
        else:
            resid_overlay_ext.append(None)
    if all(r is None for r in resid_overlay_ext):
        resid_overlay_ext = None

    kwargs.setdefault('suptitle', 'Charge Spectrum Nonlinear Term')
    kwargs.setdefault('overlay_kwargs', dict(color='red'))
    kwargs.setdefault('resid_overlay_kwargs', dict(color='purple'))
    kwargs.setdefault('legend_text_ext', _g_fit_legend_text(fits_ext))
    _plot_scatter_ext(x_ext, y_ext,
                      xlabel=r'peak number $n$',
                      ylabel=r'$q(n) - (a_0 + g_\mathrm{fit}\,n)$ (ADU)',
                      base_suffix='charge_spectrum_fit',
                      overlay_ext=overlay_ext,
                      resid_overlay_ext=resid_overlay_ext, **kwargs)


def plot_charge_spectrum_ext(q_ext, fits_ext, cubic_fits_ext=None, **kwargs):
    """The charge spectrum q(n) itself, with the parabola fit Q(n) over it.

    The companion of ``plot_charge_spectrum_fit_ext`` with nothing subtracted:
    the scatter is q(n), the curve is the whole quadratic fit
    Q(n) = a0 + g_fit*n + alpha*n^2, and the residual panel carries q(n) - Q(n),
    which is where the curvature the straight-looking spectrum hides becomes
    visible. The cubic fit, when supplied, is drawn in that panel as what it
    adds on top of the quadratic -- the same curve, and the same styling, as in
    the nonlinear-term plot.
    """
    x_ext = [np.arange(len(q)) for q in q_ext]
    y_ext = [np.asarray(q, dtype=float) for q in q_ext]
    overlay_ext = []
    resid_overlay_ext = []
    for ext, fit in enumerate(fits_ext):
        n = x_ext[ext]
        if not np.isfinite(fit['alpha']):
            overlay_ext.append(None)
            resid_overlay_ext.append(None)
            continue
        quad = _polyval_asc(fit['coeffs'], n)
        overlay_ext.append(
            (n, quad,
             "\n".join((r'$Q(n) = a_0 + g_\mathrm{fit}\,n + \alpha n^2$',
                        r'$a_0 = %s$ ADU, $\alpha = %s$ ADU'))
             % (_sci(fit['a0']), _sci(fit['alpha']))))
        if cubic_fits_ext is not None and np.isfinite(cubic_fits_ext[ext]['beta']):
            cubic = cubic_fits_ext[ext]
            resid_overlay_ext.append(
                (n, _polyval_asc(cubic['coeffs'], n) - quad,
                 r'cubic: $\alpha n^2 + \beta n^3$, $\alpha = %s$ ADU, $\beta = %s$ ADU'
                 % (_sci(cubic['alpha']), _sci(cubic['beta']))))
        else:
            resid_overlay_ext.append(None)
    if all(r is None for r in resid_overlay_ext):
        resid_overlay_ext = None

    kwargs.setdefault('suptitle', 'Charge Spectrum')
    kwargs.setdefault('overlay_kwargs', dict(color='red'))
    kwargs.setdefault('resid_overlay_kwargs', dict(color='purple'))
    kwargs.setdefault('legend_text_ext', _g_fit_legend_text(fits_ext))
    _plot_scatter_ext(x_ext, y_ext,
                      xlabel=r'peak number $n$',
                      ylabel=r'$q(n)$ (ADU)',
                      base_suffix='charge_spectrum',
                      overlay_ext=overlay_ext,
                      resid_overlay_ext=resid_overlay_ext, **kwargs)


def plot_cumulative_delta_ext(S_ext, alpha_fits_ext=None, **kwargs):
    """S(n) = (1/n) * sum delta_g vs n, with the fitted line whose slope is alpha'."""
    x_ext = [np.arange(1, len(S) + 1) for S in S_ext]
    overlay_ext = None
    if alpha_fits_ext is not None:
        overlay_ext = [
            (x, fit['slope'] * x + fit['intercept'],
             "\n".join((r"$S(n) = \alpha' n + b$",
                         r"$\alpha' = %s$ ADU, $b = %s$ ADU"))
             % (_sci(fit['slope']), _sci(fit['intercept'])))
            if np.isfinite(fit['slope']) else None
            for x, fit in zip(x_ext, alpha_fits_ext)
        ]
    kwargs.setdefault('suptitle', 'Cumulative Gain Differential')
    kwargs.setdefault('sharey', False)
    kwargs.setdefault('overlay_kwargs', dict(color='red'))
    # The direct-gain S(n) reuses this function, so the caller can override the
    # suffix to keep the two sets of saved figures apart.
    kwargs.setdefault('base_suffix', 'cumulative_delta')
    _plot_scatter_ext(x_ext, S_ext,
                      xlabel=r'peak number $n$',
                      ylabel=r'$S(n) = I(n)/n$ (ADU)',
                      overlay_ext=overlay_ext, **kwargs)


def plot_running_gain_ext(running_ext, g_ext, window_fits_ext=None,
                          residual_ext=None, **kwargs):
    """g_fit(i), the average gain over the first i peaks, vs i, against the standard gain.

    ``window_fits_ext`` (from ``fit_running_gain_window_ext``) adds, over the
    window of fit ranges it was built with, the average of g_fit(i) as a second
    horizontal line and the least-squares line through it in the residual panel.

    ``residual_ext`` is g_fit(i) - g, which rides along in the residual panel
    under the scatter rather than in a figure of its own. The window fit's
    trend line goes in that same panel, shifted by -g so it lies over the
    residual it was fitted to.
    """
    kwargs.setdefault('suptitle', 'Fitted Gain from ADU Charge Spectrum')
    # A dense thousand-point running curve, so smaller and fainter markers than
    # the per-peak plots use.
    kwargs.setdefault('s', 4)
    kwargs.setdefault('alpha', 0.25)
    # g_fit(i) settles at a different level in each extension, so a shared y
    # axis would flatten all but the widest-swinging one.
    kwargs.setdefault('sharey', False)
    kwargs.setdefault('resid_overlay_kwargs', dict(color='purple'))
    hline_ext = [[(g, r'$g(1) = %.3f$ ADU/$e^-$' % g,
                   dict(color='#e38527', x0=0))] for g in g_ext]
    resid_overlay_ext = None
    if window_fits_ext is not None:
        resid_overlay_ext = []
        for win, hlines, g in zip(window_fits_ext, hline_ext, g_ext):
            if win['i'].size:
                hlines.append((win['mean'],
                               r'$\langle g_\mathrm{fit}\rangle_{%d}^{%d} = %.3f$ ADU/$e^-$'
                               % (win['i'][0], win['i'][-1], win['mean']),
                               dict(color='red')))
            resid_overlay_ext.append(
                (win['i'], win['slope'] * win['i'] + win['intercept'] - g,
                 r'$y = (%s)\,n %s$ (ADU/$e^-$)'
                 % (_sci(win['slope']), _term(win['intercept'] - g, sci=False)))
                if np.isfinite(win['slope']) else None)
    _plot_scatter_ext([run['i'] for run in running_ext],
                      [run['g_fit'] for run in running_ext],
                      xlabel=r'peaks included in fit $i$',
                      ylabel=r'$g_\mathrm{fit}(i)$ (ADU/$e^-$)',
                      base_suffix='running_gain',
                      hline_ext=hline_ext,
                      resid_ext=residual_ext,
                      resid_overlay_ext=resid_overlay_ext,
                      resid_ylabel='$g_\\mathrm{fit}(i) - g(1)$\n(ADU/$e^-$)', **kwargs)


def plot_gain_fit_residual_ext(running_ext, residual_ext, **kwargs):
    """g_fit(i) - g vs i: the fitted-spectrum counterpart of delta_g."""
    kwargs.setdefault('suptitle', 'Average Gain Residual')
    _plot_scatter_ext([run['i'] for run in running_ext], residual_ext,
                      xlabel=r'peaks included in fit $i$',
                      ylabel=r'$g_\mathrm{fit}(i) - g$ (ADU/$e^-$)',
                      base_suffix='gain_fit_residual',
                      hline_ext=[(0.0, None)] * len(residual_ext), **kwargs)


def plot_gain_fit_residual_cumsum_ext(running_ext, cumsum_ext, parabola_fits_ext=None,
                                      **kwargs):
    """Cumulative sum of g_fit(i) - g over i: another nonlinearity estimate, in ADU.

    ``parabola_fits_ext`` (from ``fit_gain_fit_residual_cumsum_ext``) is drawn
    over the scatter with its equation in the legend.
    """
    kwargs.setdefault('suptitle', 'Integrated Average Gain Residual')
    # Same dense running curve as the average-gain plot, so the same small markers.
    kwargs.setdefault('s', 4)
    # The four extensions span very different integrated residuals, so this one
    # plot keeps its own y scale rather than the shared one the others use.
    kwargs['sharey'] = False
    x_ext = [run['i'] for run in running_ext]
    overlay_ext = None
    if parabola_fits_ext is not None:
        overlay_ext = [
            (x, _polyval_asc(fit['coeffs'], x),
             r'$y = (%s)\,i^2 (%s)\,i %s$ (ADU)'
             % (_sci(fit['alpha']), _term(fit['g_fit']), _term(fit['a0'], sci=False)))
            if np.isfinite(fit['alpha']) else None
            for x, fit in zip(x_ext, parabola_fits_ext)
        ]
    _plot_scatter_ext(x_ext, cumsum_ext,
                      xlabel=r'peaks included in fit $i$',
                      ylabel=r'$\sum_i (g_\mathrm{fit}(i) - g)$ (ADU)',
                      base_suffix='gain_fit_residual_cumsum',
                      overlay_ext=overlay_ext, overlay_kwargs=dict(color='red'),
                      **kwargs)


def plot_direct_gain_ext(g_dir_ext, g_ext, running_ext=None, **kwargs):
    """g_dir(n) = (q(n) - q(0))/n vs n, against g and (optionally) g_fit(i)."""
    x_ext = [np.arange(1, len(g_dir) + 1) for g_dir in g_dir_ext]
    overlay_ext = None
    if running_ext is not None:
        overlay_ext = [(run['i'], run['g_fit'], r'$g_\mathrm{fit}(n)$') for run in running_ext]
    kwargs.setdefault('suptitle', 'Direct Gain per Peak')
    # g_fit(i) is a measured value at every i just like the scatter it rides
    # over, not a fitted curve, so it is drawn as points too -- same marker at
    # the same size, since scatter's s is an area and plot's ms a diameter.
    kwargs.setdefault('s', 8)
    kwargs.setdefault('alpha', 0.5)
    kwargs.setdefault('overlay_kwargs',
                      dict(ls='none', marker='o', ms=np.sqrt(kwargs['s']),
                           color='purple', alpha=kwargs['alpha']))
    _plot_scatter_ext(x_ext, g_dir_ext,
                      xlabel=r'peak number $n$',
                      ylabel=r'$g_\mathrm{dir}(n)$ (ADU/$e^-$)',
                      base_suffix='direct_gain',
                      overlay_ext=overlay_ext,
                      scatter_label=r'$g_\mathrm{dir}(n)$',
                      # The residual is taken about the g_fit(i) overlay, not
                      # about g(1); narrow panel, so the unit gets its own line.
                      resid_ylabel='$g_\\mathrm{dir}(n) - g_\\mathrm{fit}(n)$\n(ADU/$e^-$)',
                      hline_ext=[(g, r'$g(1) = %.3f$ ADU/$e^-$' % g,
                                  dict(color='#e38527', x0=0))
                                 for g in g_ext], **kwargs)


def plot_direct_delta_ext(delta_dir_ext, line_fits_ext=None, **kwargs):
    """delta_dir(n) = g - g_dir(n) vs n, with the fitted line through it."""
    x_ext = [np.arange(1, len(d) + 1) for d in delta_dir_ext]
    overlay_ext = None
    if line_fits_ext is not None:
        overlay_ext = [
            (x, fit['slope'] * x + fit['intercept'],
             r'$y = (%s)\,n %s$ (ADU/$e^-$)'
             % (_sci(fit['slope']), _term(fit['intercept'], sci=False)))
            if np.isfinite(fit['slope']) else None
            for x, fit in zip(x_ext, line_fits_ext)
        ]
    kwargs.setdefault('suptitle', 'Direct Gain Differential')
    kwargs.setdefault('overlay_kwargs', dict(color='red'))
    _plot_scatter_ext(x_ext, delta_dir_ext,
                      xlabel=r'peak number $n$',
                      ylabel=r'$g - g_\mathrm{dir}(n)$ (ADU/$e^-$)',
                      base_suffix='direct_delta',
                      overlay_ext=overlay_ext, **kwargs)


def plot_gain_scan_ext(scan_ext, eval_charge=500.0, **kwargs):
    """Standard nonlinearity at eval_charge vs the gain used to calibrate, with
    the minimum-|nonlinearity| gain marked."""
    kwargs.setdefault('suptitle', f'Nonlinearity at {eval_charge:g} $e^-$ vs Gain')
    kwargs.setdefault('sharex', False)
    vline_ext = [
        (scan['best_gain'],
         r'$g_\mathrm{fit}(i=%d) = %.3f$ ADU/$e^-$, nonlin($%g\,e^-$) $= %+.3f\,e^-$'
         % (scan['best_i'], scan['best_gain'], eval_charge, scan['best_nonlin']),
         dict(color='red', alpha=0.5, ls='--'))
        if np.isfinite(scan['best_gain']) else None
        for scan in scan_ext
    ]
    # The direct gain at the same charge is not one of the scanned candidates, so
    # the nonlinearity it leaves behind is quoted in the legend rather than drawn.
    legend_text_ext = [
        r'$g_\mathrm{dir}(%g\,e^-) = %.3f$ ADU/$e^-$, nonlin($%g\,e^-$) $= %+.3f\,e^-$'
        % (eval_charge, scan['g_dir_eval'], eval_charge, scan['nonlin_dir'])
        if np.isfinite(scan.get('nonlin_dir', np.nan)) else None
        for scan in scan_ext
    ]
    # Each scanned gain is the g_fit(i) of one peak number, so the top axis
    # carries that i alongside the gain the panel is actually plotted against.
    top_axis_ext = [
        (ticks[0], ticks[1], r'peaks included in fit $i$') if ticks else None
        for ticks in (_companion_axis_ticks(scan['gain'], scan['i'])
                      for scan in scan_ext)
    ]
    _plot_scatter_ext([scan['gain'] for scan in scan_ext],
                      [scan['nonlin'] for scan in scan_ext],
                      xlabel=r'calibration gain $g_\mathrm{fit}(i)$ (ADU/$e^-$)',
                      ylabel=r'nonlinearity at %g $e^-$ ($e^-$)' % eval_charge,
                      base_suffix='gain_scan', top_axis_ext=top_axis_ext,
                      vline_ext=vline_ext, legend_text_ext=legend_text_ext, **kwargs)

# New plotting code written by Abby!! They will probably error because I suck at writing code :} <-- curly braces = smiley face with tοngue sticking out. I hope someone finds this so they know I'm still in here.
def plot_alpha_alpha_prime_vs_gain_and_index_ext(index_alpha_gfit_ext, index_alpha_prime_ext, individual_figsize=(9,6.5), together_figsize=(16,9), plot_together=True, plot_individual=False):
    """Plots alpha and/or alpha' against the index (peak number) used as the maximum for each fit. Both are optional but raises value error if neither are given (why would you do that though)."""

    if index_alpha_gfit_ext is None and index_alpha_prime_ext is None:
        raise ValueError(f"Neither alpha nor alpha' were supplied (both were None), but at least one is needed.")

    n_ext = len(index_alpha_gfit_ext if index_alpha_gfit_ext is not None
                else index_alpha_prime_ext)

    def _draw(ax, ext):
        index_alpha = alpha = gain = None
        if index_alpha_gfit_ext is not None:
            index_alpha = index_alpha_gfit_ext[ext]["i"]
            alpha = index_alpha_gfit_ext[ext]["alpha"]
            gain = index_alpha_gfit_ext[ext]["g_fit"]

        index_alpha_prime = alpha_prime = None
        if index_alpha_prime_ext is not None:
            index_alpha_prime = index_alpha_prime_ext[ext]["i"]
            alpha_prime = index_alpha_prime_ext[ext]["alpha_prime"]

        if alpha is not None:
            ax.scatter(index_alpha, alpha, c='blue', alpha=0.7, s=6, label=r'$\alpha$')

        if alpha_prime is not None:
            ax.scatter(index_alpha_prime, alpha_prime, c='pink', alpha=0.7, s=6,
                       label=r"$\alpha'$")

        ax.set_xlabel(r'Peak Number $n$')
        ax.set_ylabel('Nonlinearity parameter (ADU)')

        handles, labels = ax.get_legend_handles_labels()

        # g_fit(i) shares the peak-number axis but not the units, so it goes on
        # a twinned y axis rather than a secondary x axis: g_fit(i) is not
        # monotonic in i, so there is no well-defined gain -> index inverse for
        # secondary_xaxis to call.
        if gain is not None:
            ax_gain = ax.twinx()
            ax_gain.scatter(index_alpha, gain, c='green', alpha=0.7, s=6,
                            label=r'$g_\mathrm{fit}(i)$')
            ax_gain.set_ylabel('Gain (ADU/e-)')
            gain_handles, gain_labels = ax_gain.get_legend_handles_labels()
            handles += gain_handles
            labels += gain_labels

        if handles:
            ax.legend(handles, labels)

    if plot_together:
        fig, axs = plt.subplots(2, 2, figsize=together_figsize, constrained_layout=True)
        axs = axs.flatten()
        fig.suptitle('Nonlinearity Parameter vs Peak Number / Gain')

        for ext in range(n_ext):
            _draw(axs[ext], ext)
        plt.show()

    if plot_individual:
        for ext in range(n_ext):
            fig, ax = plt.subplots(1, 1, figsize=individual_figsize, constrained_layout=True)
            fig.suptitle('Nonlinearity Parameter vs Peak Number / Gain')
            _draw(ax, ext)
            plt.show()

def plot_nonlin_vs_gain(nonlinearity_at_gains_ext, gain_range, q, **kwargs):
    """Plots nonlinearity at (the standard measured charge - peak#) against gain."""

    kwargs.setdefault('suptitle', f'Nonlinearity at {q:g}e- vs Gain')
    kwargs.setdefault('s', 6)
    overlay_ext = None
    # Every extension was scanned over the same gains, but _plot_scatter_ext
    # indexes x per extension just like y, so the single gain axis is repeated
    # rather than handed over as one bare array. q is formatted '%g' so a float
    # charge does not put a '.' in the file stem, where the with_suffix('.jpeg')
    # that saves the figure would read it as an extension and drop the rest.
    _plot_scatter_ext([np.asarray(gain_range, dtype=float)] * len(nonlinearity_at_gains_ext),
                          nonlinearity_at_gains_ext,
                          xlabel='Gain (ADU/e-)',
                          ylabel=f'Nonlinearity at {q:g}e-',
                          base_suffix=f'nonlin_at_{q:g}_vs_gain',
                          overlay_ext=overlay_ext,
                          **kwargs)