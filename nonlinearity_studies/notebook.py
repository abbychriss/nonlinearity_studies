"""
notebook — helpers for the example notebook (scripts/nonlinearity_analysis.ipynb).

The pipelines draw and close their figures themselves, so the notebook runs them as they
are and collects what they draw (capture_figures) instead of reimplementing any step.
The figures are then shown behind one tab per extension (show_panel).
"""
import contextlib
import copy
import glob
import io
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pedestal_subtract.calibrate import full_data_stats
from pedestal_subtract.core import convert_to_electrons

from .nonlinearity_studies import get_zero_one_peaks_ext, plot_all_peaks, plot_zero_one_peaks
from .plotting import _xlim_window

CONFIG_DIR = Path(__file__).resolve().parents[1] / 'config'

# The standard pipeline's plots, as named by their --plot_<name>_individual flags and the
# functions that draw them (looked up on run_nonlinearity_studies).
# (The all-peaks and zero/one peak plots are left out: all_peaks_explorer and zero_one_explorer
# draw them, for whatever window / binning is asked for.)
STANDARD_PLOTS = ('charge_per_column', 'nonlinearity')
STANDARD_PLOT_FUNCTIONS = ('plot_charge_per_column', 'plot_nonlinearity')


def _set_plot_switches(args, plots):
    """Turn on the individual (one figure per extension) plots named in ``plots`` (every one if None)
    and turn everything else, including all the combined 2x2 plots, off."""
    for key in [k for k in vars(args) if k.startswith('plot_') and k.endswith(('_individual', '_together'))]:
        name, kind = key[len('plot_'):].rsplit('_', 1)
        setattr(args, key, kind == 'individual' and (plots is None or name in plots))


def load_args(init_argparse, base_config_file, images, stitch_fits, plots=None):
    """Parse ``base_config_file`` with ``init_argparse``, taking every key from it except:

    - the image to analyze (``images``, a path or glob; a single file unless ``stitch_fits``)
      and ``stitch_fits`` itself, which come from the notebook;
    - the plot switches: only the per-extension figures named in ``plots`` (all of the
      parser's if None) are drawn, since the notebook shows one extension at a time;
    - ``show_plots`` / ``save_output``: the notebook shows the figures itself and writes
      nothing beyond the pipeline's own caches (no plots, summary CSV or meta-study entry);
    - ``verbose``, which is on, and ``extra_plot_title``, which is blanked so no title is added to the plots.

    ``base_config_file`` is a path, or the name of a file in the repo's config/ folder.
    """
    config = Path(base_config_file)
    config = config if config.exists() else CONFIG_DIR / config
    matches = sorted(glob.glob(images))
    if not matches:
        raise FileNotFoundError(f'No files match {images}')
    if not stitch_fits and len(matches) > 1:
        raise ValueError(f'stitch_fits=False expects a single image, but {images} matches {len(matches)} files:\n  '
                         + '\n  '.join(matches))

    args = init_argparse(['-j', str(config), '--stitch_fits' if stitch_fits else '--no-stitch_fits',
                          images if stitch_fits else matches[0]])
    args.show_plots = args.save_output = False
    args.verbose, args.extra_plot_title = True, ''
    _set_plot_switches(args, plots)
    return args


@contextlib.contextmanager
def capture_figures(*targets, dpi=100):
    """Collect, as PNG bytes, the figures drawn inside the block instead of showing them.

    ``targets`` are (owner, name) pairs: a module and the name of a plot function looked up
    on it while the block runs. Each target's figures are filed under ``figures[name]`` in the
    order they were drawn; figures drawn by anything else land in ``figures['figures']``.
    The plot functions end by closing their figure (or, for a few, by plt.show()), which is
    where each one is grabbed.
    """
    figures, current = defaultdict(list), ['figures']
    close, show = plt.close, plt.show
    originals = [(owner, name, getattr(owner, name)) for owner, name in targets]

    def grab(fig):
        buffer = io.BytesIO()
        fig.savefig(buffer, format='png', dpi=dpi)
        figures[current[-1]].append(buffer.getvalue())

    def close_and_grab(*args, **kwargs):
        if not args and plt.get_fignums():
            grab(plt.gcf())
        close(*args, **kwargs)

    def tag(name, function):
        def wrapper(*args, **kwargs):
            before = set(plt.get_fignums())
            current.append(name)
            try:
                return function(*args, **kwargs)
            finally:
                for num in [n for n in plt.get_fignums() if n not in before]:  # drawn but never closed
                    grab(plt.figure(num))
                    close(num)
                current.pop()
        return wrapper

    plt.close, plt.show = close_and_grab, lambda *args, **kwargs: None
    for owner, name, function in originals:
        setattr(owner, name, tag(name, function))
    try:
        yield figures
    finally:
        plt.close, plt.show = close, show
        for owner, name, function in originals:
            setattr(owner, name, function)


def figures_of(plot_function, *args, **kwargs):
    """The figures (PNG bytes) drawn by ``plot_function(*args, **kwargs)``."""
    with capture_figures() as figures:
        plot_function(*args, **{'show_plots': False, **kwargs})
    return figures['figures']


def show_panel(pngs, label=None, n_ext=4):
    """Display figures (PNG bytes, in extension order) behind one tab per extension.

    More than ``n_ext`` figures are split into consecutive panels of ``n_ext``, for plot
    functions that draw several sets of per-extension figures in one call.
    """
    import ipywidgets as widgets
    from IPython.display import display

    panels = [pngs[i:i + n_ext] for i in range(0, len(pngs), n_ext)]
    for k, panel in enumerate(panels, 1):
        tabs = widgets.Tab([widgets.Image(value=png, format='png') for png in panel])
        tabs.titles = [f'EXT {i + 1}' for i in range(len(panel))]
        if label is not None:
            title = label if len(panels) == 1 else f'{label} ({k}/{len(panels)})'
            tabs = widgets.VBox([widgets.HTML(f'<b>{title}</b>'), tabs])
        display(tabs)


def _explorer(controls, draw, connect):
    """Show ``controls`` above an output that ``draw()`` fills, drawn once now and again whenever the
    controls fire; ``connect(redraw)`` attaches the redraw callback to them."""
    import ipywidgets as widgets
    from IPython.display import display

    output = widgets.Output()

    def redraw(*_):
        with output:
            output.clear_output(wait=True)
            draw()

    connect(redraw)
    display(widgets.HBox(controls), output)
    redraw()


def _scale_toggle(value):
    """Linear / log y-scale buttons, starting at ``value`` (the run's configured y-scale)."""
    import ipywidgets as widgets
    return widgets.ToggleButtons(options=list(dict.fromkeys(['linear', 'log', value])), value=value, description='y scale')


def all_peaks_explorer(results):
    """The all-peaks plot (one tab per extension) with boxes for the x-limits, in electrons, and
    a linear / log y-scale toggle.

    Drawn first for the run's own ``plot_all_peaks_xlim`` and ``plot_all_peaks_yscale``; "Plot"
    redraws it for the window typed in (the peak lines and numbers, and the histogram, follow the
    window), and the toggle redraws it at once.
    """
    import ipywidgets as widgets

    kwargs = results['all_peaks_plot_kwargs']
    left, right = _xlim_window(kwargs.get('xlim')) or (1000, 1010)
    inputs = [widgets.FloatText(value=value, description=name, layout=widgets.Layout(width='190px'))
              for name, value in (('xlim left', left), ('xlim right', right))]
    button, scale = widgets.Button(description='Plot'), _scale_toggle(kwargs['yscale'])

    def draw():
        window = tuple(sorted(box.value for box in inputs))
        if window[0] == window[1]:
            print('The x-limits must differ.')
            return
        show_panel(figures_of(plot_all_peaks, results['counts_ext'], results['peaks_ext'], results['centers_ext'],
                              **{**kwargs, 'xlim': window, 'yscale': scale.value,
                                 'plot_individual': True, 'plot_together': False}))

    _explorer([*inputs, button, scale], draw, lambda redraw: (button.on_click(redraw), scale.observe(redraw, names='value')))


def zero_one_explorer(results, unit):
    """The zero/one peak fit plots in ``unit`` ('ADU' or 'electrons'), one tab per extension, with a
    slider for the number of bins of the histogram being fit (the ``zero_one_n_bins`` option) and
    a linear / log y-scale toggle (starting at ``plot_zero_one_yscale``).

    Each slider position repeats the double-Gaussian fit at that binning and redraws, starting from
    the run's own ``zero_one_n_bins``. It is for looking only: the pedestals and gains the rest of
    the notebook uses stay those of the run. Use one explorer per unit for separate sliders.
    """
    import ipywidgets as widgets

    if unit not in ('ADU', 'electrons'):
        raise ValueError(f"unit must be 'ADU' or 'electrons', got {unit!r}")
    fit_kwargs, data_ext = results['zero_one_fit_kwargs'], results['data_ext']
    fit_kwargs = {**fit_kwargs, 'data_stats': [full_data_stats(data) for data in data_ext]}  # same at every binning
    slider = widgets.IntSlider(value=fit_kwargs['n_bins'], min=10, max=max(1000, fit_kwargs['n_bins']),
                               description=f'{unit} bins', continuous_update=False,
                               style={'description_width': 'initial'}, layout=widgets.Layout(width='650px'))
    scale = _scale_toggle(results['zero_one_plot_kwargs']['yscale'])
    plot_kwargs = {**results['zero_one_plot_kwargs'], 'plot_adu_individual': unit == 'ADU',
                   'plot_electrons_individual': unit == 'electrons',
                   'plot_adu_together': False, 'plot_electrons_together': False}

    def draw():
        fit = get_zero_one_peaks_ext(data_ext, **{**fit_kwargs, 'n_bins': slider.value})
        show_panel(figures_of(plot_zero_one_peaks, data_ext, *fit, **{**plot_kwargs, 'yscale': scale.value}))

    _explorer([slider, scale], draw, lambda redraw: [w.observe(redraw, names='value') for w in (slider, scale)])


def run_standard(base_config_file, images, stitch_fits):
    """Run the nonlinearity pipeline; returns (args, results, figures).

    ``results`` is run_pipeline's dict and ``figures`` holds the PNGs of each plot function
    (keyed by function name), one figure per extension. The zero/one peaks come as the ADU
    figures followed by the electron figures.
    """
    from . import run_nonlinearity_studies as run
    args = load_args(run.init_argparse, base_config_file, images, stitch_fits, plots=STANDARD_PLOTS)
    with capture_figures(*[(run, name) for name in STANDARD_PLOT_FUNCTIONS]) as figures:
        results = run.run_pipeline(args)
    return args, results, figures


def run_gain_independent(base_config_file, image):
    """Run the gain-independent pipeline on one (already stitched) image, drawing every one of its
    plots; returns (args, results, figures) with the PNGs keyed by plot function name."""
    from .gain_independent_nonlinearity import plots, plot_nonlinearity_curves, run
    args = load_args(run.init_argparse, base_config_file, str(image), stitch_fits=False)
    plot_functions = [(plots, name) for name, function in vars(plots).items()
                      if name.startswith('plot_') and getattr(function, '__module__', None) == plots.__name__]
    with plt.rc_context(), capture_figures((run, 'plot_all_peaks'), (plot_nonlinearity_curves, 'plot_nonlinearity_curves'),
                                           *plot_functions) as figures:
        results = run.main(args)
    return args, results, figures


def correct_nonlinearity(args, results):
    """Apply the fitted nonlinearity to the analyzed image and fit the corrected image in turn.

    Returns (original_ext, corrected_ext, results_corrected): the original and corrected frames
    in electrons (lists of 2-D arrays, one per extension) and run_pipeline's dict for the
    corrected image. The corrected image is saved in the run's output folder, as with
    --correct_nonlinearity.
    """
    from . import run_nonlinearity_studies as run
    from .nonlinearity_correction import correct_nonlinearity_ext

    corrected_ext, corrected_path = correct_nonlinearity_ext(
        results['full_data_ext'], results['pedestals'], results['gains'], results['parabola_coeffs'],
        results['fits_path'], results['fig_path'],
        data_nl=results['charge_minus_npeak_ext'] if args.nonlinearity_correction_model == 'data' else None)
    original_ext = [convert_to_electrons(data, pedestal, gain, flatten=False)
                    for data, pedestal, gain in zip(results['full_data_ext'], results['pedestals'], results['gains'])]

    # A corrected image is already pedestal-subtracted and in electrons: the pipeline recognizes it.
    args_corrected = copy.copy(args)
    args_corrected.file_string, args_corrected.stitch_fits = str(corrected_path), False
    _set_plot_switches(args_corrected, ())
    return original_ext, corrected_ext, run.run_pipeline(args_corrected)


def plot_nonlinearity_comparison(args, results, results_corrected, **kwargs):
    """Nonlinearity curve of the original image (blue points, red fit) with the corrected image's on the
    same axes (purple points, orange fit), one figure per extension, styled by the run's config."""
    from .cli_config import _parse_lim
    from .plotting import plot_nonlinearity

    n_ext = len(results['gains'])
    plot_nonlinearity(results['peaks_ext'], results['parabola_coeffs'], results['peak_charge_e_ext'],
                      results['charge_minus_npeak_ext'], results['fit_range_right_ext'],
                      overlay=[results_corrected[key] for key in ('peak_charge_e_ext', 'charge_minus_npeak_ext', 'parabola_coeffs')],
                      xlim=_parse_lim(args.plot_nonlinearity_xlim, n_ext=n_ext),
                      ylim=_parse_lim(args.plot_nonlinearity_ylim, n_ext=n_ext),
                      individual_figsize=tuple(args.plot_nonlinearity_individual_figsize),
                      additional_title=args.extra_plot_title, show_titles=args.show_titles,
                      suptitle='Nonlinearity Curve' if args.default_titles else None,
                      nimages=results['nimages'], line_color='r', scatter_color='b', s=6, alpha=1,
                      sharey=args.plot_nonlinearity_sharey, plot_individual=True, **kwargs)
