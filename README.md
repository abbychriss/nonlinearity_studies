# Nonlinearity Studies

A Python package for analyzing nonlinearity in CCDs.

## Overview

- **Image Stitching**: Combine multi-extension FITS images row-wise and run analysis on stitched image (via the `analysis_tools` package)
- **Pedestal Subtraction**: Independently computes and subtracts pedestal along specified axis (using my pedestal_subtract package)
- **Noise & Gain Calculation**: Determines noise and gain from data and converts to e-
- **Peak Finder**: Finds all electron peaks in charge distribution
- **Nonlinearity Computation**: Quantifies nonlinearity as a function of charge by fitting nonlinearity curve to parabola
- **Single-Electron Resolution**: Quantifies single electron resolution at a given charge *q* by fitting *n* Gaussians in window around *q* and reporting the shared peak width σ (e-), goodness of fit, and a resolved/marginal/unresolved verdict
- **Visualization**: Plots zero-one electron peaks in ADU and e-, histograms of charge distribution at any q with peak labels, and nonlinearity curves. Can be plotted per extension or together as 2×2 subplot. Control figure sizes, axis sharing, scales, limits, titles, etc. using config.
- **Meta analysis**: Plot results of nonlinearity analysis (noise, gain, resolution, nonlinearity) as a function of an image parameter using the summary CSV tables produced during individual analyses (now provided by the `analysis_tools` package via the `meta-studies` command)

## Code structure

The analysis library is organized by responsibility rather than as one large module:

- `fit_data.py` — fit models and primitives (parabola, Gaussians, peak finding)
- `fit_range.py` — automatic estimation of the nonlinearity fit range
- `nonlinearity.py` — nonlinearity computation per extension / at a given charge
- `nonlinearity_correction.py` — applies the fitted nonlinearity model to an image and saves the corrected FITS
- `resolution.py` — single-electron resolution fitting and classification
- `plotting.py` — image, peak, nonlinearity, and resolution figures
- `summary.py` — per-extension and wide summary tables (CSV)
- `run_nonlinearity_studies.py` — command-line interface; `cli_config.py` holds its config/argument helpers. `run_pipeline(args)` runs the whole analysis and returns the intermediate results as a dict (`main` is a thin wrapper around it that returns nothing, so the `run-nonlinearity-studies` console script exits cleanly)
- `notebook.py` — helpers for the example notebook `scripts/nonlinearity_analysis.ipynb` (see [Example notebook](#example-notebook))
- `nonlinearity_studies.py` — backwards-compatibility shim that re-exports the public API, so existing `from nonlinearity_studies.nonlinearity_studies import ...` imports keep working
- `gain_independent_nonlinearity/` — the gain-independent analysis, which leaves the spectrum in ADU instead of calibrating it with the zero/one gain (see [Gain-independent nonlinearity](#gain-independent-nonlinearity))

Zero/one-electron peak fitting, pedestal subtraction, and FITS I/O are provided by the `pedestal_subtract` package; FITS stitching and meta-analysis plotting by `analysis_tools`.

## Installation

### Clone the repository

Navigate to your root project directory then run:
```bash
git clone https://github.com/abbychriss/nonlinearity_studies.git
cd nonlinearity_studies
```

### Install package 

#### Using Conda environment.yml

The recommended way to install is to use the `environment.yml` file to create a conda environment called `nonlinearity-studies` configured with the right version of Python and required dependencies. It also pip installs the nonlinearity_studies package in development mode. To set up environment, run:

```bash
conda env create -f environment.yml
conda activate nonlinearity-studies
```

#### Direct pip install in development mode

Install package in development mode (the -e flag makes files editable):

```bash
pip install -e .
```

####  regular installation

```bash
pip install .
```

### Requirements for direct pip installation

- Python >= 3.8
- numpy >= 1.26.0
- 3.11 > matplotlib >= 3.8.0
- scipy >= 1.11.0
- astropy >= 5.3
- tqdm
- [`pedestal_subtract`](https://github.com/abbychriss/pedestal_subtract) (installed automatically from git)

## Usage

### Command-line interface

There are three ways to run the analysis:

**1. As a console script** (after pip installation):

```bash
run-nonlinearity-studies [OPTIONS] <file_string>
```
You can also put `file_string` and any option below in a JSON file and run with `-j path/to/config.json`. Every boolean flag has a matching `--no-<flag>` variant for overriding a JSON-config `true` from the command line.

**2. As a Python module** (works during development):
```bash
python -m nonlinearity_studies.run_nonlinearity_studies [OPTIONS] <file_string>
```

**3. As a direct executable** (requires executable permissions - should be activated automatically but if not run `chmod +x nonlinearity_studies/run_nonlinearity_studies.py`):

```bash
./nonlinearity_studies/run_nonlinearity_studies.py [OPTIONS] <file_string>
```

### Gain-independent nonlinearity

The standard analysis calibrates the whole spectrum with one gain measured from the zero and one electron peaks, then reports how far each peak's charge falls from its peak number. That assumes a constant gain, which is the very thing the nonlinearity measurement contradicts. The `gain_independent_nonlinearity` package runs the same FITS loading, pedestal subtraction and zero/one fit, but histograms the pedestal-subtracted spectrum in **ADU** (`do_convert_to_electrons=False`) and derives the nonlinearity without ever dividing by that gain.

```bash
gain-independent-nonlinearity -j config/gain_indep_nonlinearity.json
python -m nonlinearity_studies.gain_independent_nonlinearity.run [OPTIONS] <file_string>
```

It takes the same JSON config mechanism and the same file/pedestal-subtraction/zero-one/peak-finding and text-style (`font_size`, `legend_font_size`, `font_family`, `use_tex`) options as `run-nonlinearity-studies`, plus its own. Peak finding here is **two steps at two binnings**, with one config key each:

- **Coarse** — `scipy.signal.find_peaks`, keyed on the shared `--nonlinearity_peakfinder_bin_factor` (default `10`), in bins per **electron**, exactly as in the standard analysis. This spectrum is histogrammed in ADU, so the number is converted to bins per ADU (divided by `ADU_PER_E`, a module-level constant in `gain_independent_nonlinearity/run.py`, nominally `1.3`) for the histogram `find_peaks` actually runs on — the physical bin width is the same as upstream either way. `--peak_finder_widths` (electrons) and `--peak_finder_buffers` (normalized to 10 bins/e-) are converted the same way, so all three knobs mean what they mean in the standard analysis.
- **Sub-bin** — the log-parabola refinement of each of those peaks, keyed on `--sub_bin_peakfinder_bin_factor` (default `40`), in bins per **ADU**. This step has no counterpart in the standard analysis, so it is keyed in the units the spectrum is actually in, and it re-bins the same ADU spectrum over the same range rather than reusing the coarse histogram. `--sub_bin_half_width` counts bins of *this* histogram.

`--peak_finder_buffers` is applied in one place only, the coarse step's `find_peaks` distance filter; the sub-bin step has no distance filter. `ADU_PER_E` never rescales the data, and is a fixed number rather than the per-extension zero/one gains so the peak finder stays independent of the fit. The all-peaks legend names both factors, each in `bins/ADU` with its `bins/e-` equivalent beside it (computed from that extension's real gain, display only), and — when `--plot_all_peaks_window_bin_factor` re-histograms the displayed window — a third entry for the binning the bars are actually drawn at.

The all-peaks plot also overlays the sub-bin fits themselves, each drawn over exactly the window it was fit on: **purple** for the log-parabola fits and **green** for the peaks that fell back to the seeded Gaussian refit (`--sub_bin_fit_alpha`, default `0.6`). Their legend entries count how many peaks each method localized, e.g. `Log parabola fit (1436 / 1553 peaks found)` and `Gaussian fits (117 / 1553 peaks found)`. The fitted amplitudes are rescaled from the sub-bin histogram's bin width to the drawn one, so the curves sit on the bars whatever binning is displayed. Writing q(n) for the ADU charge of peak number n (n = 0 is the zero-electron peak) and g = q(1) - q(0) for the standard gain, it computes:

| quantity | module | meaning |
| --- | --- | --- |
| `delta_g(n) = q(n) - q(n-1) - g` | `deltas.py` | differential gain, with a fitted line through the scatter |
| `I(n) = sum delta_g`, `S(n) = I(n)/n` | `deltas.py` | the nonlinearity in ADU; the slope of S(n) is the nonlinearity parameter alpha' |
| `Q(n) = a0 + g_fit*n + alpha*n^2` | `spectrum_fit.py` | parabola (and optional cubic) fit to the charge spectrum; `G = g_fit + alpha` is the fit's own gain |
| `g_fit(i)` | `spectrum_fit.py` | the same fit using only the first i peaks — the average gain over that range — and its residual `g_fit(i) - g` |
| `g_dir(n) = (q(n) - q(0))/n` | `direct_gain.py` | the gain read straight off peak n, and `delta_dir(n) = g - g_dir(n)` |
| gain scan | `gain_scan.py` | recomputes the **standard** nonlinearity with every candidate gain `g_fit(i)` and reports the one minimizing \|nonlinearity\| at `--gain_indep_scan_eval_charge` |
| gain comparison | `compare_nonlinearity_methods.py`, `plot_nonlinearity_curves.py` | overlays the **standard** nonlinearity curve computed with every gain this analysis produces (zero/one, sub-bin, `g_fit(i)` at each `--gain_indep_compare_test_i` and at the scan's optimum, and `g_dir`); also sweeps a bare calibration gain over [1.0, 1.4] at fixed charge, and tracks alpha / alpha' against the peak number they were fit up to, with `g_fit(i)` on a twinned axis |

Peak numbering is anchored on the tallest detected peak, so a low-count bump in the pedestal's left tail cannot take index 0 and shift every peak number by one. `--gain_indep_fit_n_min` (default 2) is where the running fit starts — 2 is a two-point line fit, which reproduces the standard gain exactly, so the scan covers the whole space including the standard method — and `--gain_indep_fit_n_max` caps how far up the spectrum any fit reaches. Left unset it is read off the data per extension by the same changepoint detector the standard analysis uses to place `fit_range_right`, run on the raw ADU spectrum (`x = n`, `y = q(n)`); both of the detector's stages difference out the linear and quadratic parts of the curve, so no gain division or line subtraction is needed and the estimate stays gain independent. The chosen domain is reported as `fit_n_max` in the summary. Wherever a plot draws a fitted curve, a short panel beneath the scatter carries the residuals about that curve, sharing the x axis, in both the individual and the 2x2-together figures. Each calculation has its own `--plot_<name>_individual` / `--plot_<name>_together` pair; they share the `plot_gain_indep_nonlin_*` figure-size and shared-axis settings. The comparison plots are keyed the same way, as `--plot_compare_nonlinearity_curves_*`, `--plot_alpha_alpha_prime_*`, and `--plot_nonlin_vs_gain_*`; the gain sweep behind the last of those is only run when that plot is switched on, since nothing else consumes it. With `--save_output`, every scalar lands in `gain_independent_summary.csv` next to the plots.

### Example notebook

`scripts/nonlinearity_analysis.ipynb` walks through the pipeline one stage at a time: the stitched image, hot columns and zero/one calibration, the all-peaks spectrum, the nonlinearity curve, and the optional nonlinearity-corrected image and gain-independent plots. Every plot is shown in a panel with one tab per extension.

Set three variables at the top:

- `base_config_file` — a JSON config (a path, or a file name in `config/`). Every option comes from it, for both pipelines.
- `images` — the image(s) to analyze; accepts glob wildcards.
- `stitch_fits` — if `True`, all the files `images` matches are stitched into one image first. If `False`, `images` must match exactly one file (already stitched, or a single image) and the notebook raises an error otherwise; nothing is stitched, but the image is still plotted under the stitching heading.

A few config options are overridden because the notebook displays the plots itself: only the per-extension (`..._individual`) plots are drawn, `show_plots` and `save_output` are off (so nothing is written beyond the stitched file, the pedestal-subtraction cache and the corrected image), `verbose` is on, `extra_plot_title` is blanked, and `file_string`/`stitch_fits` come from the notebook. The all-peaks section has two boxes for the x-limits (in electrons) and a Plot button that redraws the plot for that window; it starts at the config's `plot_all_peaks_xlim`. The two zero/one-peak sections (ADU and electrons) each have a slider for the number of histogram bins (`zero_one_n_bins`, starting at the config's value): each position repeats the double-Gaussian fit at that binning and redraws, for looking only — the pedestals and gains used by the rest of the notebook are those of the run. The notebook runs `run_pipeline` and the gain-independent `main` as they are and collects the figures they draw (`notebook.capture_figures`), so it stays in step with the command-line interface.

### As a Python module

```python
from nonlinearity_studies import (
    get_fits,
    pedestal_subtract_ext_cached,
    get_zero_one_peaks_ext,
    get_all_peaks_ext,
    get_nonlinearity_ext,
    plot_nonlinearity,
)

fits_path = "path/to/stitched.fits"

# Load all 4 extensions
data_ext = get_fits(fits_path)

# (Optional) pedestal-subtract with FITS caching
data_ext = pedestal_subtract_ext_cached(
    data_ext, source_path=fits_path,
    n_std_to_mask=1.5, axis="row",
    use_biweight_loc=True, use_biweight_midvar=True,
)

# Fit zero/one peaks to get pedestals, gains, double-Gaussian params per extension
zero_one_counts_ext, zero_one_edges_ext, pedestals, gains, double_gauss_popts, zero_one_ranges = \
    get_zero_one_peaks_ext(data_ext)

# Find every electron peak in each extension
counts_ext, edges_ext, peaks_ext, centers_ext, hist_ranges, data_electrons_ext = get_all_peaks_ext(
    data_ext, widths=0.9, buffers=[3, 3, 3, 3],
    pedestals=pedestals, double_gauss_popts=double_gauss_popts, gains=gains,
)

# Fit the parabolic nonlinearity model
peak_charge_e_ext, charge_minus_npeak_ext, parabola_coeffs, parabola_pcovs = \
    get_nonlinearity_ext(peaks_ext, centers_ext, pedestals, gains, fit_range_right_ext=500)

plot_nonlinearity(
    peaks_ext, parabola_coeffs, peak_charge_e_ext, charge_minus_npeak_ext,
    fit_range_right_ext=500,
)
```


## Examples

In the following examples, I will assume you are working from the base `nonlinearity_studies` directory. For the first example, we will fit the zero and one electron peaks for a single 250x3500 image with 1x1 binning and 500 skips. Navigate to project directory in terminal and run:

```bash
run-nonlinearity-studies \
    "examples/images/VR-3/avg_img_CV_250x3500x500_bin1x1_125_20260317_213403_0.fz" \
    --plot_zero_one_ADU_together
```

Next, let's stitch 10 images together from `examples/images/VR-3` and run the analysis script on the stitched output:

```bash
run-nonlinearity-studies \
    "examples/images/VR-3/avg*.fz" \
    --stitch_fits \
    --plot_zero_one_ADU_together \
    --plot_nonlinearity_together
```

Now every time we want to analyze the stitched image again, we can pass the stitched image directly into the script (and drop `--stitch_fits`). The pedestal-subtracted result is cached next to the source as `<stem>.pedsub.fits` and reused on subsequent runs if the pedsub params match. Run the stitched image and save the plots:

```bash
run-nonlinearity-studies \
    "examples/images/VR-3/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits" \
    --plot_zero_one_ADU_together --plot_zero_one_electrons_together \
    --plot_nonlinearity_together \
    --save_output
```

If we want only the nonlinearity at specific charge values:

```bash
run-nonlinearity-studies \
    "examples/images/VR-3/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits" \
    --get_nonlinearity_at 10 50 500 1000
```

Note that the parabola fit range per extension is chosen automatically using the `changepoint` method, which cross-checks each result against `var(a)` to assign a per-extension confidence label. You can force a different choice for fit range right per extension by:

```bash
run-nonlinearity-studies \
    "examples/images/VR-3/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits" \
    --fit_range_right 600 1000 800 1000  --plot_nonlinearity_together --save_output
```

With `--save_output`, the full per-extension estimate (chosen value, cross-check, confidence) is written to `fit_range_estimate.json` in the run's output folder. Add `--verbose` to also print the confidence label per extension and a warning for any that need review.

To quantify single-electron resolution at one or more charges, pass `--resolution_at`. This prints a per-charge/per-extension table of σ, reduced χ², ΔAIC, peak counts, and verdict, and — with `--save_output` — merges the same quantities into `extension_summary.csv` as `<field>_at_<charge>_e` columns (e.g. `sigma_e_at_500_e`) alongside gain, noise, and nonlinearity, one row per extension. With a plot flag it also draws the charge distribution in the window with Gaussian fits:

```bash
run-nonlinearity-studies \
    "examples/images/VR-3/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits" \
    --resolution_at 250 500 1000 --resolution_window 10 \
    --plot_resolution_together --save_output
```

Or can use the json config in `examples/config/resolution_study.json`:

```bash
run-nonlinearity-studies \
    "examples/images/VR-3/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits" \
    -j examples/config/resolution.json
```

Now say you are studying the effect of a variable such as deltaV on nonlinearity. Then the script `meta_studies.py` in the `scripts` directory can compare the results of those studies. First run the nonlinearity studies to get an `extension_summary.csv` table for each deltaV, either by changing the image parameter and titles in `default_nonlinearity.json` appropriately and running

```bash
run-nonlinearity-studies -j examples/config/nonlinearity.json
```

or changing them on the command line

```bash
run-nonlinearity-studies "examples/images/VR-6_standard_deltaV/*" --stitch_fits -j "config/default_nonlinearity.json" --extra_plot_title "VR = -6, standard deltaV"
```
```bash
run-nonlinearity-studies "examples/images/VR-6_increased_deltaV/*" --stitch_fits -j "config/default_nonlinearity.json" --extra_plot_title "VR = -6, increased deltaV"
```

Then you can plot the noise/gain/nonlinearity versus deltaV per extension using `config/VR6_deltaV_vs_ext_study.json` by running

```bash
python scripts/meta_studies.py -j examples/config/VR6_deltaV_study.json
```
or compare these values versus extension

```bash
python scripts/meta_studies.py -j examples/config/VR6_deltaV_vs_ext_study.json
```

To compare the resolution at different deltaV:

```bash
run-nonlinearity-studies \
    "examples/images/VR-6_increased_deltaV/stitched-fits/cds_avg_img_L2_250x3500x500_1x1_L2_125_10_stitched.fits"  \
    -j examples/config/resolution.json \
run-nonlinearity-studies \
    "examples/images/VR-6_standard_deltaV/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits"  \
    -j examples/config/resolution.json \
python scripts/meta_studies.py -j examples/config/VR6_deltaV_study.json
```

### `meta_studies.py` config options

`meta_studies.py` is a general purpose script to compare the results of multiple datasets in a study, driven by a JSON config. For each dataset that you are comparing, you must specify an `x` value (the parameter the dataset represents) and a `table` (the path to a csv file containing the results, which may contain glob wildcards, e.g. `plots/avg*/extension_summary.csv`). `table` may instead be a **list** of such paths to attach several tables to one `x` value (typically averaged together — see `use_average`). A dataset may also carry `"exclude_ext": [...]` to drop specific extensions from that dataset only.

Required keys: `study_name`, `x_label`, `quantities`, `datasets`. Everything else is optional with the defaults shown below.

Common options (apply to both use cases):

- `study_name`: short name used in output filenames and default titles.
- `x_label`: axis label for the varied parameter; a `"NAME (UNIT)"` form (e.g. `"VR (V)"`) attaches the unit to values in titles.
- `quantities`: map of CSV column → spec. Each spec may set `ylabel`, `title`, `title_vs_ext`, `fit_line`, `plot_fit_line`, `fit_legend` and `connect_points` (per-quantity overrides), `xscale`/`yscale` (matplotlib axis scale, e.g. `"linear"` (default) or `"log"`), `xlim`/`ylim` (a `[min, max]` pair fixing that axis range; either entry may be `null` to leave that side autoscaled), and `column` (absolute 0-based CSV index, to pull a column regardless of its header name).
- `datasets`: list of `{"x": ..., "table": ..., (optional) "series": ..., "exclude_ext": [...]}`. `table` is a glob string or a list of glob strings.
- `x_axis`: `"value"` (default) plots the varied parameter on the x-axis, one line per extension; `"extension"` flips it so extension is the x-axis, each series a coloured line, and the x value pinned in the title (one figure per distinct x).
- `invert_x` (default `false`): reverse the x-axis (e.g. high VR → low VR).
- `output_dir` (default: lowest common parent directory of the table paths, when omitted or `null`): where figures, the report, and fit-stats CSV are written.
- `colors` (default Tableau 10): list (cycled per extension) or `{ext: color}` map.
- `dpi` (default `350`), `show_plots` (default `true`).
- `save_output` (default `true`): master switch for writing anything to disk (figures, the text report, and the fit-stats CSV). Set it `false` to compute and preview without leaving files behind — the summary table is still printed and plots still appear on screen per `show_plots`.
- `save_report` (default `true`): when `save_output` is on, also write the text report and fit-stats CSV (not just the figures).
- `use_average` (default `false`): collapse rows that share the same `(x, series, ext)` — whether listed under one dataset's `table` or spread across datasets — into their mean (drawn as the opaque marker). The contributing raw values and their standard deviation are retained so the spread can be shown (see `show_multiple_values` / `show_error_bars`).
- `show_multiple_values` (default `true`): when averaging pooled several tables, overlay each actual table value as a transparent marker at the same x (with a faint vertical line spanning their range, per `error_bar_line`). Points from a single table add nothing.
- `show_error_bars` (default `false`): instead of (or in addition to) the actual values, draw mean ± std error bars with transparent caps. Kept as a separate option from `show_multiple_values`.
- `error_bar_line` (default `true`): whether the faint vertical line is drawn for whichever spread is shown (`show_multiple_values` range and/or `show_error_bars` bar). Set `false` to keep just the transparent markers/caps without the connecting line.
- `plot_individual` (default `true`): one figure per quantity. `plot_quantities_together` (default `false`): all quantities as subplots in one figure. `plot_ext_together` (default `false`): one figure per quantity with one subplot per extension, all series overlaid in each subplot, or one figure per (quantity, series) when `series_separate_plot` splits the series (always a roughly square grid, e.g. 2x2).
- `individual_figsize` (default `[9, 6.5]`), `together_subplots_figsize` (default `[13, 10]`), `together_subplots_ncols` (default: roughly square grid).
- `fit_line` (default `false`): global toggle for fitting `y = slope·x + intercept` per line (percentage change is always reported).
- `plot_fit_line` (default `true`): draw the fitted line on the plots. Set `false` to keep the fit in the stats report only. Only matters when `fit_line` is on.
- `fit_legend` (default `true`): give the fitted line a legend entry (`… fit ($R^2 = …$)`). Set `false` to draw the line without one. If `plot_fit_line` is `false` and `inset_legend` is `false`, the marker's legend entry shows `$R^2 = …$` in place of `EXT<n>` instead, and in the subplot grids that quantity's cell keeps its own legend (left out of the shared one).
- `connect_points` (default `null`): draw a plain line joining adjacent data points (separate from the fit line). Left `null` it follows the historical behaviour — points are connected only when they aren't being fit — so set it `true` to also connect points under a fit line, or `false` to show bare markers. Honoured globally or per quantity.

Series options (how a second dimension is split into separate lines):

- `series_label`: human-readable name for the series dimension, used in legends/titles.
- `series_same_plot` (default `true`): overlay all series on one figure vs. split onto one figure per series.
- `suppress_series_plot` (default `false`): hide the series distinction in the plots only — every series is drawn on one shared figure with a uniform style and a single `EXT{n}` legend entry per extension (overriding `series_same_plot` and `uniform_series_style` for plotting). The series are still kept in the printed table, fits, and report.
- `uniform_series_style` (default `false`): when split, force every per-series figure to use the same solid line and circle marker so the figures don't look gratuitously different.
- `series_plot_styles`: matplotlib line options per series — either a `{series: {option: value}}` map (plus an optional shared `"default"` entry applied to every series first) or a list of option dicts cycled over the series in order. Any `Line2D` keyword works (`linestyle`, `marker`, `alpha`, `linewidth`, `markersize`, ...), and anything matplotlib won't accept is rejected when the config loads. Unset options fall back to the built-in look: a linestyle *and* a marker cycled per series, so series stay separable even when a single x value makes the linestyle invisible.

**Use case 1 — `extension_summary.csv` (vary a parameter across runs).** The series is a per-dataset tag (`"series"`), constant across each table — e.g. standard vs. increased deltaV at the same VR. See `config/VR6_deltaV_study.json`:

```json
{
  "study_name": "VR6_deltaV",
  "x_label": "VR (V)",
  "invert_x": true,
  "series_label": "deltaV",
  "series_same_plot": true,
  "series_plot_styles": {
    "Standard": {"linestyle": "-", "marker": "o"},
    "Increased": {"linestyle": "--", "marker": "s"}
  },
  "quantities": {
    "gain_adu_per_e": {"ylabel": "Gain (ADU/e-)", "title": "Gain vs VR"},
    "noise_e": {"ylabel": "Noise (e-)", "title": "Noise vs VR"},
    "nonlinearity_at_500_e": {"ylabel": "Nonlinearity at 500 e- (e-)", "title": "Nonlinearity at 500 e- vs VR"}
  },
  "datasets": [
    {"x": -6, "series": "Standard",  "table": "./examples/images/VR-6_standard_deltaV/plots/avg*/extension_summary.csv"},
    {"x": -6, "series": "Increased", "table": "./examples/images/VR-6_increased_deltaV/plots/cds*/extension_summary.csv"}
  ]
}
```

**Use case 2 — resolution at several charges (one extension_summary.csv, one quantity per charge).** Resolution results are merged into `extension_summary.csv` as `<field>_at_<charge>_e` columns (e.g. `sigma_e_at_500_e`), the same way `nonlinearity_at_<charge>_e` already works, so no `series_column` is needed for this case -- just list one quantity per (field, charge) pair. `plot_mega` then lays them out as a grid (one column per quantity) for an at-a-glance comparison across charges. See `config/VR-8_VDD_resolution_study.json`:

```json
{
  "study_name": "VR8_VDD_resolution",
  "x_label": "VDD (V)",
  "invert_x": true,
  "series_label": "Charge (e-)",
  "series_same_plot": false,
  "uniform_series_style": true,
  "plot_quantities_together": false,
  "quantities": {
    "sigma_e_at_500_e": {"ylabel": "Resolution sigma (e-)", "title": "Resolution sigma at 500 e- vs VDD"},
    "sigma_e_at_1000_e": {"ylabel": "Resolution sigma (e-)", "title": "Resolution sigma at 1000 e- vs VDD"},
    "sigma_e_at_1200_e": {"ylabel": "Resolution sigma (e-)", "title": "Resolution sigma at 1200 e- vs VDD"},
    "reduced_chi2_at_500_e": {"ylabel": "Reduced $\\chi^2$", "title": "Fit reduced $\\chi^2$ at 500 e- vs VDD"}
  },
  "datasets": [
    {"x": -19.0, "table": ".../VR-8_VDD-19/.../extension_summary.csv"},
    {"x": -20.5, "table": ".../VR-8_VDD-20.5/.../extension_summary.csv"},
    {"x": -21.5, "table": ".../VR-8_VDD-21.5/.../extension_summary.csv"}
  ]
}
```

`series_column` (split a table into one series per distinct value of a given CSV column) is still available for tables of your own that vary a second parameter within one file -- it's just no longer how resolution-at-charge tables are read.

### JSON config

It's easier to control all options using a JSON file instead of changing parameters on the command line. Check out the `/config` directory for examples of JSON files for different use cases, making sure to change any hard-coded file names that are not referenced in `\examples`. A config for `run-nonlinearity-studies` looks like:

```json
{
  "file_string": "examples/images/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits",
  "stitch_fits": false,
  "stitch_fits_axis": "vertical",
  "save_output": true,
  "output_dir": null,
  "show_plots": true,

  "plot_charge_per_column_together": false,
  "plot_charge_per_column_individual": false,
  "plot_zero_one_ADU_individual": false,
  "plot_zero_one_ADU_together": true,
  "plot_zero_one_electrons_individual": false,
  "plot_zero_one_electrons_together": true,
  "get_nonlinearity_at": 500,

  "resolution_at": [200, 600, 1000],
  "resolution_window": 10,
  "resolution_sigma_well": 0.25,
  "resolution_sigma_limit": 0.5,
  "resolution_min_peak_frac": 0.6,
  "plot_resolution_individual": false,
  "plot_resolution_together": true,
  "plot_resolution_individual_figsize": [8, 6.5],
  "plot_resolution_subplots_figsize": [13, 9],
  "plot_resolution_sharex": false,
  "plot_resolution_sharey": false,

  "do_pedestal_subtraction": true,
  "n_std_to_mask": 1.5,
  "pedsub_max_iter": 5,
  "pedestal_subtraction_axis": "row",
  "use_biweight_loc": true,
  "use_biweight_midvar": true,
  "pedsub_cache_dir": null,
  "force_pedsub": false,

  "peak_finder_widths": 0.9,
  "peak_finder_buffers": [3, 3, 3, 3],
  "peak_finder_prominences": null,
  "fit_range_right": "auto",
  "correct_nonlinearity": false,
  "nonlinearity_correction_model": "parabola",
  "nonlinearity_fit_fix_intercept": true,
  "nonlinearity_fit_known_points": null,

  "plot_zero_one_individual_figsize": [8, 6],
  "plot_zero_one_subplots_figsize": [12, 8],
  "plot_zero_one_yscale": "linear",
  "plot_zero_one_xlim": "default",
  "plot_zero_one_ylim": "default",
  "plot_zero_one_ylim_electrons": "default",
  "electron_fit_mode": "transform",
  "plot_zero_one_sharex": true,
  "plot_zero_one_sharey": true,

  "plot_all_peaks_individual": false,
  "plot_all_peaks_together": true,
  "plot_all_peaks_individual_figsize": [7, 6],
  "plot_all_peaks_subplots_figsize": [9, 7],
  "plot_all_peaks_xlim": [1000, 1050],
  "plot_all_peaks_ylim": [0, 40],
  "plot_all_peaks_yscale": "linear",
  "plot_all_peaks_sharex": true,
  "plot_all_peaks_sharey": true,

  "plot_nonlinearity_individual": false,
  "plot_nonlinearity_together": true,
  "plot_nonlinearity_individual_figsize": [6, 5],
  "plot_nonlinearity_subplots_figsize": [9, 7],
  "plot_nonlinearity_xlim": [-50, 1300],
  "plot_nonlinearity_ylim": [-40, 5],
  "plot_nonlinearity_sharex": true,
  "plot_nonlinearity_sharey": true,

  "show_titles": true,
  "extra_plot_title": "VR=-7: ",
  "nimages": 10
}
```

Then run:

```bash
run-nonlinearity-studies -j config/default_nonlinearity.json
```

Explicit command-line arguments override JSON values, so this also works:

```bash
run-nonlinearity-studies -j config/default_nonlinearity.json --no-save_output
```

## Options

##### Pipeline / I/O
- `-j`, `--json PATH`: Load command-line arguments from a JSON config file (CLI arguments override JSON values)
- `-f`, `--stitch_fits`: Stitch multi-extension FITS files before analysis
- `--stitch_fits_axis`: `vertical` (default, stack along y) or `horizontal` (stack along x). Only used with `--stitch_fits`.
- `-s`, `--save_output` / `--no-save_output`: Save all output — generated plots as JPEGs and the per-extension summary table as `extension_summary.csv` in the output directory
- `-o`, `--output_dir DIR`: Directory for all saved output (plots, summary CSV, config snapshot). Defaults to a `plots/` folder alongside the source FITS.
- `--show_plots` / `--no-show_plots`: Display plots interactively via `plt.show()` (default `true`). When disabled, figures are simply not shown (and are closed to free memory) instead of switching matplotlib to a non-interactive backend — useful for headless/batch runs, and still works alongside `--save_output`.
- `-v`, `--verbose`: Print verbose output
- `--nimages N`: Number of stitched images (used for plot labeling). Auto-detected from filenames matching `_N_stitched` or `_N_<axis>_stitched` (the name `stitch_fits` gives its output)
- `--extra_plot_title TITLE`: Additional title text prepended to every plot title

##### Plotting
Each plot type is controlled by its `_individual` and `_together` toggles (see [Plot layout](#plot-layout-per-plot-type)) — a plot is generated whenever either is `true`.
- `--plot_charge_per_column_together` / `--plot_charge_per_column_individual`: Plot the median charge per column for each extension (a combined 2×2 grid, or one figure per extension) on the **raw**, pre-pedestal-subtraction data — a diagnostic for anomalous columns. A column whose median sits ≥ `--n_std_to_mask` biweight-SDs from the extension's biweight location is flagged "hot" (drawn red). Every column is plotted; columns excluded by `--fit_cols` are shaded light grey. `--plot_charge_per_column_figsize W H` sets the combined figure size (default `13 9`); `--plot_charge_per_column_individual_figsize W H` sets the individual figure size (default `8 6`).
- `--plot_zero_one_ADU_individual` / `--plot_zero_one_ADU_together`: Render the zero/one electron peak fits in ADU, one figure per extension / as a combined 2×2 grid. The two are independent — set either, both, or neither.
- `--plot_zero_one_electrons_individual` / `--plot_zero_one_electrons_together`: Same, but in electron units. Fully independent of the ADU flags above, so ADU and electron-units plots, and individual vs. combined layout, can each be toggled on their own.
- `-g`, `--get_nonlinearity_at CHARGE...`: Evaluate the nonlinearity polynomial at one or more charge values
- `-r`, `--resolution_at CHARGE...`: Quantify single-electron resolution at one or more charge values (see [Single-electron resolution](#single-electron-resolution))

##### Per-extension summary table
Every run prints a per-extension summary table to stdout (no flag required) with the **gain** in ADU/e- (the separation between the zero- and one-electron peaks from the double-Gaussian fit), the **noise** in e- (the standard deviation of the zero-electron peak, converted from ADU to electrons by dividing by the gain), and one **nonlinearity** column per charge passed to `--get_nonlinearity_at` (falling back to 500 e- when none is given). With `--save_output`, the same table is written as a CSV (`extension_summary.csv`) in the output directory — one row per extension with a header row, also containing the parabola coefficients `nl_a`, `nl_b`, `nl_c` (nonlinearity NL(q) = a·q² + b·q + c in e-, per extension) so the correction can be re-evaluated at any charge, ready to load with `numpy.genfromtxt(..., delimiter=',', names=True)` or `numpy.loadtxt(..., delimiter=',', skiprows=1)`.

##### Pedestal subtraction
- `--do_pedestal_subtraction` / `--no-do_pedestal_subtraction`: Toggle pedestal subtraction (default `true`)
- `--n_std_to_mask FLOAT`: Mask threshold (in standard deviations) for the pedestal estimator
- `--pedsub_max_iter INT`: Max sigma-clip iterations for pedestal subtraction (default `5`). Each pass re-estimates the pedestal from the clipped zero-peak core and stops early once it converges.
- `--pedestal_subtraction_axis {row, col, row_then_col, col_then_row}`: Axis to compute the pedestal across
- `--use_biweight_loc` / `--use_biweight_midvar`: Use Tukey biweight location/midvariance instead of mean/std
- `--pedsub_cache_dir DIR`: Directory for the pedestal-subtracted FITS cache (default: alongside the source FITS)
- `--force_pedsub`: Recompute pedestal subtraction even if a matching cache exists

##### Fitting
- `--peak_finder_widths W [W ...]`: Minimum peak width required by `scipy.signal.find_peaks`, in **electrons** (internally multiplied by `nonlinearity_peakfinder_bin_factor` to convert to bins). Scalar or one per extension. Larger = filters narrow noise spikes more aggressively.
- `--peak_finder_buffers B [B ...]`: Buffer SUBTRACTED from `nonlinearity_peakfinder_bin_factor` to compute the minimum neighbor-peak distance: `d = nonlinearity_peakfinder_bin_factor - buffer`, in bins. The value is quoted **normalized to 10 bins/e-** and scaled to whatever binning is actually in use (`buffer * nonlinearity_peakfinder_bin_factor / 10`), so retuning the bin factor keeps the same physical spacing, and the same number means the same peak separation in the gain-independent module (see above). At 10 bins/e-: `buffer=0` → 10 bins = 1 electron spacing (physical), `buffer=3` → 7 bins ≈ 0.7 electron (loose), `buffer=-2` → 12 bins = 1.2 electron (strict). Larger buffer = looser; smaller/negative = stricter. Fractional values are allowed (e.g. `3.5`).
- `--peak_finder_prominences P [P ...]`: Minimum peak prominence in **histogram counts** (same units as the y-axis of `plot_all_peaks`). Scalar, one per extension, or `null` to disable. Often the most robust filter — measures how far a peak sticks up above its surrounding baseline, regardless of width.
- `--nonlinearity_peakfinder_bin_factor N`: Number of histogram bins per electron in the all-peaks histogram (default `10`). Also drives the `peak_finder_widths` electron-to-bin conversion (`width_in_bins = width_in_electrons * nonlinearity_peakfinder_bin_factor`) and the buffer math (`distance_in_bins = nonlinearity_peakfinder_bin_factor - buffer`). Higher = finer histogram, but small-width/large-buffer values become more permissive.
- `--peak_finder_range_right E`: Right edge, in electrons, of the all-peaks histogram the peak finder runs on (default `1500`). Peak n sits near n e-, so this is roughly the maximum number of peaks found. Raised automatically if `resolution_at` needs more headroom. Not used by the gain-independent module (see `--all_peaks_n_max_electrons`).
- `--fit_range_right CHARGE [CHARGE ...] | auto`: Right charge bound (in electrons) for the parabolic nonlinearity fit. Accepts a single int applied to all extensions, one int per extension (e.g. `600 850 750 1050`), or the literal `auto` to pick it per extension with a data-driven estimator (see `--auto_fit_range_method`). `auto` can also stand in for individual extensions in a per-extension list (e.g. `auto 685 auto auto`): the estimator runs everywhere and the ints then override it, with those extensions recorded as `"method": "manual"` in `fit_range_estimate.json`.
- `--nonlinearity_correction_model`: `parabola` (default) or `data`. Selects the `NL(q)` used by `--correct_nonlinearity`. `data` uses the measured nonlinearity (peak charge minus peak number) at the integer electron value nearest each pixel's charge, clipped to the range of detected peaks, instead of the parabola. Assumes peak `n` is the `n`-electron peak, as the nonlinearity fit does.
- `--correct_nonlinearity`: Subtract the fitted nonlinearity from every pixel, `q_true = q - NL(q)` with `q` the measured pixel charge in e- and `NL` the extension's parabola, and save the corrected image (in electrons, input headers kept) in the run's hash output folder (alongside the plots, since the fit determines the model) as `NL_corrected_<file_name>.fits`. The frame corrected is the one the gain and pedestal were measured on, i.e. the pedestal-subtracted frame when `do_pedestal_subtraction` is on. Pixels above `fit_range_right` are corrected with the extrapolated parabola. Off by default.
- `--nonlinearity_fit_fix_intercept` / `--no-nonlinearity_fit_fix_intercept`: Hold the nonlinearity parabola's constant term `c` at 0 (default on). `Nonlin(0) = 0` identically — peak number 0 is the zero-electron peak, so its charge is the pedestal and `charge_minus_npeak[0] = 0` by construction — but a free `c` does not have to respect that. When `fit_range_right` reaches past the charge where the curve stops being quadratic, `c` absorbs the high-charge misfit and the drawn parabola ends up offset from a scatter that is pinned to the origin (measured at ~0.5 e- at 0 e- on a 1500-peak fit). Turning it off restores the free three-parameter fit. Applies to every parabolic nonlinearity fit: the production curve, the `auto` fit-range estimators, and the gain-independent module's gain scan and comparison curves. The returned coefficients keep their `(a, b, c)` shape with `c = 0`, so downstream evaluation is unchanged.
- `--nonlinearity_fit_known_points`: Known nonlinearity values to include in the parabola fit alongside the peaks up to `fit_range_right`, as `[charge, nonlinearity]` pairs in electrons (config: `[[1000, 3.2], [1500, 5.1]]`; CLI: flat `1000 3.2 1500 5.1`). Applied to every extension. In a JSON config it can instead be a list with one entry per extension, each `[[charge, nonlinearity], ...]` or `null` for none (e.g. `[[[1577.27, 4.1]], [[1577.27, 3.8]], null, [[1577.27, 4.4]]]`); the command line form is always the same for every extension. Default `null` (peaks only). They are drawn on the nonlinearity plots, one color per point (orange first), labeled `known NL(<charge>e-) = <NL value>` in the legend. A black dashed curve shows the same fit without the known points. Only the production curve fit uses them; the `auto` fit-range estimators and the gain-independent module do not.
- `--auto_fit_range_method {changepoint, var_a, noise_onset}`: Estimator used when `--fit_range_right auto` (default `changepoint`). The nonlinearity curve is a clean parabola up to some charge, then becomes noisy; the goal is to fit only the clean part.
  - `changepoint` (**default, recommended**): a two-stage, drift-free changepoint detector. It computes *local* roughness (the robust residual scatter of a local quadratic fit in a sliding window), which stays flat through the whole clean parabola regardless of its curvature and steps up sharply at the noise onset; it then refines the exact cut on the points that are guaranteed clean. Immune to the two `var_a` failure modes below.
  - `var_a`: picks the value minimizing `var(a)` (the variance of the parabola's curvature coefficient, `pcov[0,0]`). Works well when the score curve has a clear minimum, but can fail when the curve is near-linear (score decreases monotonically into the noisy tail) or flat (score collapses to the smallest range).
  - `noise_onset` (experimental): walks the curve forward and returns the first charge where a rolling median-absolute-deviation of residuals around a local quadratic exceeds `factor × baseline`, where the baseline is the median MAD in the early "safe" region (peaks within the bottom 30% of the above-`min_charge` range). Sensitivity is tuned via `--noise_onset_window` (sliding-window size in peaks; default `30`) and `--noise_onset_factor` (default `2.5`; lower = picks the onset earlier). Note: for globally noisy extensions the baseline MAD itself is elevated, so this estimator can floor well above the true visual onset — `changepoint` is usually preferable.
- `--changepoint_window N`: (changepoint) window size in peaks for the local-roughness quadratic fit / MAD (default `25`; odd values recommended).
- `--changepoint_factor FLOAT`: (changepoint) local roughness must exceed `factor × baseline` (or the absolute floor, whichever is larger) to mark the noise onset (default `4.0`; lower = more sensitive).
- `--changepoint_floor FLOAT`: (changepoint) absolute roughness floor in **electrons** (default `0.15`). Keeps an ultra-quiet early region from setting an impossibly tight threshold that would trip on sub-noise waviness. Sits between clean scatter (~0.02–0.05) and noisy-region roughness (~0.3–0.5); the first knob to revisit if a dataset has very different peak-location noise.
- `--changepoint_persist N`: (changepoint) number of consecutive points that must exceed the threshold to count as the noise onset (default `4`; higher = more robust to isolated outliers / precursor bumps).
- `--fit_range_confidence_tol FLOAT`: (changepoint) relative tolerance for the `var(a)` cross-check (default `0.15` = 15%). After estimating each extension, the independent `var(a)` value is computed; extensions where the two disagree by more than this fraction are flagged `LOW` confidence — printed at runtime with a `<-- REVIEW` marker and recorded in `fit_range_estimate.json` in the run's output folder (when `--save_output`).
- `--auto_fit_range_tolerance FLOAT` / `--auto_fit_range_max_charge_percentile FLOAT` / `--auto_fit_range_min_peaks_in_fit INT`: (`var_a` method only) optional guards — accept the smallest candidate within `(1+tol)×min_score`, cap candidates at a charge percentile, and reject candidates enclosing too few peaks, respectively.
- `--noise_onset_window N` / `--noise_onset_factor FLOAT`: (`noise_onset` method only) sliding-window size and the MAD-over-baseline factor for noise onset.

##### Single-electron resolution
The resolution step runs only when `--resolution_at` is given. For each requested charge *q* it takes a window `[q - n/2, q + n/2]` of the all-peaks histogram (n peaks wide), fits a constrained sum of Gaussians — a *comb* with means **fixed** at the already-detected peak locations, a single **shared** σ across all components, and free amplitudes — and scores how well that comb describes the window. The shared σ, in electrons, is the resolution figure of merit. Goodness of fit is reported as the comb's reduced χ² and as ΔAIC versus a single broad Gaussian "unresolved" null (positive favors the comb). Each window also gets a three-tier verdict (`well resolved` / `marginal` / `unresolved`). When the all-peaks histogram's right edge would fall short of the largest requested window, it is automatically extended to cover it.

- `-r`, `--resolution_at CHARGE [CHARGE ...]`: Charge value(s) in electrons to evaluate. Omitted (default `null`) → the resolution step is skipped entirely.
- `--resolution_window FLOAT`: Window width *n* in electrons (≈ number of peaks per window). Default `10`.
- `--resolution_sigma_well FLOAT`: σ (e-) below which peaks are **well resolved** (separation > 4σ). Default `0.25`.
- `--resolution_sigma_limit FLOAT`: σ (e-) at/above which peaks are **unresolved** (no central valley, separation < 2σ). Default `0.5`. Between `sigma_well` and `sigma_limit` the verdict is `marginal`.
- `--resolution_min_peak_frac FLOAT`: Detected/expected peak fraction below which the window is **unresolved** regardless of σ. Default `0.6`. This catches windows where the peak finder recovered far fewer than the ~`window+1` peaks the window should hold (e.g. 2 of 11) — too smeared to detect, so the under-determined σ is meaningless.

The verdict, σ, reduced χ², ΔAIC, and detected/expected peak counts are printed in a per-charge/per-extension table. With `--save_output` the same quantities are pivoted to one row per extension and merged into `extension_summary.csv` in the run's output folder, as `<field>_at_<charge>_e` columns (e.g. `sigma_e_at_500_e`) alongside gain, noise, and nonlinearity.

##### Plot layout (per plot type)
For each of `plot_all_peaks`, `plot_nonlinearity`, `plot_resolution`, `plot_charge_per_column`:
- `--<plot>_individual` / `--no-<plot>_individual`: Render one figure per extension
- `--<plot>_together` / `--no-<plot>_together`: Render a single 2×2 subplot grid for all extensions
- `--<plot>_individual_figsize W H`: Figure size for individual mode
- `--<plot>_subplots_figsize W H`: Figure size for the 2×2 together mode
- `--<plot>_sharex` / `--no-<plot>_sharex`: Share x-axis range across the 2×2 grid
- `--<plot>_sharey` / `--no-<plot>_sharey`: Share y-axis range across the 2×2 grid

`plot_zero_one` follows the same `_individual_figsize`/`_subplots_figsize`/`_sharex`/`_sharey` knobs, but instead of one shared pair of `_individual`/`_together` toggles it has four independent ones split by unit: `--plot_zero_one_ADU_individual`, `--plot_zero_one_ADU_together`, `--plot_zero_one_electrons_individual`, `--plot_zero_one_electrons_together` (see [Plotting](#plotting)) — so ADU vs. electron units and individual vs. combined layout can each be toggled on their own.

Top-row x-axis labels and right-column y-axis labels are always hidden in the 2×2 grids (tick labels are kept regardless of `sharex`/`sharey`).

For `plot_resolution` the two modes differ slightly because there can be several charges: `individual` renders one figure per (extension, charge), while `together` renders one 2×2 grid (one extension per panel) **per charge**. Its `_individual`/`_together` toggles default to `false`/`true`, and `--plot_resolution_subplots_figsize` defaults to `[13, 9]` (wider than the other plots to fit the wide windows). Resolution plots are only produced when `--resolution_at` is also set.

##### Plot styling
- `--plot_zero_one_yscale {linear, log}`
- `--plot_zero_one_xlim LOW HIGH` / `--plot_zero_one_ylim LOW HIGH` / `--plot_zero_one_ylim_electrons LOW HIGH`: limits for the zero/one peak plots (`ylim` for the ADU version, `ylim_electrons` for the electron-units version, which is separate because the electron peaks are taller). Omit or set `"default"` to use the fit range / auto range; `"none"` to autoscale.
- `--electron_fit_mode {transform, refit}`: how the electron-unit double-Gaussian curve is drawn on the zero/one plots — `transform` (default) analytically rescales the converged ADU fit; `refit` fits the double Gaussian directly to the electron-unit histogram.
- `--plot_all_peaks_xlim LEFT RIGHT` / `--plot_all_peaks_ylim BOTTOM TOP`
- `--plot_all_peaks_yscale {linear, log}`
- `--plot_nonlinearity_xlim ...` / `--plot_nonlinearity_ylim ...`: 2 values for a shared limit across extensions, or 8 values (`L1 R1 L2 R2 L3 R3 L4 R4`) for per-extension limits
- `--show_titles` / `--no-show_titles`: Toggle all `fig.suptitle` and `ax.set_title` calls across every plot

### Pedestal-subtraction caching

Pedestal subtraction is the slowest step in the pipeline. The first time it runs for a given source FITS file, the result is written to `<stem>.pedsub.fits` alongside the source (or to `--pedsub_cache_dir`). The pedestal-subtraction parameters (`axis`, `n_std_to_mask`, `use_biweight_loc`, `use_biweight_midvar`, `max_iter`) plus an algorithm-version tag are written into the FITS header as `PEDSUB_A`, `PEDSUB_N`, `PEDSUB_L`, `PEDSUB_V`, `PEDSUB_I`, and `PEDSUB_R`. The cache is reused only when all of them match the current run, so changing any parameter (or the algorithm) triggers a recompute.

On subsequent runs:
- If the cache exists AND the header params match the current params → cached arrays are loaded instantly.
- If params differ → cache is overwritten with a fresh computation.
- Pass `--force_pedsub` to bypass the cache without deleting it.

The cache is **not** automatically invalidated if the source FITS file itself changes (e.g. you re-stitch). In that case, delete the `.pedsub.fits` file or use `--force_pedsub`.

## Functions

### Analysis Functions

- `convert_to_electrons(data, pedestal, gain, flatten=True)`: Convert ADU values to electron counts
- `calculate_noise_gain(data, zero_one_test_range='auto', n_bins=75, fit_bounds='default', ..., data_stats=None)`: Determine noise and gain from charge data with data-driven zero/one peak initialization
- `pedestal_subtract(data, n_std_to_mask, axis='row', use_biweight_loc=True, use_biweight_midvar=True, max_iter=5, tol=0.01)`: Vectorized per-row/column pedestal subtraction for a single extension. Iteratively sigma-clips to the zero-peak core (recomputing both location and scale from the clipped pixels each pass) so overlapping one-electron pixels don't bias the pedestal; stops early once the median per-line shift falls below `tol`.
- `pedestal_subtract_ext_cached(data_ext, source_path, n_std_to_mask, axis='row', use_biweight_loc=True, use_biweight_midvar=True, cache_dir=None, force=False, verbose=True)`: Wraps `pedestal_subtract` across all extensions and caches the result to `<source_stem>.pedsub.fits`. Reuses the cache when the pedestal parameters and algorithm version recorded in the cached FITS header match.
- `get_fits(file_path)`: Load FITS file data (returns list of 4 extension arrays)
- `find_all_peaks(data, width, buffer, pedestal, noise, gain, ...)`: Single-extension peak finder
- `fit_nonlinearity(peaks, centers, pedestal, gain, fit_range_right, ..., fix_intercept=True, known_points=None)`: Single-extension parabolic fit of the nonlinearity curve. `fix_intercept` holds the constant term at 0 (see `--nonlinearity_fit_fix_intercept`); the returned coefficients are `(a, b, c)` either way, with `c = 0` and a zero-padded row/column in `pcov` when it is on.
- `get_zero_one_peaks_ext(data_ext, n_bins=75, fit_bounds='default', zero_one_test_range='auto', ..., data_stats=None)`: Per-extension zero/one peak fits; returns counts, edges, pedestals, gains, double-Gaussian popts, fit ranges. To refit at several `n_bins` (as the notebook's slider does), pass `data_stats=[pedestal_subtract.calibrate.full_data_stats(d) for d in data_ext]` so the binning-independent robust statistics are computed once.
- `get_all_peaks_ext(data_ext, widths, buffers, pedestals, double_gauss_popts, gains, ...)`: Per-extension peak finder
- `get_nonlinearity_ext(peaks_ext, centers_ext, pedestals, gains, fit_range_right_ext, ...)`: Per-extension nonlinearity fit
- `estimate_fit_range_right_changepoint_ext(peaks_ext, centers_ext, pedestals, gains, win=25, factor=4.0, floor=0.15, persist=4, cross_check=True, confidence_rel_tol=0.15, ...)`: Per-extension `fit_range_right` estimate via the default changepoint detector. Returns `(values, diagnostics)`, where `diagnostics` holds the chosen value, the `var(a)` cross-check, and a per-extension `OK`/`LOW` confidence label. (Single-extension core: `estimate_fit_range_right_changepoint`.)
- `estimate_optimal_fit_range_right_ext(...)` / `estimate_fit_range_right_by_noise_onset_ext(...)`: Alternative `fit_range_right` estimators (`var_a` and experimental `noise_onset`; see `--auto_fit_range_method`).
- `correct_nonlinearity_ext(data_ext, pedestals, gains, parabola_coeffs, fits_path, out_dir, data_nl=None)`: Apply `q_true = q - NL(q)` to every pixel of every extension and write `NL_corrected_<file_name>.fits` into `out_dir` (with `data_nl`, the per-extension measured charge-minus-peak-number, NL is looked up at the nearest integer e- instead of using the parabola); returns the corrected frames (e-) and the output path
- `get_nonlinearity_at_ext(q, parabola_coeffs, parabola_pcovs, fit_range_right_ext)`: Evaluate the parabolic fit at one or more charge values, per extension
- `resolution_at_charge(counts, centers, peaks, q, window, gain, s0)`: Single-extension single-electron resolution at charge `q`. Fits a constrained *n*-Gaussian comb (means fixed at the detected peaks, one shared σ, free amplitudes) in `[q - window/2, q + window/2]` and returns a result dict with `sigma_e`, `sigma_e_err`, `reduced_chi2`, `delta_aic`, `n_components`, `expected_peaks`, and the windowed data/fit curves for plotting.
- `resolution_at_charge_ext(counts_ext, centers_ext, peaks_ext, gains, double_gauss_popts, charges, window=10.0, sigma_well=0.25, sigma_limit=0.5, min_peak_frac=0.6, ...)`: Per-extension resolution at one or more charges. Returns a nested list (per extension, per charge) of result dicts, each augmented with `ext` and a `verdict`.
- `classify_resolution(res, sigma_well=0.25, sigma_limit=0.5, min_peak_frac=0.6)`: Map a result dict to a three-tier verdict (`well resolved` / `marginal` / `unresolved`) from the fitted σ relative to the 1 e- spacing and the detected/expected peak fraction.
- `summarize_resolution(results_ext, save_path=None)`: Print the per-charge/per-extension resolution table (σ, reduced χ², ΔAIC, peak counts, verdict); optionally writes that same long-format (one row per (ext, q)) table as CSV when `save_path` is given -- this is a diagnostic/console view, separate from the merged `extension_summary.csv` output. (`format_resolution_table(results_ext)` returns the printable table string.)
- `build_resolution_summary_wide(results_ext, charges)`: Pivot `resolution_at_charge_ext`'s per-(ext, charge) results to one row per extension, with a `<field>_at_<charge>_e` column per (field, charge) pair (`sigma_e`, `sigma_e_err`, `reduced_chi2`, `delta_aic`, `n_peaks_fit`, `expected_peaks`, `verdict`) -- mirrors `nonlinearity_at_<charge>_e`. Returns `(rows, charges)`.
- `summarize_extensions(gains, double_gauss_popts, parabola_coeffs, nonlinearity_charges=500, resolution_rows=None, resolution_charges=None, save_path=None)`: Print the per-extension summary table of gain (ADU/e-), noise (e-), nonlinearity at one or more charges (`nonlinearity_charges` accepts a scalar or list), and -- when `resolution_rows`/`resolution_charges` (from `build_resolution_summary_wide`) are given -- the merged resolution columns; writes it as one combined CSV (`extension_summary.csv`) when `save_path` is given. Returns the summary rows. (Builders: `build_extension_summary(..., resolution_rows=None)` returns `(rows, charges)` as dicts; `format_extension_summary(rows, charges, resolution_charges=None)` turns them into the printable table string.)

### Plotting Functions

Each plotting function has `plot_individual` and `plot_together` options to plot each extension separately or as one subplot, individual + 2x2 subplot sizes, `sharex`/`sharey` for the 2×2 subplot, `show_titles`, `show_plots`, and `save_plots`/`fig_path` for output.

- `plot_zero_one_peaks(data_ext, zero_one_counts_ext, zero_one_edges_ext, pedestals, gains, double_gauss_popts, zero_one_ranges, plot_adu_individual=False, plot_adu_together=True, plot_electrons_individual=False, plot_electrons_together=True, yscale='linear', ...)`: Visualize zero/one electron peak fits in ADU and/or electrons, individually and/or combined.
- `plot_all_peaks(counts_ext, peaks_ext, centers_ext, xlim, ylim='none', yscale='log', draw_lines=True, ...)`: Visualize the full charge distribution with a marker at each identified peak.
- `plot_nonlinearity(peaks_ext, parabola_coeffs, peak_charge_e_ext, charge_minus_npeak_ext, fit_range_right_ext, xlim='default', ylim='default', overlay=None, ...)`: Plot the nonlinearity curve and parabolic fit. `xlim`/`ylim` accept `'default'`, `'none'`, a single `(left, right)` applied to all extensions, or a list of 4 per-extension tuples. `overlay=(peak_charge_e_ext, charge_minus_npeak_ext, parabola_coeffs)` draws a second dataset (e.g. the nonlinearity-corrected image's) on the same axes, in `overlay_scatter_color`/`overlay_line_color` (purple/orange), with both datasets named by `overlay_names` in the legend.
- `set_plot_style(font_size, legend_font_size, use_tex=False, font_family=None)`: Apply the text style of the `font_size` / `legend_font_size` / `use_tex` / `font_family` config options to matplotlib's rcParams; both pipelines call it at the start of a run.
- `plot_images(panels, panel_labels=None, stack='horizontal', vmin=None, vmax=None, norm=LogNorm, norm_kwargs=None, cmap='gray', ...)`: One figure per extension showing 2-D pixel arrays on a greyscale with a shared colorbar. `panels` is a list of per-extension image lists; each extension's panels are drawn side by side (`stack='horizontal'`) or one above the other (`'vertical'`) on one shared norm. `norm` is a matplotlib normalization class: `LogNorm` (default), or e.g. `SymLogNorm` with `norm_kwargs=dict(linthresh=1)` to also show values <= 0, or `PowerNorm` with `norm_kwargs=dict(gamma=0.3)`. Unless `vmin`/`vmax` are given, the norm spans the 1st-99.9th percentile of the pixels (the positive ones for `LogNorm`). Individual-only: it has no `plot_together` option.
- `plot_resolution(results_ext, charges, sigma_well=0.25, sigma_limit=0.5, plot_individual=False, plot_together=True, ...)`: Plot the resolution windows from `resolution_at_charge_ext` — data histogram, the individual comb components, the comb sum, the single-Gaussian null, and a verdict annotation box. `plot_individual` makes one figure per (extension, charge); `plot_together` makes one 2×2 grid per charge.

### Utility

- `stitch_fits(file_path, **kwargs)`: Stitch multiple FITS files across each extension. All files must have the same shape.

## Structure

```
nonlinearity_studies/                # Repository root
├── config/                          # Example JSON configs
│   ├── default_nonlinearity.json
│   ├── resolution_study.json
│   └── ...
├── examples/                        # Examples directory
│   └── images/                      # Sample images
├── nonlinearity_studies/            # Main package directory
│   ├── __init__.py                  # Package initialization
│   ├── nonlinearity_studies.py      # Core analysis + plotting functions
│   ├── stitch_fits.py               # FITS image stitching utility
│   ├── notebook.py                  # Helpers for scripts/nonlinearity_analysis.ipynb
│   └── run_nonlinearity_studies.py  # Command-line interface
├── pyproject.toml                   # Build system configuration
├── setup.py                         # Package configuration
├── MANIFEST.in                      # Packaging manifest (non-code files to include)
├── environment.yml                  # Conda environment configuration file
├── .gitignore                       # File telling git which files not to track
└── README.md                        # This file
```

## License

This project is part of the Privitera research group at the University of Chicago for the DAMIC-M collaboration.

## Authors

- Abby Chriss
