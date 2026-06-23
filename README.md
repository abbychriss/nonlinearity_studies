# Nonlinearity Studies

A Python package for analyzing nonlinearity in CCDs.

## Overview

- **Image Stitching**: Combine multi-extension FITS images row-wise and run analysis on stitched image
- **Pedestal Subtraction**: Independently computes and subtracts pedestal along specified axis (using my pedestal_subtract package)
- **Noise & Gain Calculation**: Determines noise and gain from data and converts to e-
- **Peak Finder**: Finds all electron peaks in charge distribution
- **Nonlinearity Computation**: Quantifies nonlinearity as a function of charge by fitting nonlinearity curve to parabola
- **Single-Electron Resolution**: Quantifies single electron resolution at a given charge *q* by fitting *n* Gaussians in window around *q* and reporting the shared peak width σ (e-), goodness of fit, and a resolved/marginal/unresolved verdict
- **Visualization**: Plots zero-one electron peaks in ADU and e-, histograms of charge distribution at any q with peak labels, and nonlinearity curves. Can be plotted per extension or together as 2×2 subplot. Control figure sizes, axis sharing, scales, limits, titles, etc. using config.
- **Meta analysis**: Plot results of nonlinearity analysis (noise, gain, resolution, nonlinearity) as function of image parameter using summary csv tables produced during individual analyses

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
cd nonlinearity_studies
pip install -e .
```

####  regular installation

```bash
pip install .
```

### Requirements for direct pip installation

- Python >= 3.8
- numpy >= 1.26.0
- matplotlib >= 3.8.0
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
counts_ext, edges_ext, peaks_ext, centers_ext, hist_ranges = get_all_peaks_ext(
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

#### Options

##### Pipeline / I/O
- `-j`, `--json PATH`: Load command-line arguments from a JSON config file (CLI arguments override JSON values)
- `-f`, `--stitch_fits`: Stitch multi-extension FITS files before analysis
- `-s`, `--save_plots`: Save generated plots as JPEGs
- `--save_csv`: Save the per-extension summary table as `extension_summary.csv` in the output directory (independent of `--save_plots`)
- `-o`, `--output_dir DIR`: Directory for all saved output (plots, summary CSV, config snapshot). Defaults to a `plots/` folder alongside the source FITS.
- `--show_plots` / `--no-show_plots`: Display plots interactively via `plt.show()` (default `true`). When disabled, figures are simply not shown (and are closed to free memory) instead of switching matplotlib to a non-interactive backend — useful for headless/batch runs, and still works alongside `--save_plots`.
- `-v`, `--verbose`: Print verbose output
- `--nimages N`: Number of stitched images (used for plot labeling). Auto-detected from filenames matching `_N_stitched`
- `--extra_plot_title TITLE`: Additional title text prepended to every plot title

##### What to plot
Each plot type is controlled by its `_individual` and `_together` toggles (see [Plot layout](#plot-layout-per-plot-type)) — a plot is generated whenever either is `true`.
- `-z`, `--plot_zero_one_adu`: Render the zero/one electron peak fits in ADU (units selector for the zero/one plot)
- `--plot_zero_one_electrons`: Render the zero/one electron peak fits in e- (independent of the ADU flag)
- `-g`, `--get_nonlinearity_at CHARGE...`: Evaluate the nonlinearity polynomial at one or more charge values
- `-r`, `--resolution_at CHARGE...`: Quantify single-electron resolution at one or more charge values (see [Single-electron resolution](#single-electron-resolution))

##### Per-extension summary table
Every run prints a per-extension summary table to stdout (no flag required) with the **gain** in ADU/e- (the separation between the zero- and one-electron peaks from the double-Gaussian fit), the **noise** in e- (the standard deviation of the zero-electron peak, converted from ADU to electrons by dividing by the gain), and one **nonlinearity** column per charge passed to `--get_nonlinearity_at` (falling back to 500 e- when none is given). With `--save_csv`, the same table is written as a CSV (`extension_summary.csv`) in the output directory — one row per extension with a header row, ready to load with `numpy.genfromtxt(..., delimiter=',', names=True)` or `numpy.loadtxt(..., delimiter=',', skiprows=1)`.

##### Pedestal subtraction
- `--do_pedestal_subtraction` / `--no-do_pedestal_subtraction`: Toggle pedestal subtraction (default `true`)
- `--n_std_to_mask FLOAT`: Mask threshold (in standard deviations) for the pedestal estimator
- `--pedsub_max_iter INT`: Max sigma-clip iterations for pedestal subtraction (default `5`). Each pass re-estimates the pedestal from the clipped zero-peak core and stops early once it converges.
- `--pedestal_subtraction_axis {row, col, row_then_col, col_then_row}`: Axis to compute the pedestal across
- `--use_biweight_loc` / `--use_biweight_midvar`: Use Tukey biweight location/midvariance instead of mean/std
- `--pedsub_cache_dir DIR`: Directory for the pedestal-subtracted FITS cache (default: alongside the source FITS)
- `--force_pedsub`: Recompute pedestal subtraction even if a matching cache exists

##### Fitting
- `--peak_finder_widths W [W ...]`: Minimum peak width required by `scipy.signal.find_peaks`, in **electrons** (internally multiplied by `bin_factor` to convert to bins). Scalar or one per extension. Larger = filters narrow noise spikes more aggressively.
- `--peak_finder_buffers B [B ...]`: Buffer (in **bins**) SUBTRACTED from `bin_factor` to compute the minimum neighbor-peak distance: `d = bin_factor - buffer`. With the default `bin_factor=10`: `buffer=0` → 10 bins = 1 electron spacing (physical), `buffer=3` → 7 bins ≈ 0.7 electron (loose), `buffer=-2` → 12 bins = 1.2 electron (strict). Larger buffer = looser; smaller/negative = stricter.
- `--peak_finder_prominences P [P ...]`: Minimum peak prominence in **histogram counts** (same units as the y-axis of `plot_all_peaks`). Scalar, one per extension, or `null` to disable. Often the most robust filter — measures how far a peak sticks up above its surrounding baseline, regardless of width.
- `--bin_factor N`: Number of histogram bins per electron in the all-peaks histogram (default `10`). Also drives the `peak_finder_widths` electron-to-bin conversion (`width_in_bins = width_in_electrons * bin_factor`) and the buffer math (`distance_in_bins = bin_factor - buffer`). Higher = finer histogram, but small-width/large-buffer values become more permissive.
- `--fit_range_right CHARGE [CHARGE ...] | auto`: Right charge bound (in electrons) for the parabolic nonlinearity fit. Accepts a single int applied to all extensions, one int per extension (e.g. `600 850 750 1050`), or the literal `auto` to pick it per extension with a data-driven estimator (see `--auto_fit_range_method`).
- `--auto_fit_range_method {changepoint, var_a, noise_onset}`: Estimator used when `--fit_range_right auto` (default `changepoint`). The nonlinearity curve is a clean parabola up to some charge, then becomes noisy; the goal is to fit only the clean part.
  - `changepoint` (**default, recommended**): a two-stage, drift-free changepoint detector. It computes *local* roughness (the robust residual scatter of a local quadratic fit in a sliding window), which stays flat through the whole clean parabola regardless of its curvature and steps up sharply at the noise onset; it then refines the exact cut on the points that are guaranteed clean. Immune to the two `var_a` failure modes below.
  - `var_a`: picks the value minimizing `var(a)` (the variance of the parabola's curvature coefficient, `pcov[0,0]`). Works well when the score curve has a clear minimum, but can fail when the curve is near-linear (score decreases monotonically into the noisy tail) or flat (score collapses to the smallest range).
  - `noise_onset` (experimental): walks the curve forward and returns the first charge where a rolling median-absolute-deviation of residuals around a local quadratic exceeds `factor × baseline`, where the baseline is the median MAD in the early "safe" region (peaks within the bottom 30% of the above-`min_charge` range). Sensitivity is tuned via `--noise_onset_window` (sliding-window size in peaks; default `30`) and `--noise_onset_factor` (default `2.5`; lower = picks the onset earlier). Note: for globally noisy extensions the baseline MAD itself is elevated, so this estimator can floor well above the true visual onset — `changepoint` is usually preferable.
- `--changepoint_window N`: (changepoint) window size in peaks for the local-roughness quadratic fit / MAD (default `25`; odd values recommended).
- `--changepoint_factor FLOAT`: (changepoint) local roughness must exceed `factor × baseline` (or the absolute floor, whichever is larger) to mark the noise onset (default `4.0`; lower = more sensitive).
- `--changepoint_floor FLOAT`: (changepoint) absolute roughness floor in **electrons** (default `0.15`). Keeps an ultra-quiet early region from setting an impossibly tight threshold that would trip on sub-noise waviness. Sits between clean scatter (~0.02–0.05) and noisy-region roughness (~0.3–0.5); the first knob to revisit if a dataset has very different peak-location noise.
- `--changepoint_persist N`: (changepoint) number of consecutive points that must exceed the threshold to count as the noise onset (default `4`; higher = more robust to isolated outliers / precursor bumps).
- `--fit_range_confidence_tol FLOAT`: (changepoint) relative tolerance for the `var(a)` cross-check (default `0.15` = 15%). After estimating each extension, the independent `var(a)` value is computed; extensions where the two disagree by more than this fraction are flagged `LOW` confidence — printed at runtime with a `<-- REVIEW` marker and recorded in `fit_range_estimate.json` in the run's output folder (when `--save_plots`).
- `--auto_fit_range_tolerance FLOAT` / `--auto_fit_range_max_charge_percentile FLOAT` / `--auto_fit_range_min_peaks_in_fit INT`: (`var_a` method only) optional guards — accept the smallest candidate within `(1+tol)×min_score`, cap candidates at a charge percentile, and reject candidates enclosing too few peaks, respectively.
- `--noise_onset_window N` / `--noise_onset_factor FLOAT`: (`noise_onset` method only) sliding-window size and the MAD-over-baseline factor for noise onset.

##### Single-electron resolution
The resolution step runs only when `--resolution_at` is given. For each requested charge *q* it takes a window `[q - n/2, q + n/2]` of the all-peaks histogram (n peaks wide), fits a constrained sum of Gaussians — a *comb* with means **fixed** at the already-detected peak locations, a single **shared** σ across all components, and free amplitudes — and scores how well that comb describes the window. The shared σ, in electrons, is the resolution figure of merit. Goodness of fit is reported as the comb's reduced χ² and as ΔAIC versus a single broad Gaussian "unresolved" null (positive favors the comb). Each window also gets a three-tier verdict (`well resolved` / `marginal` / `unresolved`). When the all-peaks histogram's right edge would fall short of the largest requested window, it is automatically extended to cover it.

- `-r`, `--resolution_at CHARGE [CHARGE ...]`: Charge value(s) in electrons to evaluate. Omitted (default `null`) → the resolution step is skipped entirely.
- `--resolution_window FLOAT`: Window width *n* in electrons (≈ number of peaks per window). Default `10`.
- `--resolution_sigma_well FLOAT`: σ (e-) below which peaks are **well resolved** (separation > 4σ). Default `0.25`.
- `--resolution_sigma_limit FLOAT`: σ (e-) at/above which peaks are **unresolved** (no central valley, separation < 2σ). Default `0.5`. Between `sigma_well` and `sigma_limit` the verdict is `marginal`.
- `--resolution_min_peak_frac FLOAT`: Detected/expected peak fraction below which the window is **unresolved** regardless of σ. Default `0.6`. This catches windows where the peak finder recovered far fewer than the ~`window+1` peaks the window should hold (e.g. 2 of 11) — too smeared to detect, so the under-determined σ is meaningless.

The verdict, σ, reduced χ², ΔAIC, and detected/expected peak counts are printed in a per-charge/per-extension table. With `--save_csv` the same table is written to `resolution_summary.csv` in the run's output folder.

##### Plot layout (per plot type)
For each of `plot_zero_one`, `plot_all_peaks`, `plot_nonlinearity`, `plot_resolution`:
- `--<plot>_individual` / `--no-<plot>_individual`: Render one figure per extension
- `--<plot>_together` / `--no-<plot>_together`: Render a single 2×2 subplot grid for all extensions
- `--<plot>_individual_figsize W H`: Figure size for individual mode
- `--<plot>_subplots_figsize W H`: Figure size for the 2×2 together mode
- `--<plot>_sharex` / `--no-<plot>_sharex`: Share x-axis range across the 2×2 grid
- `--<plot>_sharey` / `--no-<plot>_sharey`: Share y-axis range across the 2×2 grid

Top-row x-axis labels and right-column y-axis labels are always hidden in the 2×2 grids (tick labels are kept regardless of `sharex`/`sharey`).

For `plot_resolution` the two modes differ slightly because there can be several charges: `individual` renders one figure per (extension, charge), while `together` renders one 2×2 grid (one extension per panel) **per charge**. Its `_individual`/`_together` toggles default to `false`/`true`, and `--plot_resolution_subplots_figsize` defaults to `[13, 9]` (wider than the other plots to fit the wide windows). Resolution plots are only produced when `--resolution_at` is also set.

##### Plot styling
- `--plot_zero_one_yscale {linear, log}`
- `--plot_all_peaks_xlim LEFT RIGHT` / `--plot_all_peaks_ylim BOTTOM TOP`
- `--plot_all_peaks_yscale {linear, log}`
- `--plot_nonlinearity_xlim ...` / `--plot_nonlinearity_ylim ...`: 2 values for a shared limit across extensions, or 8 values (`L1 R1 L2 R2 L3 R3 L4 R4`) for per-extension limits
- `--show_titles` / `--no-show_titles`: Toggle all `fig.suptitle` and `ax.set_title` calls across every plot

## Examples

For the first example, we will fit the zero and one electron peaks for a single 250x3500 image with 1x1 binning and 500 skips. Navigate to project directory in terminal and run:

```bash
run-nonlinearity-studies \
    "examples/images/ten-images/avg_img_CV_250x3500x500_bin1x1_125_20260317_213403_0.fz" \
    --plot_zero_one_adu --plot_zero_one_together
```

Next, let's stitch 10 images together from `examples/images/ten-images/` and run the analysis script on the stitched output:

```bash
run-nonlinearity-studies \
    "examples/images/ten-images/*" \
    --stitch_fits \
    --plot_zero_one_adu --plot_zero_one_together \
    --plot_nonlinearity_together
```

Now every time we want to analyze the stitched image again, we can pass the stitched image directly into the script (and drop `--stitch_fits`). The pedestal-subtracted result is cached next to the source as `<stem>.pedsub.fits` and reused on subsequent runs if the pedsub params match. Run the stitched image and save the plots:

```bash
run-nonlinearity-studies \
    "examples/images/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits" \
    --plot_zero_one_adu --plot_zero_one_electrons --plot_zero_one_together \
    --plot_nonlinearity_together \
    --save_plots
```

If we want only the nonlinearity at specific charge values:

```bash
run-nonlinearity-studies \
    "examples/images/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits" \
    --get_nonlinearity_at 10 50 500 1000
```

To quantify single-electron resolution at one or more charges, pass `--resolution_at`. This prints (and, with `--save_csv`, saves) a per-charge/per-extension table of σ, reduced χ², ΔAIC, peak counts, and verdict, and with a plot flag draws the windowed comb fits:

```bash
run-nonlinearity-studies \
    "examples/images/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits" \
    --resolution_at 200 600 1000 --resolution_window 10 \
    --plot_resolution_together --save_plots --save_csv
```

A ready-made config is provided at `config/resolution_study.json`.

To let the pipeline choose the parabola fit range per extension, pass `--fit_range_right auto`. This uses the `changepoint` estimator and cross-checks each result against `var(a)` to assign a per-extension confidence label:

```bash
run-nonlinearity-studies \
    "examples/images/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits" \
    --fit_range_right auto --plot_nonlinearity_together --save_plots
```

With `--save_plots`, the full per-extension estimate (chosen value, cross-check, confidence) is written to `fit_range_estimate.json` in the run's output folder. Add `--verbose` to also print the confidence label per extension and a warning for any that need review:

```
EXT 3: changepoint fit_range_right = 693 e-  (var(a) cross-check = 100 e-, rel diff = 86%, confidence = LOW)  <-- REVIEW
  WARNING: low-confidence fit_range_right on EXT [3] (changepoint and var(a) disagree by > 15%); inspect the nonlinearity plot(s).
```

The same options can be written in JSON using argparse destination names. A config might look like:

```json
{
  "file_string": "examples/images/stitched-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits",
  "stitch_fits": false,
  "save_plots": true,
  "save_csv": true,
  "output_dir": null,
  "show_plots": true,

  "plot_zero_one_adu": true,
  "plot_zero_one_electrons": true,
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

  "plot_zero_one_individual": false,
  "plot_zero_one_together": true,
  "plot_zero_one_individual_figsize": [8, 6],
  "plot_zero_one_subplots_figsize": [12, 8],
  "plot_zero_one_yscale": "linear",
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
run-nonlinearity-studies -j config/default_nonlinearity.json --no-save_plots
```

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
- `calculate_noise_gain(data, zero_one_test_range='auto', n=200, fit_bounds='default')`: Determine noise and gain from charge data with data-driven zero/one peak initialization
- `pedestal_subtract(data, n_std_to_mask, axis='row', use_biweight_loc=True, use_biweight_midvar=True, max_iter=5, tol=0.01)`: Vectorized per-row/column pedestal subtraction for a single extension. Iteratively sigma-clips to the zero-peak core (recomputing both location and scale from the clipped pixels each pass) so overlapping one-electron pixels don't bias the pedestal; stops early once the median per-line shift falls below `tol`.
- `pedestal_subtract_ext_cached(data_ext, source_path, n_std_to_mask, axis='row', use_biweight_loc=True, use_biweight_midvar=True, cache_dir=None, force=False, verbose=True)`: Wraps `pedestal_subtract` across all extensions and caches the result to `<source_stem>.pedsub.fits`. Reuses the cache when the pedestal parameters and algorithm version recorded in the cached FITS header match.
- `get_fits(file_path)`: Load FITS file data (returns list of 4 extension arrays)
- `find_all_peaks(data, width, buffer, pedestal, noise, gain, ...)`: Single-extension peak finder
- `fit_nonlinearity(peaks, centers, pedestal, gain, fit_range_right, ...)`: Single-extension parabolic fit of the nonlinearity curve
- `get_zero_one_peaks_ext(data_ext, n=200, fit_bounds='default', zero_one_test_range='auto')`: Per-extension zero/one peak fits; returns counts, edges, pedestals, gains, double-Gaussian popts, fit ranges
- `get_all_peaks_ext(data_ext, widths, buffers, pedestals, double_gauss_popts, gains, ...)`: Per-extension peak finder
- `get_nonlinearity_ext(peaks_ext, centers_ext, pedestals, gains, fit_range_right_ext, ...)`: Per-extension nonlinearity fit
- `estimate_fit_range_right_changepoint_ext(peaks_ext, centers_ext, pedestals, gains, win=25, factor=4.0, floor=0.15, persist=4, cross_check=True, confidence_rel_tol=0.15, ...)`: Per-extension `fit_range_right` estimate via the default changepoint detector. Returns `(values, diagnostics)`, where `diagnostics` holds the chosen value, the `var(a)` cross-check, and a per-extension `OK`/`LOW` confidence label. (Single-extension core: `estimate_fit_range_right_changepoint`.)
- `estimate_optimal_fit_range_right_ext(...)` / `estimate_fit_range_right_by_noise_onset_ext(...)`: Alternative `fit_range_right` estimators (`var_a` and experimental `noise_onset`; see `--auto_fit_range_method`).
- `get_nonlinearity_at_ext(q, parabola_coeffs, parabola_pcovs, fit_range_right_ext)`: Evaluate the parabolic fit at one or more charge values, per extension
- `resolution_at_charge(counts, centers, peaks, q, window, gain, s0)`: Single-extension single-electron resolution at charge `q`. Fits a constrained *n*-Gaussian comb (means fixed at the detected peaks, one shared σ, free amplitudes) in `[q - window/2, q + window/2]` and returns a result dict with `sigma_e`, `sigma_e_err`, `reduced_chi2`, `delta_aic`, `n_components`, `expected_peaks`, and the windowed data/fit curves for plotting.
- `resolution_at_charge_ext(counts_ext, centers_ext, peaks_ext, gains, double_gauss_popts, charges, window=10.0, sigma_well=0.25, sigma_limit=0.5, min_peak_frac=0.6, ...)`: Per-extension resolution at one or more charges. Returns a nested list (per extension, per charge) of result dicts, each augmented with `ext` and a `verdict`.
- `classify_resolution(res, sigma_well=0.25, sigma_limit=0.5, min_peak_frac=0.6)`: Map a result dict to a three-tier verdict (`well resolved` / `marginal` / `unresolved`) from the fitted σ relative to the 1 e- spacing and the detected/expected peak fraction.
- `summarize_resolution(results_ext, save_path=None)`: Print the per-charge/per-extension resolution table (σ, reduced χ², ΔAIC, peak counts, verdict); writes it as CSV (`resolution_summary.csv`) when `save_path` is given. (`format_resolution_table(results_ext)` returns the printable table string.)
- `summarize_extensions(gains, double_gauss_popts, parabola_coeffs, nonlinearity_charges=500, save_path=None)`: Print the per-extension summary table of gain (ADU/e-), noise (e-), and nonlinearity at one or more charges (`nonlinearity_charges` accepts a scalar or list); writes it as CSV when `save_path` is given. Returns the summary rows. (Builders: `build_extension_summary` returns `(rows, charges)` as dicts; `format_extension_summary(rows, charges)` turns them into the printable table string.)

### Plotting Functions

Each plotting function accepts `plot_individual` and `plot_together` toggles, individual + 2×2 figure sizes, `sharex`/`sharey` for the 2×2 grid, `show_titles`, `show_plots` (display interactively, or close the figure when `False`), and `save_plots`/`fig_path` for output.

- `plot_zero_one_peaks(data_ext, zero_one_counts_ext, zero_one_edges_ext, pedestals, gains, double_gauss_popts, zero_one_ranges, do_plot_adu=True, do_convert_to_electrons=False, yscale='linear', ...)`: Visualize zero/one electron peak fits in ADU and/or electrons.
- `plot_all_peaks(counts_ext, peaks_ext, centers_ext, xlim, ylim='none', yscale='log', draw_lines=True, ...)`: Visualize the full charge distribution with a marker at each identified peak.
- `plot_nonlinearity(peaks_ext, parabola_coeffs, peak_charge_e_ext, charge_minus_npeak_ext, fit_range_right_ext, xlim='default', ylim='default', ...)`: Plot the nonlinearity curve and parabolic fit. `xlim`/`ylim` accept `'default'`, `'none'`, a single `(left, right)` applied to all extensions, or a list of 4 per-extension tuples.
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
│   └── images/                      # Sample FITS/FZ images
├── nonlinearity_studies/            # Main package directory
│   ├── __init__.py                  # Package initialization
│   ├── nonlinearity_studies.py      # Core analysis + plotting functions
│   ├── stitch_fits.py               # FITS image stitching utility
│   └── run_nonlinearity_studies.py  # Command-line interface
├── pyproject.toml                   # Build system configuration
├── setup.py                         # Package configuration
├── MANIFEST.in                      # Packaging manifest (non-code files to include)
├── environment.yml                  # Conda environment configuration file
├── .gitignore                       # File telling git which files not to track
└── README.md                        # This file
```

Generated output (pedestal-subtracted FITS caches, the `plots/` directory, and `*.jpeg` figures) is git-ignored and not shown above.


## License

This project is part of the Privitera research group at the University of Chicago for the DAMIC-M collaboration. Please contact the authors for licensing information.

## Authors

- Abby Chriss
