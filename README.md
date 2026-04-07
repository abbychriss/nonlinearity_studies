# Nonlinearity Studies

A Python package for analyzing nonlinearity in CCDs.

## Overview

- **Noise & Gain Calculation**: Determines pedestal, noise and gain from data and converts to e-
- **Peak Finder**: Finds all electron peaks in charge distribution
- **Nonlinearity Computation**: Quantifies detector response linearity as a function of charge
- **Visualization**: Creates histograms of charge distribution and plots nonlinearity curve
- **Image Stitching**: Combine multi-extension FITS images and run analysis on stitched image

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

## Usage

### As a Python module

```python
from nonlinearity_studies import (
    convert_to_electrons,
    calculate_noise_gain,
    get_fits,
    plot_nonlinearity,
)

# Load FITS data
data = get_fits("path/to/fits/file.fits")

# Calculate noise and gain
noise, gain, pedestal = calculate_noise_gain(data)

# Convert to electrons
data_electrons = convert_to_electrons(data, pedestal, gain)

# Plot results
plot_nonlinearity(data_electrons)
```

### Command-line interface

There are three ways to run the analysis:

**1. As a Python module** (works during development):
```bash
python -m nonlinearity_studies.run_nonlinearity_studies [OPTIONS] <file_string>
```

**2. As a direct executable** (requires executable permissions - should be activated automatically but if not run `chmod +x nonlinearity_studies/run_nonlinearity_studies.py`):

```bash
./nonlinearity_studies/run_nonlinearity_studies.py [OPTIONS] <file_string>
```

**3. As a console script** (after pip installation):

```bash
run-nonlinearity-studies [OPTIONS] <file_string>
```

#### Options

- `-f`, `--stitch_fits`: Stitch multi-extension FITS files before analysis
- `-z`, `--plot_zero_one`: Plot zero and one electron peaks
- `-a`, `--plot_all_peaks`: Plot all identified peaks
- `-g`, `--get_nonlinearity_at CHARGES`: Calculate nonlinearity at specific charge values
- `-n`, `--plot_nonlinearity`: Plot nonlinearity curve
- `-s`, `--save_plots`: Save generated plots to local computer
- `-o`, `--output_dir`: (Optional) output directory for saved plots
- `-v`, `--verbose`: Print verbose output
- `--nimages N`: Number of stitched images (used for plot labeling)

## Examples

For the first example, we will fit the zero and one electron peaks for a single 250x3500 image with 1x1 binning and 500 skips. Navigate to project directory in terminal and run:

```bash
run-nonlinearity-studies \
    "examples/images/ten-images/avg_img_CV_250x3500x500_bin1x1_125_20260317_213403_0.fz" \
    --plot_zero_one \
```

Next, let's stitch 10 images together from examples/images/ten-images folder and run the analysis script on these
 images:

```bash
run-nonlinearity-studies \
    "examples/images/ten-images/*" \
    --stitch_fits \
    --plot_zero_one \
    --plot_nonlinearity \
```
  
Now every time we want to analyze the stitched image again, we can pass the stitched image directly into the script and remove the --stitch_fits flag, instead of restitching and overwriting the stitched image. Run the stitched image and save the plots:

```bash
run-nonlinearity-studies \
    "examples/images/combined-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits" \
    --plot_zero_one \
    --plot_nonlinearity \
    --save_plots
```

If we want to just get the nonlinearity at specific charge values we can run

```bash
run-nonlinearity-studies \
    "examples/images/combined-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits" \
    --get_nonlinearity_at 10 50 500 1000
```

## Functions

### Analysis Functions

- `convert_to_electrons(data, pedestal, gain, flatten=True)`: Convert ADU values to electron counts
- `calculate_noise_gain(data, zero_one_test_range=[8,15], n=200, fit_bounds='default')`: Determine noise and gain from charge data (single file)
- `get_fits(file_path)`: Load FITS file data
- `get_zero_one_peaks_ext(data, extension=None)`: Identify zero and one electron peaks for all extensions
- `get_all_peaks_ext(data, extension=None)`: Find all electron peaks for everu extension
- `get_nonlinearity_ext(data, extension=None)`: Calculate nonlinearity curve for all extensions
- `get_nonlinearity_at_ext(data, charge_values, extension=None)`: Calculate nonlinearity at specific charges for all extensions

### Plotting Functions

- `plot_zero_one_peaks(data, **kwargs)`: Visualize zero and one electron peaks
- `plot_all_peaks(data, **kwargs)`: Visualize all identified peaks
- `plot_nonlinearity(data, **kwargs)`: Plot nonlinearity curve

### Utility

- `stitch_fits`(file_path, **kwargs): stitching multiple FITS files across each extension. Data be same shape.

## Structure

```
nonlinearity_studies/
├── __init__.py                      # Package initialization
├── examples/                        # Examples directory (images, scripts, etc)
├── nonlinearity_studies/            # Main package directory
    ├── nonlinearity_studies.py      # Core analysis functions
    ├── stitch_fits.py               # FITS image stitching utility
    └── run_nonlinearity_studies.py  # Command-line interface
├── setup.py                         # Package configuration
├── environment.yml                  # Conda environment configuration file
├── .gitignore                       # File telling git which files not to track
└── README.md                        # This file
```


## License

This project is part of the Privitera research group at the University of Chicago for the DAMIC-M collaboration. Please contact the authors for licensing information.

## Authors

- Abby Chriss
