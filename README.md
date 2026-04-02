# Nonlinearity Studies

A Python package for analyzing nonlinearity in CCDs.

## Overview

This package provides tools to:
- Load and analyze FITS format detector data
- Calculate noise and gain characteristics from zero and one electron peaks
- Measure nonlinearity across different charge ranges
- Generate publication-quality plots of charge distributions and nonlinearity curves
- Process data using command-line tools for batch analysis

## Features

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

#### Development mode

Install package in development mode (the -e flag makes files editable):

```bash
cd nonlinearity_studies
pip install -e .
```

#### OR regular installation

```bash
pip install .
```

### Requirements

- Python >= 3.8
- numpy >= 1.19.0
- matplotlib >= 3.3.0
- scipy >= 1.5.0
- astropy >= 4.0

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

**2. As a direct executable** (requires executable permissions - should be activated automatically but if not run ):
```bash
chmod +x nonlinearity_studies/run_nonlinearity_studies.py
./nonlinearity_studies/run_nonlinearity_studies.py [OPTIONS] <file_string>
```

**3. As a console script** (after pip installation):
```bash
run_nonlinearity_studies [OPTIONS] <file_string>
```

#### Options

- `--stitch_fits`: Stitch multi-extension FITS files before analysis
- `--plot_zero-one-peaks`: Plot zero and one electron peaks
- `--plot_all_peaks`: Plot all identified peaks
- `--get_nonlinearity_at CHARGES`: Calculate nonlinearity at specific charge values
- `--plot_nonlinearity`: Plot nonlinearity curve
- `--save_plots`: Save generated plots to disk

## Examples

First let's stitch 10 images together from examples/images/ten-images folder and run the analysis script on these images. Navigate to project directory in terminal and run:
```bash
./nonlinearity_studies/run_nonlinearity_studies.py \
    "examples/images/ten-images/*" \
    --stitch_fits \
    --plot_zero_one_peaks \
    --plot_all_peaks \
    --plot_nonlinearity \
    --save_plots
```
  
Now every time we want to analyze the stitched image again, we can pass the stitched image directly into the script instead of restitching and overwriting the stitched image, as follows:
```bash
./nonlinearity_studies/run_nonlinearity_studies.py \
    "combined-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits" \
    --plot_zero_one-peaks \
    --plot_nonlinearity \
    --save_plots
```
where the key change was removing the --stitch-fits flag.

If we want to just get the nonlinearity at specific charge values we can run
```bash
./nonlinearity_studies/run_nonlinearity_studies.py \
    "combined-fits/avg_img_CV_250x3500x500_bin1x1_125_10_stitched.fits" \
    --get_nonlinearity_at 10 50 500 100
```

## Core Functions

### Analysis Functions

- `convert_to_electrons(data, pedestal, gain, flatten=True)`: Convert ADU values to electron counts
- `calculate_noise_gain(data, zero_one_test_range=[8,15], n=200, fit_bounds='default')`: Determine noise and gain from charge data
- `get_fits(file_path)`: Load FITS file data
- `get_zero_one_peaks_ext(data, extension=None)`: Identify zero and one electron peaks
- `get_all_peaks_ext(data, extension=None)`: Identify all charge peaks
- `get_nonlinearity_ext(data, extension=None)`: Calculate nonlinearity across charge range
- `get_nonlinearity_at_ext(data, charge_values, extension=None)`: Calculate nonlinearity at specific charges

### Plotting Functions

- `plot_zero_one_peaks(data, **kwargs)`: Visualize zero and one electron peaks
- `plot_all_peaks(data, **kwargs)`: Visualize all identified peaks
- `plot_nonlinearity(data, **kwargs)`: Plot nonlinearity curve

## Module Structure

```
nonlinearity_studies/
├── __init__.py                  # Package initialization
├── nonlinearity_studies.py      # Core analysis functions
├── stitch_fits.py               # FITS image stitching utility
└── run_nonlinearity_studies.py  # Command-line interface
```

## Dependencies

This package requires:
- **numpy** >= 1.19.0
- **matplotlib** >= 3.3.0
- **scipy** >= 1.5.0
- **astropy** >= 4.0

The `stitch_fits` utility is included in the package for stitching multi-extension FITS files.

## Project Structure

The project structure includes:

```
nonlinearity_studies/
├── examples/                    # Example directory (images, scripts, etc)
├── nonlinearity_studies/        # Main package directory
├── setup.py                     # Package configuration
└── README.md                    # This file
```

## License

This project is part of the Privitera research group at the University of Chicago for the DAMIC-M collaboration. Please contact the authors for licensing information.

## Authors

- Abby Chriss