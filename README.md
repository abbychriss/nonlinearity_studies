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

- **Noise & Gain Calculation**: Determines detector noise and gain from charge data
- **Peak Analysis**: Identifies and gets properties of zero and one electron peaks
- **Nonlinearity Measurement**: Quantifies detector response linearity vs. charge
- **Visualization**: Creates histograms and plots of charge distributions
- **Batch Processing**: Command-line interface for processing multiple datasets
- **Image Stitching**: Combine multi-extension FITS images

## Installation

### Clone the repository

Navigate to directory you want to store project in.
```bash
mkdir nonlinearity_studies
git clone https://github.com/abbychriss/nonlinearity_studies.git
cd nonlinearity_studies
```

### From source (development mode)

Clone the repository and install in development mode:

```bash
cd /path/to/nonlinearity_studies
pip install -e .
```

### Regular installation

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

**2. As a direct executable** (requires executable permissions):
```bash
chmod +x nonlinearity_studies/run_nonlinearity_studies.py
./nonlinearity_studies/run_nonlinearity_studies.py [OPTIONS] <file_string>
```

**3. As a console script** (after pip installation):
```bash
run-nonlinearity-studies [OPTIONS] <file_string>
```

#### Options

- `--stitch-fits`: Stitch multi-extension FITS files before analysis
- `--plot-zero-one-peaks`: Plot zero and one electron peaks
- `--plot-all-peaks`: Plot all identified peaks
- `--get-nonlinearity-at CHARGE`: Calculate nonlinearity at specific charge values
- `--plot-nonlinearity`: Plot nonlinearity curve
- `--save-plots`: Save generated plots to disk

#### Example

```bash
python -m nonlinearity_studies.run_nonlinearity_studies \
    --stitch-fits \
    --plot-nonlinearity \
    --save-plots \
    "my_dataset.fits"
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

The original project structure includes:

```
nonlinearity_studies/
├── nonlinearity_studies/        # Main package directory
├── extra/                       # Extra analysis scripts and older implementations
├── setup.py                     # Package configuration
└── README.md                    # This file
```

## Examples

See the `extra/` directory for additional example scripts and analysis utilities, including:

- `plot_sigma_xy.py`: Plot size of clusters vs energy to visualize nonlinearity in front side and back side Fe-55 events (work in progress)
- `vijay_linearity_study.py`: Alternative linearity analysis approach

## License

This project is part of the Privitera research group at the University of Chicago for the DAMIC-M collaboration. Please contact the authors for licensing information.

## Authors

- Abby Chriss

## References

For information on detector linearity studies and gain/noise calculations, see the extra analysis scripts in the `extra/` directory.

## FAQ

**Q: What FITS format does this support?**
A: The package supports standard FITS files with multi-extension capability. See `astropy.io.fits` documentation for supported formats.

**Q: Can I process multiple files?**
A: Yes, use the command-line interface with the `--stitch-fits` option for batch processing.

**Q: Is the stitch_fits utility included?**
A: Yes! The `stitch_fits` module is included in the package for stitching multi-extension FITS files. You can use it directly via `from nonlinearity_studies import stitch_fits` or via the CLI with `--stitch-fits`.

## Support

For issues, questions, or contributions, please contact the package author or check the project repository.
