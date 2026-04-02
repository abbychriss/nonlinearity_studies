# Package Structure Changes

## Summary
Your `nonlinearity_studies` project has been successfully converted into a pip-installable Python package!

## New Structure

```
nonlinearity_studies/
├── nonlinearity_studies/           ← The main package directory
│   ├── __init__.py                 ← Package initialization (NEW)
│   ├── nonlinearity_studies.py     ← Core analysis functions (moved from src/nonlinearity/)
│   ├── stitch_fits.py              ← FITS stitching utility (copied from tools/)
│   └── run_nonlinearity_studies.py ← CLI script (moved from src/)
├── extra/                          ← Old scripts and utilities (unchanged)
│   ├── plot_combined_charge_dist.py
│   ├── plot_sigma_xy.py
│   ├── vijay_linearity_study.py
│   └── old_nonlinearity_studies.py
├── src/                            ← OLD STRUCTURE (can be deleted)
│   ├── nonlinearity/
│   │   ├── __init__.py
│   │   └── nonlinearity_studies.py
│   └── run_nonlinearity_studies.py
├── setup.py                        ← Package configuration (NEW)
├── pyproject.toml                  ← Modern packaging format (NEW)
├── MANIFEST.in                     ← Include additional files (NEW)
├── .gitignore                      ← Git ignore rules (NEW)
├── README.md                       ← Comprehensive documentation (NEW)
└── .vscode/                        ← VS Code settings (unchanged)
```

## What Was Created/Modified

### New Files:
1. **`setup.py`** - Traditional package configuration. Includes:
   - Package metadata (name, version, author)
   - Dependencies (numpy, matplotlib, scipy, astropy)
   - Console script entry point for CLI

2. **`pyproject.toml`** - Modern Python packaging standard (PEP 517/518)
   - Alternative to setup.py (setuptools reads this)
   - Specifies build system requirements
   - Tool configurations

3. **`README.md`** - Comprehensive project documentation with:
   - Project overview and features
   - Installation instructions
   - Usage examples (as module and CLI)
   - API documentation
   - FAQ section

4. **`.gitignore`** - Git ignore patterns for:
   - Python cache and compiled files
   - Virtual environments
   - Distribution packages
   - IDE directories
   - Project-specific output files

5. **`MANIFEST.in`** - Include non-code files in packages
   - Ensures README.md is included
   - Includes Python files from package

6. **`nonlinearity_studies/__init__.py`** - Package initialization that:
   - Exposes main functions for easy importing (including `stitch_fits`)
   - Defines `__version__` and `__author__`
   - Provides the public API

### Modified Files:
1. **`nonlinearity_studies/run_nonlinearity_studies.py`** - Updated to:
   - Use relative imports (works with installed package)
   - Import bundled `stitch_fits` from local package
   - Fix argparse to use dashes in argument names
   - Add proper `if __name__ == '__main__'` block
   - Use `pathlib.Path` for all path operations
   - Improve documentation

### Copied Files:
- `src/nonlinearity/nonlinearity_studies.py` → `nonlinearity_studies/nonlinearity_studies.py`
- `src/run_nonlinearity_studies.py` → `nonlinearity_studies/run_nonlinearity_studies.py`
- `tools/stitch_fits.py` → `nonlinearity_studies/stitch_fits.py` (bundled for self-contained package)

## Installation

### Option 1: Development Installation (Recommended for Development)
```bash
cd /Users/abbychriss/Desktop/Privitera_335/scripts/nonlinearity_studies
pip install -e .
```
The `-e` flag installs in "editable" mode, so changes to the code are immediately reflected.

### Option 2: Regular Installation
```bash
cd /Users/abbychriss/Desktop/Privitera_335/scripts/nonlinearity_studies
pip install .
```

### Option 3: Install with Development Dependencies
```bash
pip install -e ".[dev]"
```
Installs testing and linting tools (pytest, black, flake8)

## Usage After Installation

### As a Python Module
```python
from nonlinearity_studies import (
    convert_to_electrons,
    calculate_noise_gain,
    get_fits,
    plot_nonlinearity,
    stitch_fits,  # Now bundled in the package!
)

# Now you can use these functions directly
data = get_fits("path/to/fits/file.fits")
noise, gain, pedestal = calculate_noise_gain(data)
```

### As a Command-Line Tool
```bash
# After installing the package:
run-nonlinearity-studies data.fits --plot-nonlinearity --save-plots

# Or run as a module:
python -m nonlinearity_studies.run_nonlinearity_studies data.fits --help
```

## Next Steps (Optional)

1. **Delete the old src/ directory** (no longer needed):
   ```bash
   rm -rf src/
   ```

2. **Initialize a git repository** and push to GitHub:
   ```bash
   git init
   git add .
   git commit -m "Initial commit: Convert to pip-installable package"
   ```

3. **Update the setup.py URL**:
   - Change `url=` to point to your actual GitHub repository
   - Update `author=` if needed

4. **Test the installation**:
   ```bash
   pip install -e .
   python -c "import nonlinearity_studies; print(nonlinearity_studies.__version__)"
   ```

## Dependencies

The package requires:
- **numpy** >= 1.19.0
- **matplotlib** >= 3.3.0
- **scipy** >= 1.5.0
- **astropy** >= 4.0

**All dependencies are included!** The `stitch_fits` utility is now bundled in the package, so no external tools dependencies are needed.

## Notes

- The old `src/` directory structure is still present. You can safely delete it after verifying everything works.
- All imports have been updated to work with the installed package.
- The `stitch_fits` utility is now bundled in the package, making it completely self-contained with no external tool dependencies.
- All code refactored to use `pathlib.Path` for clean path handling instead of string concatenation.
- The package now follows Python packaging best practices (PEP 427, PEP 517, PEP 518).
- The console script entry point `run-nonlinearity-studies` will be automatically installed when you install the package.

## Verification

To verify everything is working:

```bash
# Install in development mode
cd /Users/abbychriss/Desktop/Privitera_335/scripts/nonlinearity_studies
pip install -e .

# Test imports
python -c "from nonlinearity_studies import get_fits, plot_nonlinearity, stitch_fits; print('Success!')"

# Check version
python -c "import nonlinearity_studies; print(f'Version: {nonlinearity_studies.__version__}')"

# View help
python -m nonlinearity_studies.run_nonlinearity_studies --help
```

Questions or issues? Check the README.md for more detailed documentation!
