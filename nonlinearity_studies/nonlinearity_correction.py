"""
nonlinearity_correction — apply the fitted nonlinearity model to an image.
"""
import re
from pathlib import Path

import numpy as np

from astropy.io import fits
from pedestal_subtract.core import convert_to_electrons

from .fit_data import parabola

NL_CORRECTED_KEY = 'NLCORR'  # primary-header flag: image is already in electrons, gain == 1


def is_nl_corrected(fits_path):
    """True if the file was written by correct_nonlinearity_ext (NLCORR flag)."""
    return bool(fits.getheader(fits_path, 0).get(NL_CORRECTED_KEY, False))


def correct_nonlinearity_ext(data_ext, pedestals, gains, parabola_coeffs, fits_path, out_dir, data_nl=None):
    """Subtract the modeled nonlinearity from every pixel: q_true = q - NL(q).

    ``data_ext`` must be the frames the pedestals/gains were measured on (ADU);
    the corrected frames are in electrons. Saved in ``out_dir`` (the run-hash
    output folder, since the fit determines the model) as
    NL_corrected_<file_name>.fits, keeping the input's headers.
    """
    fits_path = Path(fits_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / (f"NL_corrected_{re.sub(r'[.]fits([.]fz)?$', '', fits_path.name)}.fits")
    corrected_ext = []
    with fits.open(fits_path) as hdul:
        hdul[0].header[NL_CORRECTED_KEY] = (True, 'image is NL-corrected, in electrons')
        for ext, data in enumerate(data_ext):
            q = convert_to_electrons(data, pedestals[ext], gains[ext], flatten=False)
            if data_nl is None:
                nl = parabola(q, *parabola_coeffs[ext])
            else:
                nl_n = np.asarray(data_nl[ext], dtype=float)
                nl = nl_n[np.clip(np.rint(q), 0, len(nl_n) - 1).astype(int)]
            corrected_ext.append(q - nl)
            hdul[ext + 1].data = corrected_ext[ext]
        hdul.writeto(out_path, overwrite=True)
    return corrected_ext, out_path
