"""
summary — split out of nonlinearity_studies.py.
"""
import csv
import numpy as np
from pathlib import Path

_RESOLUTION_WIDE_FIELDS = [
    ('sigma_e', 'sigma_e', 'Res sigma [e-]'),
    ('sigma_e_err', 'sigma_e_err', 'Res sigma err [e-]'),
    ('reduced_chi2', 'reduced_chi2', 'Res red chi2'),
    ('delta_aic', 'delta_aic', 'Res dAIC'),
    ('n_peaks_fit', 'n_components', 'Res peaks fit'),
    ('expected_peaks', 'expected_peaks', 'Res exp peaks'),
    ('verdict', 'verdict', 'Res verdict'),
]


_RESOLUTION_TERSE_FIELDS = {'sigma_e', 'verdict'}


def _normalize_charges(nonlinearity_charges):
    """Return nonlinearity_charges as a flat list of floats (accepts a scalar or a list)."""
    if nonlinearity_charges is None:
        return []
    if isinstance(nonlinearity_charges, (list, tuple, np.ndarray)):
        return [float(q) for q in nonlinearity_charges]
    return [float(nonlinearity_charges)]


def _nonlin_column_name(charge):
    """CSV/column name for the nonlinearity evaluated at a charge, e.g. 'nonlinearity_at_500_e'."""
    return f'nonlinearity_at_{charge:g}_e'


def _resolution_column_name(field, charge):
    """CSV/column name for a resolution quantity at a charge, e.g. 'sigma_e_at_500_e'."""
    return f'{field}_at_{charge:g}_e'


def _resolution_wide_fieldnames(charges):
    return [_resolution_column_name(field, q) for q in charges for field, _, _ in _RESOLUTION_WIDE_FIELDS]


def build_resolution_summary_wide(results_ext, charges):
    """Pivot per-(ext, charge) resolution results to one row per extension.

    Mirrors build_extension_summary's nonlinearity_at_<charge>_e columns: instead of one row per
    (ext, q) (the long format previously written by summarize_resolution), this produces one row
    per extension with a '<field>_at_<charge>_e' column per (field, charge) pair, so the result
    can be merged into the extension summary table by 'ext'.

    Args:
        results_ext: nested list from resolution_at_charge_ext -- one list per extension of one
            result dict (or None, if that window couldn't be fit) per charge, in the same order
            as `charges`.
        charges: the charges (in e-) results_ext is ordered by (as passed to
            resolution_at_charge_ext).

    Returns:
        (rows, charges) where charges is the normalized list of charges, and rows is a list of
        dicts (one per extension), keyed by 'ext' plus one '<field>_at_<charge>_e' key per
        (field, charge) pair present. A charge whose fit failed for a given extension (result is
        None) simply has no keys for that charge in that row -- same as a blank cell in the CSV.
    """
    charges = _normalize_charges(charges)
    rows = []
    for ext, charge_results in enumerate(results_ext):
        row = {'ext': ext + 1}
        for q, res in zip(charges, charge_results):
            if res is None:
                continue
            for col_field, res_key, _ in _RESOLUTION_WIDE_FIELDS:
                row[_resolution_column_name(col_field, q)] = res[res_key]
        rows.append(row)
    return rows, charges


def build_extension_summary(gains, double_gauss_popts, parabola_coeffs, nonlinearity_charges=500,
                            resolution_rows=None):
    """Build a per-extension summary of gain, noise in electrons, nonlinearity at charge(s), and
    (optionally) single-electron resolution at charge(s).

    Args:
        gains: list of per-extension gains (ADU/e-), i.e. m1 - m0 from the double-Gaussian fit.
        double_gauss_popts: list of per-extension double-Gaussian popts (s0, m0, s1, m1, N0, N1).
            s0 (popt[0]) is the standard deviation of the zero-electron peak in ADU.
        parabola_coeffs: list of per-extension parabola coefficients (a, b, c) from the
            nonlinearity fit.
        nonlinearity_charges: charge or list of charges (in e-) at which to evaluate the
            nonlinearity (default 500).
        resolution_rows: optional list of per-extension resolution dicts (one per extension, in
            ext order) from build_resolution_summary_wide. When given, their
            '<field>_at_<charge>_e' columns are merged into each row by 'ext'. Pass None (the
            default) to omit resolution columns entirely.

    Returns:
        (rows, charges) where charges is the normalized list of nonlinearity charges, and rows is
        a list of dicts (one per extension) with flat keys: 'ext', 'gain_adu_per_e', 'noise_e',
        one 'nonlinearity_at_<charge>_e' per nonlinearity charge, and (if resolution_rows was
        given) the merged resolution columns.
    """
    charges = _normalize_charges(nonlinearity_charges)
    resolution_by_ext = {r['ext']: r for r in resolution_rows} if resolution_rows else {}
    rows = []
    for ext, gain in enumerate(gains):
        noise_adu = double_gauss_popts[ext][0]  # s0: std of zero-electron peak in ADU
        noise_e = noise_adu / gain  # convert ADU -> e- by dividing by gain (ADU/e-)
        a, b, c = parabola_coeffs[ext]
        row = {
            'ext': ext + 1,
            'gain_adu_per_e': gain,
            'noise_e': noise_e,
            'noise_adu': noise_adu,
        }
        for q in charges:
            row[_nonlin_column_name(q)] = a * q**2 + b * q + c
        res_row = resolution_by_ext.get(ext + 1)
        if res_row:
            row.update({k: v for k, v in res_row.items() if k != 'ext'})
        rows.append(row)
    return rows, charges


def _summary_fieldnames(charges, resolution_charges=None):
    fields = ['ext', 'gain_adu_per_e', 'noise_e', 'noise_adu'] + [_nonlin_column_name(q) for q in charges]
    if resolution_charges:
        fields += _resolution_wide_fieldnames(resolution_charges)
    return fields


def format_extension_summary(rows, charges, resolution_charges=None, verbose=False):
    """Format per-extension summary rows (from build_extension_summary) as an aligned table
    string. Pass the resolution_charges from build_resolution_summary_wide to also show the
    merged single-electron-resolution columns; a row missing a value (e.g. a charge whose fit
    failed for that extension) renders as a blank cell rather than raising.

    By default (verbose=False) only the sigma_e and verdict resolution columns are shown, to
    keep the console table readable; pass verbose=True to also show sigma_e_err, reduced_chi2,
    delta_aic, n_peaks_fit, and expected_peaks per charge.
    """
    headers = ['EXT', 'Gain [ADU/e-]', 'Noise [e-]', 'Noise [ADU]'] + [f'Nonlin @ {q:g} e- [e-]' for q in charges]
    fields = ['ext', 'gain_adu_per_e', 'noise_e', 'noise_adu'] + [_nonlin_column_name(q) for q in charges]
    resolution_charges = _normalize_charges(resolution_charges)
    resolution_field_specs = _RESOLUTION_WIDE_FIELDS if verbose else [
        spec for spec in _RESOLUTION_WIDE_FIELDS if spec[0] in _RESOLUTION_TERSE_FIELDS
    ]
    for q in resolution_charges:
        headers += [f'{label} @ {q:g} e-' for _, _, label in resolution_field_specs]
        fields += [_resolution_column_name(col_field, q) for col_field, _, _ in resolution_field_specs]

    def _cell(row, field):
        if field not in row:
            return ''
        val = row[field]
        if field == 'ext':
            return f'{val}'
        if isinstance(val, str):
            return val
        return f'{val:.4f}'

    cells = [[_cell(r, f) for f in fields] for r in rows]
    widths = [max(len(h), *(len(row[i]) for row in cells)) if cells else len(h)
              for i, h in enumerate(headers)]
    fmt = lambda vals: '  '.join(v.rjust(widths[i]) for i, v in enumerate(vals))
    lines = [fmt(headers), '-' * (sum(widths) + 2 * (len(widths) - 1))]
    lines += [fmt(row) for row in cells]
    return '\n'.join(lines)


def summarize_extensions(gains, double_gauss_popts, parabola_coeffs, nonlinearity_charges=500,
                         resolution_rows=None, resolution_charges=None, verbose=False, save_path=None):
    """Print a per-extension table of gain, noise in e-, nonlinearity, and (optionally)
    single-electron resolution, and optionally save it as CSV.

    resolution_rows/resolution_charges are the outputs of build_resolution_summary_wide; when
    given, their columns are merged into each row by 'ext' so one combined CSV covers both
    nonlinearity and resolution, keyed by extension. The resolution columns are already printed
    on their own by summarize_resolution, so the console table here only shows gain, noise, and
    nonlinearity (verbose has no effect on this table); the saved CSV always has every column.

    Returns the list of summary rows from build_extension_summary.
    """
    rows, charges = build_extension_summary(gains, double_gauss_popts, parabola_coeffs,
                                             nonlinearity_charges, resolution_rows=resolution_rows)
    print('\n***********************************************************')
    print('\nNonlinearity analysis summary:\n')
    print(format_extension_summary(rows, charges, resolution_charges=None, verbose=verbose))

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = _summary_fieldnames(charges, resolution_charges)
        with save_path.open('w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, restval='')
            writer.writeheader()
            writer.writerows(rows)
        print(f'\nSaved per-extension summary table (CSV) to {save_path}')

    return rows
