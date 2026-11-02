"""
Peak localization and numbering for the gain-independent (ADU) spectrum.

``sub_bin_log_parabola`` and ``gaussian_refit_failed_peaks`` refine each
scipy-found peak (a bin index) to a sub-bin location. ``find_zero_peak_ext``
then turns that list of locations into the charge-vs-peak-number arrays
``q(n)`` that every calculation in this package consumes.
"""

import numpy as np
from scipy.optimize import curve_fit

from ..fit_data import _gauss_single


# ---------------------------------------------------------------------------
# Sub-bin peak localization
# ---------------------------------------------------------------------------
def sub_bin_log_parabola(counts, centers, peak_idx, half_width=3, min_good_bins=5):
    """Sub-bin peak positions by weighted log-parabola, batched over all peaks.

    For a Gaussian sampled on a histogram, log(counts) is a parabola in the bin
    center, so the peak position is the vertex of a parabola fit to log(counts)
    over a small window. That is a 3x3 linear solve per peak, so all peaks are
    done in a single stacked np.linalg.solve -- no Python loop, no curve_fit.

    Weighting: counts are Poisson, so var(ln y) ~= var(y)/y^2 = 1/y, and the
    correct weight is w = y. (This is the Poisson case; Guo's original y^2
    weighting assumes constant variance on y, which is not our situation.)
    Unweighted log-least-squares would badly over-weight the near-empty flank
    bins.

    Returns a dict of arrays, all length len(peak_idx):
        mu        sub-bin peak position, same units as `centers`
        sigma     fitted Gaussian width
        amp       fitted peak height
        n_counts  raw counts summed over the fit window
        sigma_mu  naive position error, sigma/sqrt(n_counts)   [PRELIMINARY]
        ok        bool, False where the fit was not a usable maximum
        refit     bool, True where gaussian_refit_failed_peaks supplied the fit
                  instead (always False here; set by that function)
        fit_lo    left edge of the window the fit was made over, in `centers`
        fit_hi    right edge of that window
    plus the scalar `bin_width` of the histogram the fit was made on, which a
    caller needs to rescale `amp` onto a differently binned histogram.
    """
    counts = np.asarray(counts, dtype=float)
    centers = np.asarray(centers, dtype=float)
    peak_idx = np.asarray(peak_idx, dtype=int)
    n_pk = peak_idx.size
    if n_pk == 0:
        empty = np.zeros(0)
        return dict(mu=empty, sigma=empty, amp=empty, n_counts=empty,
                    sigma_mu=empty, ok=np.zeros(0, dtype=bool),
                    refit=np.zeros(0, dtype=bool), fit_lo=empty, fit_hi=empty,
                    bin_width=float(centers[1] - centers[0]) if centers.size > 1 else np.nan)

    bin_width = centers[1] - centers[0]
    offs = np.arange(-half_width, half_width + 1)
    idx = peak_idx[:, None] + offs[None, :]                 # (P, W)
    inside = (idx >= 0) & (idx < counts.size)
    idx_clipped = np.clip(idx, 0, counts.size - 1)

    y = counts[idx_clipped]
    # Local abscissa, measured from the peak's own bin center. Keeps the 3x3
    # normal matrix well conditioned regardless of where in the spectrum we are.
    x = centers[idx_clipped] - centers[peak_idx][:, None]

    good = inside & (y > 0)
    w = np.where(good, y, 0.0)
    lny = np.log(np.where(good, y, 1.0))

    X = np.stack([np.ones_like(x), x, x ** 2], axis=-1)     # (P, W, 3)
    M = np.einsum('pwi,pw,pwj->pij', X, w, X)               # (P, 3, 3)
    rhs = np.einsum('pwi,pw,pw->pi', X, w, lny)             # (P, 3)

    n_good = good.sum(axis=1)
    # Solve only the well-posed systems; np.linalg.solve raises on any singular
    # matrix in the stack, so substitute the identity for the bad ones and mask
    # the results afterwards.
    solvable = n_good >= min_good_bins
    det = np.linalg.det(M)
    solvable &= np.isfinite(det) & (np.abs(det) > 0)
    M_safe = np.where(solvable[:, None, None], M, np.eye(3)[None, :, :])
    rhs_safe = np.where(solvable[:, None], rhs, 0.0)
    # numpy >= 2 resolves solve(a, b) with b.ndim == a.ndim - 1 as a stack of
    # matrices rather than a stack of vectors, so make the vector axis explicit.
    coef = np.linalg.solve(M_safe, rhs_safe[..., None])[..., 0]   # (P,3) -> [c,b,a]

    c, b, a = coef[:, 0], coef[:, 1], coef[:, 2]
    with np.errstate(divide='ignore', invalid='ignore'):
        mu_local = -b / (2.0 * a)
        sigma = np.sqrt(-1.0 / (2.0 * a))
        amp = np.exp(c - b ** 2 / (4.0 * a))

    ok = (
        solvable
        & (a < 0)                                   # must be concave: a maximum
        & np.isfinite(mu_local) & np.isfinite(sigma) & np.isfinite(amp)
        # vertex must land inside the fit window, else the window is not
        # centered on the feature and the quadratic is extrapolating
        & (np.abs(mu_local) <= half_width * abs(bin_width))
        & (sigma > 0)
    )

    mu = np.where(ok, centers[peak_idx] + mu_local, centers[peak_idx])
    sigma = np.where(ok, sigma, np.nan)
    amp = np.where(ok, amp, np.nan)

    n_counts = np.where(good, y, 0.0).sum(axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        sigma_mu = np.where(n_counts > 0, sigma / np.sqrt(n_counts), np.nan)

    # The window each fit was made over, in the same units as `centers`, so a
    # plot can draw the fitted curve over exactly the range that produced it.
    half_span = half_width * abs(bin_width)
    fit_lo = np.clip(centers[peak_idx] - half_span, centers[0], centers[-1])
    fit_hi = np.clip(centers[peak_idx] + half_span, centers[0], centers[-1])

    return dict(mu=mu, sigma=sigma, amp=amp, n_counts=n_counts,
                sigma_mu=sigma_mu, ok=ok, refit=np.zeros(n_pk, dtype=bool),
                fit_lo=fit_lo, fit_hi=fit_hi, bin_width=float(bin_width))


def gaussian_refit_failed_peaks(counts, centers, peak_idx, sub_bin_result, sigma_seed,
                                 k_sigma=3.0, min_good_bins=5):
    """Fallback localization for peaks the log-parabola solve couldn't fit.

    ``sub_bin_log_parabola`` marks a peak ``ok=False`` when its local log(counts)
    window isn't concave, is singular, or its vertex lands outside the window --
    all symptoms of a window too narrow to see the actual peak shape. Rather than
    leave those few peaks (typically a few percent of the total) at the raw bin
    center, this refits a single Gaussian directly to counts vs. centers via
    scipy.optimize.curve_fit, seeded with the peak finder's raw peak location and
    ``sigma_seed`` (the shared zero/one-electron sigma), over a window sized to
    that seed (+/- k_sigma * sigma_seed) rather than a fixed bin count -- wide
    enough to constrain the fit regardless of sub_bin_peakfinder_bin_factor.
    Only the failed subset is refit, so a per-peak Python loop (curve_fit isn't
    batchable) stays cheap.

    Modifies ``sub_bin_result``'s ``mu``, ``sigma``, ``amp``, ``n_counts``,
    ``sigma_mu``, ``ok``, ``refit``, ``fit_lo`` and ``fit_hi`` in place for every
    failed peak the Gaussian refit succeeds on -- ``refit`` marking which peaks
    this function, rather than the log-parabola, supplied the fit for. Peaks where the refit also fails (too few nonzero bins, no
    convergence, or a nonsensical fit) are left as ``sub_bin_log_parabola`` set
    them: raw bin center, ``ok=False``.

    Returns ``sub_bin_result`` (same dict, mutated).
    """
    counts = np.asarray(counts, dtype=float)
    centers = np.asarray(centers, dtype=float)
    peak_idx = np.asarray(peak_idx, dtype=int)
    bin_width = centers[1] - centers[0]
    half_width_bins = max(1, int(np.ceil(k_sigma * sigma_seed / abs(bin_width))))

    for i in np.flatnonzero(~sub_bin_result['ok']):
        p = peak_idx[i]
        lo = max(0, p - half_width_bins)
        hi = min(counts.size, p + half_width_bins + 1)
        x = centers[lo:hi]
        y = counts[lo:hi]
        good = y > 0
        if good.sum() < min_good_bins:
            continue

        mu0 = centers[p]
        amp0 = counts[p] if counts[p] > 0 else y[good].max()
        try:
            popt, _ = curve_fit(
                _gauss_single, x[good], y[good], p0=(amp0, mu0, sigma_seed),
                sigma=np.sqrt(y[good]), absolute_sigma=True,
                bounds=([0.0, x.min(), 0.0], [np.inf, x.max(), np.inf]),
                maxfev=2000,
            )
        except RuntimeError:
            continue

        amp_fit, mu_fit, sigma_fit = popt
        if not (x.min() <= mu_fit <= x.max()) or sigma_fit <= 0:
            continue

        n_counts_fit = float(y[good].sum())
        sub_bin_result['mu'][i] = mu_fit
        sub_bin_result['sigma'][i] = sigma_fit
        sub_bin_result['amp'][i] = amp_fit
        sub_bin_result['n_counts'][i] = n_counts_fit
        sub_bin_result['sigma_mu'][i] = sigma_fit / np.sqrt(n_counts_fit) if n_counts_fit > 0 else np.nan
        sub_bin_result['ok'][i] = True
        sub_bin_result['refit'][i] = True
        sub_bin_result['fit_lo'][i] = x[0]
        sub_bin_result['fit_hi'][i] = x[-1]

    return sub_bin_result


# ---------------------------------------------------------------------------
# Peak numbering
# ---------------------------------------------------------------------------
def find_zero_peak_ext(counts_ext, peaks_ext, mu_ext):
    """Number the detected ADU peaks from the zero-electron peak outwards.

    ``find_all_peaks`` returns every local maximum from the histogram's left
    edge onwards, so its index 0 is not guaranteed to be the zero-electron peak:
    a low-count bump in the pedestal's left tail can clear the peak finder's
    width/distance filters and take index 0, shifting every real peak number by
    one. (Measured on real data: such a bump carried a few hundred counts
    against ~180,000 in the true zero peak.)

    The zero-electron peak is by construction the tallest feature in the
    spectrum -- by orders of magnitude -- so it is identified here as the
    detected peak with the most counts, exactly the criterion
    ``pedestal_subtract._auto_zero_one_setup`` uses to seed its own zero peak
    (``np.argmax`` of the smoothed counts). Everything to its left is dropped.

    Args:
        counts_ext: list (one per extension) of histogram counts, from
            ``get_all_peaks_ext``.
        peaks_ext: list (one per extension) of peak bin indices into those counts.
        mu_ext: list (one per extension) of sub-bin peak positions in ADU,
            aligned 1:1 with ``peaks_ext[ext]``.

    Returns:
        (q_ext, zero_idx_ext) where ``q_ext[ext]`` is the ADU charge of peak
        number n = 0, 1, 2, ... (so ``q_ext[ext][0]`` is the zero-electron
        peak), and ``zero_idx_ext[ext]`` is the index the zero peak had in the
        original ``peaks_ext[ext]`` -- i.e. how many peaks were dropped.
    """
    q_ext = []
    zero_idx_ext = []

    for counts, peaks, mu in zip(counts_ext, peaks_ext, mu_ext):
        peaks = np.asarray(peaks, dtype=int)
        mu = np.asarray(mu, dtype=float)
        if peaks.size == 0:
            q_ext.append(np.zeros(0))
            zero_idx_ext.append(0)
            continue

        zero_idx = int(np.argmax(np.asarray(counts, dtype=float)[peaks]))
        q_ext.append(mu[zero_idx:])
        zero_idx_ext.append(zero_idx)

    return q_ext, zero_idx_ext


def peak_number_mask(n, n_max=None):
    """Boolean mask selecting the peak numbers a fit is allowed to use.

    ``n_max`` counts peaks *including* the zero-electron peak (the convention
    used by --gain_indep_fit_n_min/--gain_indep_fit_n_max), so n_max=1000 keeps
    peak numbers 0..999. ``None`` keeps everything.
    """
    n = np.asarray(n)
    if n_max is None:
        return np.ones(n.shape, dtype=bool)
    return n <= int(n_max) - 1


def n_max_per_ext(n_max, n_ext):
    """Broadcast ``n_max`` to one value per extension.

    A scalar (or None) applies to every extension; a list is checked for length
    and passed through, so the per-extension changepoint estimate from
    ``spectrum_fit.estimate_fit_n_max_ext`` can be threaded through every fit
    unchanged.
    """
    if n_max is None or np.isscalar(n_max):
        return [n_max] * n_ext
    n_max = list(n_max)
    if len(n_max) != n_ext:
        raise ValueError(f'n_max has {len(n_max)} entries but there are {n_ext} extensions')
    return n_max
