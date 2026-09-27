"""Statistics of flat T/Q/U patches, identical for MHD tiles and PySM patches.

    stats = patch_stats(T, Q, U, kmin=3, kmax=8, periodic=False)

Non-periodic patches are tapered (cosine edge over 1/8 of the side) before the
FFT; E/B, spectra and correlations then use python/analysis.py.
  aTT, aEE, aBB   power-law slopes of C_k over [kmin, kmax]
  bbee            BB/EE band-power ratio, same band
  rte, rtb, reb   correlation coefficients, same band
  p_med, p_q16/84 polarization fraction P/T quantiles (p0 applied by caller)
  S_med           median polarization-angle dispersion function at lag LAG px
                  (Planck 2018 XII definition, annulus of radius LAG)
  S_p_slope       d log S / d log p over pixels (Planck: about -1)
  skew_lnT, kurt_T, kurt_Q   one-point non-Gaussianity
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import qu_to_eb, spec2d_cross, fit_powerlaw

LAG = 2


def taper(n, frac=0.125):
    w = np.ones(n)
    m = int(frac * n)
    x = 0.5 * (1 - np.cos(np.pi * np.arange(m) / m))
    w[:m] = x; w[n - m:] = x[::-1]
    return np.outer(w, w)


def angle_dispersion(Q, U, lag=LAG):
    """S(x, lag): rms of psi(x) - psi(x + d) over the circle |d| = lag."""
    psi = 0.5 * np.arctan2(U, Q)
    acc, cnt = np.zeros_like(psi), 0
    for dx in range(-lag, lag + 1):
        for dy in range(-lag, lag + 1):
            if abs(np.hypot(dx, dy) - lag) > 0.5:
                continue
            d = psi - np.roll(np.roll(psi, dx, 0), dy, 1)
            d = 0.5 * np.arctan2(np.sin(2 * d), np.cos(2 * d))
            acc += d**2; cnt += 1
    return np.degrees(np.sqrt(acc / cnt))


def moments(x):
    x = x.ravel() - x.mean()
    s = x.std()
    return (x**3).mean() / s**3, (x**4).mean() / s**4 - 3


def patch_stats(T, Q, U, kmin=3, kmax=8, periodic=False, p0=None, E=None, B=None,
                beam_sigma=0.0):
    """E, B: E/B maps already decomposed on a larger domain (full sky, or the
    full periodic box face) and cut like T; avoids E->B leakage from the patch
    edges.  Without them E/B are computed on the (tapered) patch itself.
    beam_sigma: Gaussian beam sigma in units of the patch side, divided out of
    all spectra."""
    T = np.asarray(T, np.float64); Q = np.asarray(Q, np.float64); U = np.asarray(U, np.float64)
    n = T.shape[0]
    w = 1.0 if periodic else taper(n)
    Tw = (T - T.mean()) * w
    if E is None:
        E, B = qu_to_eb(Q * w, U * w)
    else:
        E = (np.asarray(E, np.float64) - np.mean(E)) * w
        B = (np.asarray(B, np.float64) - np.mean(B)) * w
    k, CTT, _ = spec2d_cross(Tw); _, CEE, _ = spec2d_cross(E); _, CBB, _ = spec2d_cross(B)
    _, CTE, _ = spec2d_cross(Tw, E); _, CTB, _ = spec2d_cross(Tw, B); _, CEB, _ = spec2d_cross(E, B)
    if beam_sigma:
        bl2 = np.exp(-(2 * np.pi * k * beam_sigma)**2)
        CTT, CEE, CBB = CTT / bl2, CEE / bl2, CBB / bl2
        CTE, CTB, CEB = CTE / bl2, CTB / bl2, CEB / bl2
    m = (k >= kmin) & (k <= kmax)
    r = lambda a, b, c: c[m].sum() / np.sqrt(a[m].sum() * b[m].sum())
    out = dict(aTT=fit_powerlaw(k, CTT, kmin, kmax)[1], aEE=fit_powerlaw(k, CEE, kmin, kmax)[1],
               aBB=fit_powerlaw(k, CBB, kmin, kmax)[1], bbee=CBB[m].sum() / CEE[m].sum(),
               rte=r(CTT, CEE, CTE), rtb=r(CTT, CBB, CTB), reb=r(CEE, CBB, CEB))
    P = np.hypot(Q, U)
    p = P / T if p0 is None else p0 * P / T
    ok = T > 0
    out["p_med"], out["p_q16"], out["p_q84"] = np.percentile(p[ok], [50, 16, 84])
    S = angle_dispersion(Q, U)
    inner = np.zeros_like(ok); c = n // 8
    inner[c:n - c, c:n - c] = True
    sel = ok & inner & (p > 0)
    out["S_med"] = float(np.median(S[sel]))
    lp, ls = np.log(p[sel]), np.log(np.maximum(S[sel], 1e-3))
    out["S_p_slope"] = float(np.polyfit(lp, ls, 1)[0])
    out["skew_lnT"] = moments(np.log(np.maximum(T[ok], T[ok].max() * 1e-6)))[0]
    out["kurt_T"] = moments(T)[1]
    out["kurt_Q"] = moments(Q)[1]
    return out


KEYS = ["aTT", "aEE", "aBB", "bbee", "rte", "rtb", "reb", "p_med", "S_med", "S_p_slope",
        "skew_lnT", "kurt_T", "kurt_Q"]
