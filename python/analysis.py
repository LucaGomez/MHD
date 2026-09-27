"""Analysis pipeline reproducing Stalpes, Collins & Huffenberger (2024).

Given a simulation snapshot (rho, v, H on a periodic cube) this computes
  * 3d shell-averaged power spectra of rho, v, H
  * projected maps T, Q, U along a line of sight  (Eqs. 4, 7, 8 of the paper)
  * E/B decomposition in the flat-sky approximation (Eq. 9)
  * 2d annulus-averaged auto- and cross-spectra of T, E, B
  * power-law fits  C_k = A k^alpha  over the inertial window
"""
import numpy as np, os, json, re

# scipy.fft is ~2x faster than numpy.fft and, crucially, keeps single-precision
# input in single precision (numpy always promotes to complex128).  At 512^3
# that is the difference between 2 GB and 1 GB per transform.
try:
    from scipy import fft as _scipy_fft
    _NW = int(os.environ.get("MHD_FFT_WORKERS", "-1"))   # -1 = all cores

    class _fft:                                        # thin threaded wrapper
        fftn = staticmethod(lambda a: _scipy_fft.fftn(a, workers=_NW))
        fft2 = staticmethod(lambda a: _scipy_fft.fft2(a, workers=_NW))
        ifft2 = staticmethod(lambda a: _scipy_fft.ifft2(a, workers=_NW))
except ImportError:                                   # pragma: no cover
    _fft = np.fft


# ----------------------------------------------------------------- I/O
def read_snapshot(path):
    meta = {}
    with open(path.replace(".bin", ".txt")) as f:
        for line in f:
            p = line.split()
            if len(p) >= 2:
                try:    meta[p[0]] = float(p[1])
                except ValueError: meta[p[0]] = " ".join(p[1:])
    n = int(meta["n"])
    a = np.fromfile(path, dtype=np.float32).reshape(7, n, n, n)
    # keep the fields in the precision they were written in; reductions that
    # need it accumulate in float64 explicitly
    d = dict(rho=a[0], v=a[1:4], H=a[4:7], meta=meta, n=n)
    return d


# ------------------------------------------------------- power spectra
def _kgrid(n, dim):
    k1 = np.fft.fftfreq(n) * n
    if dim == 3:
        kx, ky, kz = np.meshgrid(k1, k1, k1, indexing="ij")
        return np.sqrt(kx**2 + ky**2 + kz**2), (kx, ky, kz)
    kx, ky = np.meshgrid(k1, k1, indexing="ij")
    return np.sqrt(kx**2 + ky**2), (kx, ky)


_KCACHE = {}
def kmag(n, dim):
    if (n, dim) not in _KCACHE:
        _KCACHE[(n, dim)] = _kgrid(n, dim)
    return _KCACHE[(n, dim)]


def shell_average(P2, n, dim):
    """Average power in shells of unit thickness in |k| (cosmology convention)."""
    kk, _ = kmag(n, dim)
    idx = np.rint(kk).astype(np.int64).ravel()
    w = P2.ravel()
    nb = n // 2 + 1
    keep = idx < nb
    cnt = np.bincount(idx[keep], minlength=nb)
    tot = np.bincount(idx[keep], weights=w[keep], minlength=nb)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = tot / cnt
    k = np.arange(nb)
    return k[1:], out[1:], cnt[1:]


def spec3d(fields):
    """Shell-averaged 3d spectrum; `fields` is a list of components summed."""
    n = fields[0].shape[0]
    P2 = None
    for f in fields:
        F = _fft.fftn(f)
        p = (F.real**2 + F.imag**2) / (f.size**2)
        del F
        P2 = p if P2 is None else P2 + p
    return shell_average(P2, n, 3)


def spec2d_cross(A, B=None):
    n = A.shape[0]
    FA = _fft.fft2(A) / A.size
    FB = FA if B is None else _fft.fft2(B) / B.size
    P2 = np.real(FA * np.conj(FB))
    return shell_average(P2, n, 2)


# --------------------------------------------------- projected T, Q, U
def project(rho, H, axis=2, pol_perp=True):
    """Stokes maps for optically thin dust, Eqs. (4),(7),(8).

    axis is the line of sight; the two remaining axes are (horizontal,
    vertical) on the sky, in cyclic order.

    Grains align with their long axis perpendicular to H, so the emitted
    polarization direction is perpendicular to the sky-projected field:
    psi_pol = psi_H + pi/2, i.e. cos(2 psi_pol) = -cos(2 psi_H).  Eqs. (7)-(8)
    of the paper are written with psi = field angle; taking them literally
    flips the sign of Q and U (and hence of E, B and r_TE).  We use the
    physical (polarization-angle) convention, which is what reproduces the
    positive T-E correlation quoted in the paper and measured by Planck.
    """
    sgn = -1.0 if pol_perp else 1.0
    ih, iv = [(1, 2), (2, 0), (0, 1)][axis]      # sky-plane components
    Hh, Hv, Hl = H[ih], H[iv], H[axis]
    H2 = Hh**2 + Hv**2 + Hl**2
    H2 = np.where(H2 > 0, H2, 1e-300)
    T = rho.sum(axis=axis, dtype=np.float64)
    Q = sgn * (rho * (Hh**2 - Hv**2) / H2).sum(axis=axis, dtype=np.float64)
    U = sgn * (rho * (2.0 * Hh * Hv) / H2).sum(axis=axis, dtype=np.float64)
    return T, Q, U


def qu_to_eb(Q, U):
    """Flat-sky E/B: (E~ + iB~) = (Q~ + iU~) exp(-2 i theta_k),  Eq. (9)."""
    n = Q.shape[0]
    _, (kx, ky) = kmag(n, 2)
    th = np.arctan2(ky, kx)
    FQ, FU = _fft.fft2(Q), _fft.fft2(U)
    c, s = np.cos(2 * th), np.sin(2 * th)
    FE = FQ * c + FU * s
    FB = -FQ * s + FU * c
    FE[0, 0] = FB[0, 0] = 0.0
    return np.real(_fft.ifft2(FE)), np.real(_fft.ifft2(FB))


# --------------------------------------------------------------- fits
def fit_powerlaw(k, C, kmin, kmax):
    m = (k >= kmin) & (k <= kmax) & (C > 0) & np.isfinite(C)
    if m.sum() < 3:
        return np.nan, np.nan
    p = np.polyfit(np.log(k[m]), np.log(C[m]), 1)
    return np.exp(p[1]), p[0]          # amplitude, slope


def local_slope(k, C, half=1):
    lk, lc = np.log(k), np.log(np.maximum(C, 1e-300))
    s = np.full_like(lk, np.nan, dtype=float)
    for i in range(half, len(k) - half):
        s[i] = np.polyfit(lk[i-half:i+half+1], lc[i-half:i+half+1], 1)[0]
    return s


# ------------------------------------------------- full snapshot pass
def analyse(snap, axis=2):
    rho, v, H = snap["rho"], snap["v"], snap["H"]
    out = {"n": snap["n"], "meta": snap["meta"]}
    k3, Crho, _ = spec3d([rho])
    _, Cv, _ = spec3d(list(v))
    _, CH, _ = spec3d(list(H))
    out["k3"], out["Crho"], out["Cv"], out["CH"] = k3, Crho, Cv, CH

    T, Q, U = project(rho, H, axis=axis)
    E, B = qu_to_eb(Q, U)
    out["maps"] = dict(T=T, Q=Q, U=U, E=E, B=B)
    k2, CTT, _ = spec2d_cross(T)
    _, CEE, _ = spec2d_cross(E)
    _, CBB, _ = spec2d_cross(B)
    _, CTE, _ = spec2d_cross(T, E)
    _, CTB, _ = spec2d_cross(T, B)
    _, CEB, _ = spec2d_cross(E, B)
    out.update(k2=k2, CTT=CTT, CEE=CEE, CBB=CBB, CTE=CTE, CTB=CTB, CEB=CEB)

    vrms = np.sqrt(np.float64((v.astype(np.float64)**2).sum(axis=0).mean()))
    b0 = np.sqrt((H.mean(axis=(1, 2, 3), dtype=np.float64)**2).sum())
    rhobar = rho.mean(dtype=np.float64)
    out["Ms"] = vrms / snap["meta"].get("cs", 1.0)
    out["Ma"] = vrms / (b0 / np.sqrt(rhobar))
    out["sigma_lnrho"] = np.std(np.log(rho))
    return out


def finite_snapshot(s):
    """A run that blew up (or was started from one that did) writes NaN fields;
    one such snapshot turns every average it enters into NaN."""
    return bool(np.isfinite(s["rho"]).all() and np.isfinite(s["H"]).all())
