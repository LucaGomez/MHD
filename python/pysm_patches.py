"""Flat 128x128 T/Q/U dust patches from PySM3 models at 353 GHz, plus a full-sky
reference measurement of the same model (curved-sky pseudo-C_ell).

    python3 python/pysm_patches.py --models d1 d10 d12 --out pysm_patches

Patches: gnomonic projections, 128 px of RESO arcmin (default 9.375', i.e.
20 deg), centred on HEALPix nside=4 pixel centres with LATMIN <= |b| <= LATMAX.
The map is smoothed first (Gaussian, FWHM = 2 pixels) so the nearest-pixel
projection does not alias; E and B are also decomposed on the full sphere and
projected as scalar maps, which avoids E->B leakage at the patch edges.  Q/U are rotated from HEALPix (COSMO) to the IAU
convention and to the flat patch's (x = east, y = north) axes, so that
python/analysis.qu_to_eb applies with the same convention as the MHD maps
(checked against the curved-sky result: positive TE, BB/EE ~ 0.5).

Full-sky reference: healpy.anafast on |b| >= LATMIN with an apodized mask;
slopes over ell 40-600, BB/EE and r_TE over the same range (Planck 2018 XI
measured alpha_EE = -2.42, alpha_BB = -2.54, BB/EE = 0.53, r_TE = 0.36).
"""
import argparse, os, time
import numpy as np
import healpy as hp
import astropy.units as u


def sky_map(model, nside):
    if model == "d12":
        return d12_353(nside)
    import pysm3
    sky = pysm3.Sky(nside=nside, preset_strings=[model], output_unit="uK_RJ")
    m = sky.get_emission(353 * u.GHz)
    return np.asarray(m.value, dtype=np.float64)


def d12_353(nside, nlayers=6, color_correction=0.911):
    """d12 (MKD layered dust) at its own reference frequency, 353 GHz, where
    every layer's modified black body factor is 1: the emission is just the sum
    of the layer amplitudes.  Built one layer at a time and downgraded as we go
    -- pysm3.Sky holds all six nside=2048 layers at once and needs ~5 GB."""
    from astropy.utils.data import download_file
    base = "https://portal.nersc.gov/project/cmb/pysm-data/mkd_dust/2048/thermaldust_ampl{}.fits"
    tot = np.zeros((3, hp.nside2npix(nside)))
    for i in range(1, nlayers + 1):
        lay = hp.read_map(download_file(base.format(i), cache=True), field=None, dtype=np.float32)
        tot += hp.ud_grade(np.asarray(lay, dtype=np.float64), nside)
        del lay
    return tot * color_correction          # MJy/sr; only ratios are used


def fullsky_reference(m, latmin, nside):
    npix = hp.nside2npix(nside)
    th, ph = hp.pix2ang(nside, np.arange(npix))
    b = 90.0 - np.degrees(th)
    mask = (np.abs(b) >= latmin).astype(float)
    mask = hp.smoothing(mask, fwhm=np.radians(5.0), lmax=3 * 64)   # soft edge
    mask = np.clip((mask - 0.5) * 2, 0, 1)
    fsky = mask.mean()
    cl = hp.anafast(m * mask, lmax=700, pol=True)       # TT EE BB TE EB TB
    ell = np.arange(cl.shape[1])
    s = (ell >= 40) & (ell <= 600)
    fit = lambda c: np.polyfit(np.log(ell[s]), np.log(np.abs(c[s])), 1)[0]
    # D_ell-free power law in C_ell, same as Planck's alpha (C_ell ~ ell^alpha)
    out = dict(fsky=fsky, aTT=fit(cl[0]), aEE=fit(cl[1]), aBB=fit(cl[2]),
               bbee=cl[2][s].sum() / cl[1][s].sum(),
               rte=cl[3][s].sum() / np.sqrt(cl[0][s].sum() * cl[1][s].sum()))
    return out, ell, cl


def patches(m, nside, reso, npx, latmin, latmax):
    fwhm = np.radians(2 * reso / 60.0)
    lmax = min(3 * nside - 1, int(4 * np.pi / fwhm) * 2)
    alm = hp.map2alm(m, lmax=lmax, pol=True)
    bl = hp.gauss_beam(fwhm, lmax=lmax)
    alm = [hp.almxfl(a, bl) for a in alm]
    ms = hp.alm2map(alm, nside, lmax=lmax, pol=True)
    # E and B as scalar maps, decomposed on the full sphere (no patch leakage)
    me = hp.alm2map(alm[1], nside, lmax=lmax)
    mb = hp.alm2map(alm[2], nside, lmax=lmax)
    del alm
    th, ph = hp.pix2ang(4, np.arange(hp.nside2npix(4)))
    lon, lat = np.degrees(ph), 90 - np.degrees(th)
    sel = (np.abs(lat) >= latmin) & (np.abs(lat) <= latmax)
    T, Q, U, E, B, cen = [], [], [], [], [], []
    for lo, la in zip(lon[sel], lat[sel]):
        proj = hp.projector.GnomonicProj(rot=(lo, la, 0.0), xsize=npx, ysize=npx, reso=reso)
        f = lambda x, y, z: hp.vec2pix(nside, x, y, z)
        t = proj.projmap(ms[0], f); q = proj.projmap(ms[1], f); uu = proj.projmap(ms[2], f)
        e = proj.projmap(me, f); bb = proj.projmap(mb, f)
        # true sky coordinates of every pixel, on the same grid as projmap
        x, y = proj.ij2xy()
        tth, pph = proj.xy2ang(np.ravel(x), np.ravel(y))
        plat = (90 - np.degrees(tth)).reshape(t.shape)
        plon = np.degrees(pph).reshape(t.shape)
        plon = (plon - lo + 180) % 360 - 180            # continuous around the centre
        t, q, uu, e, bb, plat, plon = orient([t, q, uu, e, bb, plat, plon], plat, plon)
        # HEALPix E is the physical E mode (dust: TE > 0), which is also what
        # our flat convention gives, so it is used as is.  The (east, north)
        # frame is a mirror image of the sky as seen from inside, so the sign
        # of B (hence of TB, EB) is not meaningful here -- nor for MHD boxes.
        # HEALPix (COSMO) -> our convention (angle from +x=east towards
        # +y=north): Q -> -Q, U -> -U.
        q, uu = -q, -uu
        # rotate from the local (east, north) basis at each pixel to the
        # patch axes: alpha = direction of local east in patch coordinates
        gx, gy = np.gradient(plat)                      # d lat / d x, d lat / d y
        alpha = np.arctan2(gy, gx) - np.pi / 2
        c, s_ = np.cos(2 * alpha), np.sin(2 * alpha)
        q, uu = c * q - s_ * uu, s_ * q + c * uu
        T.append(t); Q.append(q); U.append(uu); E.append(e); B.append(bb); cen.append((lo, la))
    f32 = lambda a: np.array(a, np.float32)
    return f32(T), f32(Q), f32(U), f32(E), f32(B), np.array(cen)


def orient(arrs, plat, plon):
    """Transpose/flip so that axis 0 = x points east (longitude increases) and
    axis 1 = y points north (latitude increases), at the patch centre."""
    n = plat.shape[0]; c = n // 2
    dlat0 = plat[c + 1, c] - plat[c - 1, c]
    dlat1 = plat[c, c + 1] - plat[c, c - 1]
    if abs(dlat0) > abs(dlat1):          # latitude runs along axis 0: swap axes
        arrs = [a.T for a in arrs]
        plat, plon = plat.T, plon.T
    if plat[c, c + 1] < plat[c, c - 1]:
        arrs = [a[:, ::-1] for a in arrs]; plon = plon[:, ::-1]
    if plon[c + 1, c] < plon[c - 1, c]:
        arrs = [a[::-1, :] for a in arrs]
    return arrs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", nargs="+", default=["d1", "d10", "d12"])
    ap.add_argument("--nside", type=int, default=1024)
    ap.add_argument("--reso", type=float, default=9.375)
    ap.add_argument("--npx", type=int, default=128)
    ap.add_argument("--latmin", type=float, default=35)
    ap.add_argument("--latmax", type=float, default=75)
    ap.add_argument("--out", default="pysm_patches")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for model in a.models:
        t0 = time.time()
        m = sky_map(model, a.nside)
        ref, ell, cl = fullsky_reference(m, a.latmin, a.nside)
        T, Q, U, E, B, cen = patches(m, a.nside, a.reso, a.npx, a.latmin, a.latmax)
        # beam sigma in units of the patch side, for patch_stats(beam_sigma=)
        bsig = (2 * a.reso / 2.3548) / (a.reso * a.npx)
        np.savez_compressed(os.path.join(a.out, f"{model}.npz"), T=T, Q=Q, U=U, E=E, B=B,
                            centers=cen, beam_sigma=bsig,
                            reso=a.reso, nside=a.nside, ell=ell, cl=cl,
                            **{f"ref_{k}": v for k, v in ref.items()})
        print(f"{model}: {len(T)} patches, full sky |b|>={a.latmin:g} (fsky {ref['fsky']:.2f}), "
              f"ell 40-600: a_TT={ref['aTT']:.2f} a_EE={ref['aEE']:.2f} a_BB={ref['aBB']:.2f} "
              f"BB/EE={ref['bbee']:.2f} r_TE={ref['rte']:.2f}   ({time.time()-t0:.0f}s)", flush=True)
        del m


if __name__ == "__main__":
    main()
