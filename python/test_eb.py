"""Unit tests of the projection and E/B decomposition against configurations
whose answer is known analytically."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import project, qu_to_eb, spec2d_cross

n = 64
x = (np.arange(n) + 0.5) / n
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")

def band(C, k, a=3, b=8):
    m = (k >= a) & (k <= b)
    return C[m].sum()

def case(name, rho, H, expect):
    T, Q, U = project(rho, H, axis=2)
    E, B = qu_to_eb(Q, U)
    k, CEE, _ = spec2d_cross(E)
    _, CBB, _ = spec2d_cross(B)
    _, CTE, _ = spec2d_cross(T, E)
    _, CTT, _ = spec2d_cross(T)
    ee, bb = band(CEE, k), band(CBB, k)
    rte = band(CTE, k) / np.sqrt(band(CTT, k) * ee)
    print(f"{name:<46s} B/E = {bb/ee:8.4f}   r_TE = {rte:+.3f}   [{expect}]")
    return bb / ee, rte

print("mean field along x-hat, line of sight z-hat, dust polarization _|_ H")
# 1. filaments elongated along x (density varies with y only), field along x:
#    polarization is vertical everywhere -> pure E, positively correlated with T
rho = 1.0 + 0.5 * np.cos(2 * np.pi * 4 * Y)
H = np.stack([np.ones_like(X), np.zeros_like(X), np.zeros_like(X)])
case("filament || H  (both horizontal)", rho, H, "B/E=0, r_TE=+1")

# 2. same filaments, field along y (perpendicular to the filament):
#    polarization horizontal -> pure E again, but anti-correlated with T
H2 = np.stack([np.zeros_like(X), np.ones_like(X), np.zeros_like(X)])
case("filament _|_ H", rho, H2, "B/E=0, r_TE=-1")

# 3. field at 45 deg to the filament: pure B modes, no E
s = 1 / np.sqrt(2)
H3 = np.stack([s * np.ones_like(X), s * np.ones_like(X), np.zeros_like(X)])
case("filament at 45 deg to H", rho, H3, "E=0 (pure B)")

# 4. field at 22.5 deg: equal E and B power
c, s2 = np.cos(np.pi / 8), np.sin(np.pi / 8)
H4 = np.stack([c * np.ones_like(X), s2 * np.ones_like(X), np.zeros_like(X)])
case("filament at 22.5 deg to H", rho, H4, "B/E=1")

# 5. field out of the plane of the sky (gamma=90): no polarization at all
H5 = np.stack([np.zeros_like(X), np.zeros_like(X), np.ones_like(X)])
T, Q, U = project(rho, H5, axis=2)
print(f"{'field along the line of sight':<46s} max|Q|,|U| = {abs(Q).max():.2e},"
      f" {abs(U).max():.2e}   [both 0]")
