"""Exact solution of the isothermal Riemann problem, used to validate the solver."""
import numpy as np
from scipy.optimize import brentq

def exact(x, t, rL=1.0, uL=0.0, rR=0.125, uR=0.0, cs=1.0, x0=0.5):
    # left rarefaction (u + cs ln rho invariant), right shock
    f = lambda rs: (uL - cs*np.log(rs/rL)) - (uR + cs*(rs-rR)/np.sqrt(rs*rR))
    rs = brentq(f, 1e-8, 100.0)
    us = uL - cs*np.log(rs/rL)
    S  = uR + cs*np.sqrt(rs/rR)          # shock speed
    xi = (x-x0)/t
    rho = np.empty_like(xi); u = np.empty_like(xi)
    head, tail = uL-cs, us-cs
    for i, s in enumerate(xi):
        if s < head:
            rho[i], u[i] = rL, uL
        elif s < tail:
            rho[i] = rL*np.exp(-(s+cs-uL)/cs); u[i] = s+cs
        elif s < S:
            rho[i], u[i] = rs, us
        else:
            rho[i], u[i] = rR, uR
    return rho, u, rs, us, S

if __name__ == "__main__":
    import sys
    for fn in sys.argv[1:]:
        d = np.loadtxt(fn)
        x, rho, vx = d[:,0], d[:,1], d[:,2]
        re, ue, rs, us, S = exact(x, 0.1)
        m = (x>0.25)&(x<0.8)
        print(f"{fn}: N={len(x)}  L1(rho)={np.abs(rho[m]-re[m]).mean():.5e} "
              f" L1(v)={np.abs(vx[m]-ue[m]).mean():.5e}  rho*={rs:.4f} u*={us:.4f} S={S:.4f}")
