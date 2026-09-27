/* Isothermal ideal MHD with Dedner GLM divergence cleaning.
 * Scheme follows Stalpes, Collins & Huffenberger (2024, arXiv:2404.02874):
 *   - piecewise-linear reconstruction
 *   - HLLD Riemann solver for isothermal MHD (Mignone 2007, JCP 225, 1427)
 *   - Dedner et al. (2002) mixed hyperbolic/parabolic divergence cleaning
 *   - stochastic OU driving (Federrath et al. 2010) at constant energy
 *     injection rate (Mac Low & Klessen 2004)
 */
#ifndef MHD_H
#define MHD_H

#define NG 2      /* ghost zones */
#define NV 8      /* rho, mx, my, mz, bx, by, bz, psi */
#define IRHO 0
#define IMX  1
#define IMY  2
#define IMZ  3
#define IBX  4
#define IBY  5
#define IBZ  6
#define IPSI 7

typedef struct {
    int n;            /* cells per side (interior) */
    int nt;           /* n + 2*NG */
    long ncell;       /* nt^3 */
    double L, dx;
    double cs;        /* isothermal sound speed */
    double cfl;
    double ch;        /* GLM hyperbolic speed (set each step) */
    double glm_alpha; /* psi damping parameter */
    double rho_floor;
    /* driving */
    double ms_target, ma_target;
    double edot;      /* energy injection rate */
    double tdrive;    /* OU correlation time */
    double zeta;      /* solenoidal fraction parameter */
    unsigned long seed;
    /* control */
    double tmax;
    double dt_snap;
    char outdir[512];
} Params;

#endif
