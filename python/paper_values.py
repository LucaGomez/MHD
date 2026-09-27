"""Published values from Stalpes, Collins & Huffenberger (2024), for comparison."""

# Table 1: linear fits  alpha_q = a + b M_S + c M_A
TABLE1 = {
    "alpha_rho": dict(a=-3.61, b=0.16, c=-0.00),
    "alpha_v":   dict(a=-3.86, b=0.02, c= 0.14),
    "alpha_H":   dict(a=-3.31, b=0.02, c=-0.28),
    "alpha_TT":  dict(a=-3.66, b=0.15, c= 0.09),
    "alpha_EE":  dict(a=-3.63, b=0.17, c= 0.28),
    "alpha_BB":  dict(a=-4.82, b=0.28, c= 0.64),
}

# ranges quoted in the text (section 3.1-3.3, 5)
RANGES = {
    "alpha_rho": (-3.5, -2.5),      # M_S = 0.5 -> 7
    "alpha_v":   (-3.9, -3.5),
    "alpha_H":   (-3.75, -3.3),     # high M_A -> low M_A
    "alpha_EE":  (-3.5, -2.3),      # M_S = 0.5 -> 6
    "alpha_BB":  (-3.5, -2.1),      # weakly magnetised
}

# amplitude ratios and correlations at high sonic Mach number
HIGH_MS = {
    "A_EE/A_TT": 0.62,
    "A_BB/A_TT": 0.34,
    "A_BB/A_EE": (0.55, 0.07),
    "r_TE": 0.3,
}

INFERRED = dict(Ms=4.7, Ma=1.5)          # from alpha_EE=-2.4, alpha_BB=-2.5

PLANCK = dict(aEE=-2.42, aBB=-2.54, BB_EE=0.53, rTE=0.355, rTB=0.055)

# sign statistics of the parity-violating correlations
PARITY = {"rTB_positive": (14, 21), "rTB_positive_highMs": (8, 9),
          "rEB_positive": (14, 21), "rEB_positive_highMs": (9, 9)}
