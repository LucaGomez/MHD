# 10-minute talk outline

One slide per minute, each pointing at a figure or table that exists in this
repository.  The live version of this story is the HTML report in `report/`.

| # | Slide | Show | Say |
|---|---|---|---|
| 1 | The problem | the 4-step chain from `report/index.html` | CMB B-modes are buried under Galactic dust; a CNN can learn to remove it, but it needs training maps, and the sky gives us only one realization |
| 2 | The plan | same chain, step 4 highlighted | simulate MHD turbulence, project it to Stokes maps, check the maps against real dust, then train |
| 3 | Does our solver reproduce the paper? | `results/01_own_solver_convergence.txt` | slopes converge toward Table 1 at 128 → 256 → 512; two amplitudes never do (BB/EE 0.90 vs 0.55) |
| 4 | Independent check | Table 1 comparison in `README.md` §2.1 | 100 snapshots of public 256³ boxes, a third code: Table 1 reproduced except α_BB, which misses the same way our solver does |
| 5 | Which turbulence looks like dust | `figures/well_vs_pysm/stats.png` | same estimator on PySM3 dust and on simulations; public boxes bracket dust but never match it; M_A ≈ 1.5 does |
| 6 | The 512³ stall | `figures/report/fig_cost.png` | time steps per dynamical time: 3.5k → 48k → 132k; a few near-empty cells set the step for the whole box |
| 7 | The fix, tested | limiter table in `README.md` §3 | Enzo's Alfvén cap; cost halved at 64³, mass +1.2e-4, statistics unchanged within the scatter |
| 8 | Is 256³ enough? | `figures/report/fig_resolution.png` | cascade ends at k ≈ 12; binning 256² → 128² costs 0.02 in slope; cropping native tiles costs 0.4–0.9 |
| 9 | The production plan | ℓ table in `README.md` §4 | 10° maps put the simulated range at ℓ = 108–414; several 256³ seeds at ~900 core-hours each instead of one 512³ box |
| 10 | What an agent-assisted workflow changed | `docs/AGENT_LOG.md` | measured metadata instead of trusting labels (the M_A ≈ 8 surprise); cheap experiments before expensive ones; two retracted claims kept in the log |

**Questions to expect**

* *Why not just use public data?*  No public dataset sits at M_A ≈ 1.5, which is
  the only field strength that matches dust (slide 5).
* *Is the B-mode discrepancy a bug?*  Possibly, but two independent codes miss the
  same way, and our E-mode slopes match; the next test is the paper's own pipeline
  on our boxes.
* *Why is the CNN not trained?*  Compute.  The blocking problem was the time-step
  collapse (slide 6); with it fixed the training set costs ~5 × 900 core-hours.
* *How much of this did the agent do?*  All the code and analysis in this
  repository, under review at each decision point; the decision log records what
  was tried and what was wrong.
