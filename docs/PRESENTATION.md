# 10-minute talk outline

One slide per minute, each pointing at a figure or table that exists in this
repository.  The live version of this story is the HTML report in `report/`.

| # | Slide | Show | Say |
|---|---|---|---|
| 1 | The problem | the 4-step chain from `report/index.html` | CMB B-modes are buried under Galactic dust; a CNN can learn to remove it, but it needs training maps, and the sky gives us only one realization |
| 2 | The plan | same chain, step 4 highlighted | simulate MHD turbulence, project it to Stokes maps, check the maps against real dust, then train |
| 3 | Does our solver reproduce the paper? | `results/01_own_solver_convergence.txt` | slopes converge toward Table 1 at 128 → 256 → 512; two amplitudes never do (BB/EE 0.90 vs 0.55) |
| 4 | Independent check | Table 1 comparison in `README.md` §2.1 | 100 snapshots of public 256³ boxes, a third code: Table 1 reproduced except α_BB, which misses the same way our solver does |
| 5 | Do simulated maps look like dust? | `figures/well_vs_pysm/stats.png` | same estimator on PySM3 dust and on simulations: **none of the boxes match** — E-modes too shallow, BB/EE ≈ 1 against 0.6, T–E correlation off in a field-dependent way |
| 6 | The 512³ stall | `figures/report/fig_cost.png` | time steps per dynamical time: 3.5k → 48k → 132k; a few near-empty cells set the step for the whole box |
| 7 | The fix, tested | limiter table in `README.md` §3 | Enzo's Alfvén cap; cost halved at 64³, mass +1.2e-4, statistics unchanged within the scatter |
| 8 | Is 256³ enough? | `figures/report/fig_resolution.png` | cascade ends at k ≈ 12; binning 256² → 128² costs 0.02 in slope; cropping native tiles costs 0.4–0.9 |
| 9 | The production plan | ℓ table in `README.md` §4 | 10° maps put the simulated range at ℓ = 108–414; several 256³ seeds at ~900 core-hours each instead of one 512³ box |
| 10 | Cross-agent review | `review/round1-claude.md` | a second agent (Codex) reviewed the work: 10 findings, 2 blocking, one of which **removed a headline claim** — the apparent M_A ≈ 1.5 match came from measuring cropped tiles in the numerically damped range |

**The one slide people will ask about**

Slide 5 used to say the opposite. The review found that the tiles feeding the Enzo
row sampled box wavenumbers 12–52 — past the dissipation knee — where slopes
steepen and happened to land near the dust values. On properly binned faces the
match disappears. Say this plainly: it is the strongest evidence in the talk that
the process has error-correction in it.

**Questions to expect**

* *Why not just use public data?*  No public dataset sits at M_A ≈ 1.5, the
  paper's Planck point; and after the review, no box we have — public or ours —
  reproduces the dust models' E/B statistics, which is now the open problem.
* *Is the B-mode excess a bug?*  It could be: all three codes were analysed with
  one pipeline, and the review already found one estimator bug. The next test
  pushes known pure-E realizations through the whole projection path.
* *Why is the CNN not trained?*  Compute.  The blocking problem was the time-step
  collapse (slide 6); with it fixed the training set costs ~5 × 900 core-hours.
* *How much of this did the agent do?*  All the code and analysis in this
  repository, under review at each decision point; the decision log records what
  was tried and what was wrong.
