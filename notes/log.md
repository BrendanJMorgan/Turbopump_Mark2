# Log

## 2026-11-07 - turbine split into stator and rotor (Claude)
- Structure only. Brendan's in-progress stator lines were moved as written (still WIP, so the pipeline stops at
  `gg.T_manifold`). Checked by running pre- and post-split copies with the same two stand-ins (`c_F = 1.50`,
  old `specific_work`): all 46 shared turbine values identical, and the plot dashboard builds.

## 2026-11-07 - turbine method checked against source (Claude)
- Read NASA Turbines 1974 eqs. 1-14 and Figs. 16-20 and ran the Emmert method in a scratch script; `turbines.py`
  not changed. Open issue 3 rewritten with the specific mismatches. Fig. 16 is phi^2, not phi or psi.
- Explainer report with a live velocity triangle bench: https://claude.ai/artifact/4bVVXCd9zZR767sGAw9sfQ

## 2026-11-06 - gap triage (Claude)
- Ran `thrust_chamber_assembly()` + `powerhead()` and ranked what's missing against the output; `open_issues.md`
  rewritten as tiers. No code changed. `plots()` and `rotordynamics.py` weren't run; of the Network model only
  `run_cycle.py` was tried, and it crashes.
- Dropped as resolved: `tca.R_gas` and GG `R_gas` scaling, GG mid-rewrite items, Zweifel radians, `t.Tt`, link
  blocking the `powerhead.py` import, turbine on ox speed, film-gas thickness forced to 0, `suction_margin_factor`
  unused. The old ~2.5% fuel inducer shortfall and ~628 K coolant exit no longer match; replaced with current numbers.

## 2026-10-23 - plotting merged into plotly
- `plots.py` is now the plotly dashboard (one function per figure); `plots_plotly.py` deleted, matplotlib
  version dropped (it was already broken - imported a missing `plot_html`). Ported its wall/film temperature
  zoom, meridional reference lines, shaft power printout, and the commented-out unrolled inducer / beta plots.

## 2026-10-23 - code review and bug-fix pass (Claude)
- Full review of `Python/Components` and `Python/Network`; clear-cut bugs fixed across cooling, mixture,
  powerhead, inducers, blades, volute, pump maps, synchronous link shear formula, NOTES.md counts.
  `pump_maps.h5` regenerated (map now hits design head at the design point).
- Species block in `combustion_chamber.py` rewritten: per-station mole fractions straight from `py_cea`.
- Cooling loop non-convergence root causes: (1) film mass flow at station 0 re-added the injection every outer
  pass; (2) latent heat evaluated at film T (0 above Tc) -> divide by zero when the film goes subcritical near the
  throat. Both fixed; loop converges in ~7 passes at dx = 1 mm.
- Tried a joint (sparse matrix) film/coolant solve after N marching passes; identical results and pass counts
  at 1x/3x/6x film flow and dx down to 0.25 mm, so it was removed. Recoverable from this session's history if
  film/coolant coupling ever stalls the loop again.
- Nothing committed at time of writing.
