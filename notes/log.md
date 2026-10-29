# Log

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
