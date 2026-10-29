# Open issues

## Bugs / likely bugs (not yet fixed)
- `combustion_chamber.py`: `tca.R_gas = 1000 * 8314.5 / wm` is 1000x too large (`8314.5 / wm` is already J/kg-K). Unused so far.
- `gas_generator.py` is mid-rewrite: uses `gg.A_throat` before it's set; `gg.R`, `gg.Tc`, `gg.gamma` no longer
  defined but `turbines.py` reads them; `R_gas = 8.3145 / wm` is 1000x too small.
- `turbines.py`: Zweifel lookup passes radians (`np.pi/2 - t.blade_angle`) into a degree table; `t.Tt` doesn't
  exist on the turbine object (probably `gg.Tt`).
- `combustion_chamber_cooling.py`: film-gas thickness is computed then forced to 0; 0/0 film resistances give
  NaN (harmless in current branches). Film properties use `gas.velocity` as film velocity ("VELOCITY PROBABLY WRONG").
- `synchronous_links.py` is unfinished (syntax errors) - blocks `powerhead.py` import, so `main.py` can't run end
  to end. Also: resistivity table values don't look like copper resistivity; wire-gauge pick compares an area
  to a diameter; `np.interp(m.temperature, curve)` has wrong arguments.
- `powerhead.py`: turbine runs at ox shaft speed but sits on the fuel shaft; `turbine.gear_efficiency` doesn't exist.
- `Network/components.py`: `SynchronousLink.electrical_frequency` reports mechanical rev/s (needs pole pairs).

## Modeling decisions waiting on Brendan
- Inducer speed selection (Brumfield) vs check (Gulich 7.23a + clearance penalty) use different models, and the
  check applies `flow_margin_factor` (1.2) only to the *required* side -> fuel inducer fails by ~2.5%.
  `suction_margin_factor` is unused. Options: same flow on both sides; select speed against the Gulich check.
- Supercritical film model: film heats to ~1240 K then flashes at the throat. Needs an entrainment / film
  effectiveness model; this dominates wall temperature far more than the solver does.
- Regen coolant exits ~628 K - near dodecane's pseudocritical region (~660 K) and coking territory.
- `pumps.py` efficiency mixes conventions unless impeller losses are also halved (see codebase.md).

## Repo hygiene
- `Python/.venv/` (~7.7k files), `__pycache__`, and `Documents/~$Diagrams.pptx` are tracked in git despite `.gitignore`.
