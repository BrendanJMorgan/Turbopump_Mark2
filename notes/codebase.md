# Codebase

## Layout
- `Python/Components/` - sizing/design pipeline. Entry: `main.py` -> `thrust_chamber_assembly()` (CEA, geometry,
  regen/film cooling) -> `powerhead()` (pumps, inducers, impellers, volute, GG, turbine, turbo-electric link).
  Inputs in `config.yaml`, loaded into module-level classes in `engine_state.py` (state is stored by mutating
  those classes, e.g. `tca.mdot = ...`).
- Turbine mirrors the pump layout: `turbines.py` holds the pitchline, velocity triangles, efficiency loop and flow
  totals; `stators.py` (`stator_flow` once, `stator_geometry` each loop pass, since the plate depends on `gg.mdot`)
  writes to `t.stator`; `rotors.py` (`rotor`, blade geometry + loss coefficients) writes to `t.rotor`. Config inputs
  sit under `turbine: stator:` / `rotor:` in `config.yaml`.
- `cad_export.py` (called from `plots()`) writes `Solidworks/Volute/<pump>_volute.step`: the volute fluid core as
  one solid (wrap lofted through `vol.contour` sections from 10 deg, fused with a ruled loft of
  `vol.nozzle_contour`). Shape-agnostic - it only reads the surface grids. Uses OpenCascade via `cadquery-ocp`.
  A single loft through wrap + cone overshoots ~4 % at the junction, hence two lofts and a fuse.
- `Python/Network/` - FullFlow cycle-balance model. Topology lives in `cycle.yaml`; `model.py` assembles it,
  `physics.py` holds the hooks. See `Python/Network/NOTES.md`. `run_cycle.py` and `primary_flow.py` are stale.
- `Python/.venv` - Python 3.11 venv with rocketcea, rocketprops, CoolProp, fullflow etc. Run with
  `Python/.venv/Scripts/python.exe`. (A bare `python` on this machine hangs - it's the Windows Store stub.)
- `Python/.github/copilot-instructions.md` is out of date (describes an `Inputs` dataclass that no longer exists).

## Units policy
Everything lives in base SI (m, kg, s, K, Pa, J, rad, rad/s). Convert at the boundary - right where a value
comes out of a library or lookup table - and nowhere else. Non-SI only for display (plots, printouts) or where
an external API demands it (FullFlow rotor speeds are rpm). **Every** numeric config value is `{value, unit}`
(use `dimensionless` for ratios/counts) and every engine_state load goes through `baseunit()`, which rejects bare
numbers. Only strings/booleans are read raw. Custom pint units in engine_state.py: `specific_speed_US`
(`Ns_US`/`Nss_US`, rpm-gpm-ft, ~1/2733) and `specific_speed_EU` (`nq`/`nss_EU`, Gulich rpm-m3/s-m); base is the
dimensionless w*sqrt(Q)/(gH)^0.75.

## Gotchas
- **rocketcea py_cea arrays:** `rockt.cstr` is **ft/s** (divide by 3.28084). `prtout.ppp` is bar, `prtout.ttt`
  is K, `prtout.wm` is g/mol (R = 8314.5 / wm J/kg-K), `therm.mw` g/mol.
- **rocketcea units:** `CEA_Obj` defaults to English units (Pc in psia, density lbm/ft3, c* ft/s, T in Rankine)
  unless `pc_units=` / `add_user_units` say otherwise. `get_full_cea_output` here passes `pc_units='bar'`.
- **rocketcea global state:** every CEA call overwrites `py_cea.prtout` / `py_cea.comp`. Read station data
  right after the run you want; don't make another CEA call in between.
- **py_cea species indexing:** `cdata.prod[j]` pairs with `comp.en[j-1]` and `therm.mw[j-1]` (off by one;
  verified against `get_SpeciesMoleFractions`). Stations are columns 0-3: injector, combustor end, throat, exit
  (finite-area run with `fac_CR`).
- **rocketprops:** `ViscAtTdegR` returns poise (x0.1 -> Pa-s), `SGLiqAtTdegR` is specific gravity.
- **Head coefficient:** `impellers.py` uses the US form psi = gH/u2^2 (half of Gulich's 2gH/u2^2). Loss
  coefficients from Gulich are normalised by u2^2 with the Gulich convention, so they need halving to pair with
  the US psi (volute losses are halved; check impeller losses in `blades.py` are treated consistently).
- **Slip factor:** stored as the fractional *reduction* (1 - gamma_Gulich), used as `u*(1-slip)`.
- **Regen coolant** runs counterflow (`flow_direction: -1`, enters at the nozzle end, index -1). The coolant
  march in `combustion_chamber_cooling.py` follows the coolant; the film march follows the gas (+x).
- The dodecane film is supercritical in the chamber (Pc 50 bar >> 18 bar) so it never boils there; it flashes
  off near the throat where gas pressure drops below the critical pressure (~x = 0.51 m).
