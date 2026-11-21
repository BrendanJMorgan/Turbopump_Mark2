# Open issues

Triage last redone 2026-11-06 from a live run of `thrust_chamber_assembly()` + `powerhead()` (config with the
3 mm chamber wall). Numbers below are from that run. Hand estimates are marked "est." and should be rechecked.

## Tier 1 - changes the design numbers
Suggested order: 1 -> 2 -> 3 -> GG flow iteration (tier 2) -> 4. Item 5 is independent.

1. **Shaft speed has no suction check.** Required Nss is ~53,000 US (ox) / ~50,000 (fuel); the Gulich 7.23a
   estimate gives ~25,000 / ~24,000 (`inducers.py:40`, `:160`). Both are computed, neither is compared.
   Brumfield sets speed from NPSH, so required Nss stays about the same if tank pressure rises: the selection
   rule is what has to change. The check also applies `flow_margin_factor` to the required side only.
   Limiting to the realistic value puts ox near 12k rpm (now 25.9k) and roughly doubles impeller radius (est.).
2. **Pump power is too low.** Ox hydraulic efficiency reads 0.963 because `vaned_diffuser.py` is unfinished and
   `total_loss = 0.0`; fuel reads 0.787. `pumps.py:58` also omits leakage and disk friction, and `inlet_loss` /
   `outlet_loss` are 0. Turbine power (106 kW) could rise ~30% (est.).
3. **Turbine work model misreads the Emmert method (NASA Turbines 1974, eqs. 1-14).** Eq. 5b is the energy available
   to one element and sets its outlet velocity (eqs. 6-7); it is not stage work. Work is eq. 10b, efficiency is eq. 11
   (work / dh_s). The efficiency loop in `turbines.py` uses 5b as work over `dh_s + KE`, giving 1,222 kJ/kg from 899 and a false 71%.
   Other mismatches: Fig. 16 reads phi^2, so eq. 3 is `2*coeff - 1`, not `2*coeff**2 - 1` (`rotors.py:48`); `2u` triangle
   forces a 44 deg nozzle (from axial); `c_tau = 1.50` should be replaced by c1 from eqs. 2/6; `K` should be w2/w1 and
   windage (eq. 14) is 5.6e-5*(U/c0)^2*rho*N*(1-F)/F, not the code's `1.539e-9` linear form. Book units for eq. 12
   aren't stated on that page; rpm and lb/ft3 gives ~1%. Scratch Emmert run: 44 deg -> eta 0.375, 0.56 kg/s;
   65 deg, U/c0 0.35 -> 0.58, 0.36 kg/s (needs ~62% admission at 6.35 mm blade); 70 deg -> 0.61, 0.345 kg/s (~73%).
   Fig. 20 shows ~60% at 31% admission, for a different low-PR turbine. GG flow depends strongly on
   `gg.c_star_efficiency` (0.85 -> 0.27 kg/s). Item 1's lower speed may hit `r_pitchline_max`, so recheck U/c0 after it.
4. **Electric link is not sized.** `synchronous_links.py` is unfinished (blank inductance terms, syntax errors)
   and commented out in `powerhead.py`, so link efficiency is effectively 1. By its own shear formula the margin
   is 0.38 (motor) / 0.41 (alternator), so it will raise as soon as it's enabled. Lower speed from item 1 likely
   fixes this (torque x2, radius x2 -> shear ~1/4, est.). Also: resistivity table values don't look like copper
   resistivity; wire-gauge pick compares an area to a diameter; `np.interp(m.temperature, curve)` has wrong args.
5. **Coolant inlet is ~100 K too hot.** The 400 K preset (`combustion_chamber_cooling.py:73`) is never replaced
   by the real inlet temperature: coolant enters at 400 K instead of ~293 K and exits at 460 K.

## Tier 2 - needed to close the design
- **Flow balance:** GG fraction guess 5%, result 3.66%, no iteration (`powerhead.py:26`). Fuel pump 3.3% oversized;
  `gg.A_throat` stays at the guessed flow.
- **Injectors:** no module; dp is a flat 25% of Pc (`engine_state.py:62`). GG pressure and feed are also a guess.
- **Feed system:** pump inlet = tank pressure, no line loss or hydrostatic head; ox downstream is a flat 10 psi
  allowance. Fuel tank at 0.5 bar is below ambient - decide whether tank values are gauge or absolute.
- **Structures:** no checks. Turbine rim ~380 m/s at ~990 K (disk ~500 MPa or more, est.). Wall went 6 -> 3 mm with
  no hoop/thermal check; hot wall 677 K, copper properties fixed. Shaft size isn't calculated in the pipeline.
- **Shaft support:** axial thrust, bearings, seals not modelled (`shaft_efficiency = 1`).
- **Rotordynamics:** `rotordynamics.py` is stale - one 1.5 m shaft at 20k rpm, hardcoded dimensions; the design is
  now two shafts joined by the link.
- **Network model:** `run_cycle.py` crashes (`mark2_cycle` missing). Design point and pump maps are old placeholders
  (ox 10.9 kg/s vs 9.1), and nothing feeds sizing results into it.

## Tier 3 - before hardware
- **Film model:** no entrainment / effectiveness model, film heats to ~1180 K, and it takes ~60% of gas-side heat
  (1.86 MW in, 0.72 MW to regen). Dominates wall temperature far more than the solver does.
- Transients and off-design: start, shutdown, throttling.
- Combustion stability; GG chamber sizing (`gg.r_chamber` unused).
- A single pass/fail margin report and a few regression tests.

## Smaller bugs
- `rotors.py:50`: incidence loss looked up at 26 deg; actual incidence is 0 (blade angle = flow angle).
- `gg.c_tau_efficiency` is never used. GG exhaust thrust is in the config comment but not in the sizing.
- `combustion_chamber_cooling.py`: 0/0 film resistances give NaN (harmless in current branches, warns on every
  run). Film properties use `gas.velocity` as film velocity ("VELOCITY PROBABLY WRONG").
- `powerhead.py`: `turbine.gear_efficiency` doesn't exist (`common_shaft: false` branch).
- `Network/components.py:80`: `SynchronousLink.electrical_frequency` reports mechanical rev/s (needs pole pairs).
- `pumps.py` efficiency mixes conventions unless impeller losses are also halved (see codebase.md). Not rechecked.
- `inducers.py`: `ind.blockage = 0.98` is a dummy; inducer exit swirl isn't passed to the impeller.

## Repo hygiene
- `Python/.venv/` (~7.7k files), `__pycache__`, and `Documents/~$Diagrams.pptx` are tracked in git despite `.gitignore`.
