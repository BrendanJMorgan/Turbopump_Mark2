# Mark 2 cycle model — scope and assumptions

FullFlow 2.0.1 + ThermoProp 2.1.0.

```
cycle.yaml        the network: parameters, states, elements, and an ordered
                  `build` manifest of physics hooks + element refs
physics.py        named `build` hooks — the equations that aren't one component
                  (pump torque, dissipation heating, CEA chain, shaft balance,
                  performance)
components.py     SynchronousLink, ChokedFlowParameter, ImpulseTurbineStage
model.py          build_cycle(config="cycle.yaml") — generic assembler + resolver
make_pump_maps.py regenerates pump_maps.h5 (the large (Q, N) grids live there,
                  not in the YAML)
diagram.py        wiring diagram straight from cycle.yaml (no import of model.py)
run_cycle.py      driver  (STALE: still imports the old `mark2_cycle` package)
```

The topology is data.  To change the cycle, edit `cycle.yaml`; only add Python
in `physics.py` when a new derived quantity needs an equation.  `inputs.py`
(`CycleInputs`) is superseded by the YAML `parameters` / `states` blocks.

## What the solve actually determines

19 unknowns, 19 residuals:

| unknown | residual |
| --- | --- |
| 7 branch mass flows (feed lines, crossover, mainline, jacket) | DarcyWeisbach momentum |
| 9 node pressures | Volume mass balances |
| 2 pump mass flows | pump-curve discharge-pressure error |
| fuel shaft speed | fuel shaft torque balance |
| turbine torque | turbine efficiency error |

Chamber and GG pressures fall out of the mass balances against nozzle and
turbine flow. Injectors and the MFV are algebraic — they produce mass flow from
ΔP and add no unknown.

The ox shaft speed is **not** an unknown. The link is pole-locked, so
`N_ox = N_fuel * pole_pair_ratio` is a kinematic constraint, and the ox torque
balance is satisfied identically by the link component. One shaft, one speed.

## Deliberate simplifications

- **Secondary flows are absent.** No bearing coolant tap, no GG film, no seal
  purge, no tank pressurisation flow. The helium bottle and press valve collapse
  into a fixed regulated ullage pressure per tank.
- **Regen jacket is one lumped path**: one ΔP, one heat load. The station
  solution stays in `combustion_chamber_cooling.py`; only ΔP and ΔT feed back
  into the cycle.
- **Film coolant** is carried through the nozzle mass balance but excluded from
  the core mixture ratio. There is no mixing model, so chamber gas properties
  are core values while the throat passes total flow.
- **GG exhaust thrust is not included.** The dump nozzle isn't modelled.
- **Pump torque comes from the map, not from `ConstantDensityPump`.** That
  component takes torque as an input and back-computes efficiency, which is
  backwards for design work — so torque is formed as `g·H·ṁ/(η·ω)` and handed in.
  Same inversion for the turbine, but there it needs a `Balance` because
  `GasTurbine` computes mass flow internally.
- **c\* efficiency is applied as η² on total temperature**, since c\* ∝ √(RT₀).
  This gives η on c\* and the correct mass flow at a given Pc.
- **`SynchronousLink` is a power-balance model only.** No torque-angle relation,
  so it cannot predict pull-out — it will transmit whatever torque is asked of
  it. Pull-out margin has to be checked separately against `synchronous_links.py`.
- **Fuel transport properties use RP-1** (ThermoProp has no liquid Jet-A
  transport data), while combustion uses `Jet-A(L)`. Same split as `config.yaml`.

## Known issues

- **Residual scaling.** Mass balances are in kg/s and momentum residuals in
  N; they go into the same least-squares vector. With the placeholder design
  point the solve stalls on `xtol` at max |residual| ≈ 0.16 kg/s (~1.5% of
  branch flow) rather than driving to tolerance. `jacobian_method="3-point"`
  with tight tolerances is prohibitively slow because every residual call runs
  two CEA equilibrium solves. Worth either normalising the node balances by a
  reference flow or accepting a looser `rtol`.
- **Map bounds.** `extrapolate=False`, so a solver excursion outside the pump
  map raises rather than silently extrapolating. Good for catching bad guesses,
  but the map has to be wide enough to contain the transient path if this model
  is later run with `Transient`.
- **Cost.** Two `Equilibrium` lookups per residual evaluation dominate runtime.
  If it gets painful, tabulate T₀(OF, p), γ, R offline and swap the `Lookup`
  for a `Map`.

## Transient

`Rotor`, `Volume`, and `DarcyWeisbach` all already expose dynamics, so the same
network runs under `Transient` once `volume=` is supplied to the nodes and
`effective_area=` to the lines. The GSE-powered start would enter as an
additional drive torque on the fuel shaft, sequenced with `Sequence`/`Actuator`.
