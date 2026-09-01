"""
PLACEHOLDER design point for the Mark 2 cycle model.

Every number here is a stand-in.  Replace each block with the corresponding
output of the forward-solve code; the `source:` comments name the module.  The
pump maps in particular are affinity-scaled from a single point and are NOT a
substitute for sweeping impellers.py over a (Q, N) grid.
"""

from __future__ import annotations

import numpy as np

from mark2_cycle import CycleInputs, Propellants

# ---------------------------------------------------------------------------
# design point (config.yaml)
# ---------------------------------------------------------------------------
CHAMBER_PRESSURE = 5.0e6      # Pa
DESIGN_SPEED = 30000.0        # rpm  (~3142 rad/s shaft speed cap)
OXIDIZER_FLOW = 10.885        # kg/s
FUEL_FLOW = 5.575             # kg/s
OXIDIZER_DENSITY = 1142.0     # kg/m^3
FUEL_DENSITY = 810.0          # kg/m^3
OXIDIZER_HEAD = 531.0         # m
FUEL_HEAD = 995.0             # m
OXIDIZER_EFFICIENCY = 0.65
FUEL_EFFICIENCY = 0.60


def placeholder_map(design_flow, design_head, design_efficiency, design_speed) -> dict:
    """
    Affinity-scaled synthetic map.  PLACEHOLDER.

        H(Q, N) = H_d * (N/N_d)^2 * psi(phi)
        eta(Q, N) = eta_d * (1 - 1.6 * (phi - 1)^2)
        phi = (Q / Q_d) / (N / N_d)
    """
    flow = np.linspace(0.15 * design_flow, 1.9 * design_flow, 25)
    speed = np.linspace(0.25 * design_speed, 1.3 * design_speed, 15)

    phi = (flow[:, None] / design_flow) / (speed[None, :] / design_speed)
    speed_ratio = (speed[None, :] / design_speed) ** 2

    head = design_head * speed_ratio * (1.15 - 0.15 * phi - 0.20 * phi**2)
    efficiency = design_efficiency * (1.0 - 1.6 * (phi - 1.0) ** 2)
    efficiency = np.clip(efficiency, 0.05, 0.95)

    return dict(
        volumetric_flow=flow,
        rotor_speed=speed,
        head_rise=head,
        efficiency=efficiency,
    )


def design_point() -> CycleInputs:
    return CycleInputs(
        ambient_pressure=93900.0,
        ambient_temperature=293.0,
        gravity=9.80665,
        propellants=Propellants(
            oxidizer_fluid="Oxygen",
            fuel_propellant="RP-1",     # transport data; ThermoProp has no Jet-A liquid props
            oxidizer_reactant="O2(L)",
            fuel_reactant="Jet-A(L)",
        ),
        # ---- tanks ----
        oxidizer_tank_pressure=2.5e5,
        oxidizer_tank_temperature=90.188,
        fuel_tank_pressure=1.0e5,
        fuel_tank_temperature=293.0,
        # ---- oxidizer side ----
        oxidizer_feed_line=dict(
            length=0.60, hydraulic_diameter=0.050,
            cross_sectional_area=np.pi / 4 * 0.050**2,
            friction_factor=0.125,   # lumped f*L/D for the whole run
        ),
        oxidizer_pump_map=placeholder_map(
            OXIDIZER_FLOW / OXIDIZER_DENSITY, OXIDIZER_HEAD, OXIDIZER_EFFICIENCY, DESIGN_SPEED
        ),
        oxidizer_crossover=dict(
            length=0.20, hydraulic_diameter=0.030,
            cross_sectional_area=np.pi / 4 * 0.030**2,
            friction_factor=0.290,
        ),
        tca_oxidizer_injector=dict(discharge_coefficient=0.75, cross_sectional_area=3.00e-4),
        gg_oxidizer_injector=dict(discharge_coefficient=0.75, cross_sectional_area=2.66e-6),
        # ---- fuel side ----
        fuel_feed_line=dict(
            length=0.60, hydraulic_diameter=0.040,
            cross_sectional_area=np.pi / 4 * 0.040**2,
            friction_factor=0.059,
        ),
        fuel_pump_map=placeholder_map(
            FUEL_FLOW / FUEL_DENSITY, FUEL_HEAD, FUEL_EFFICIENCY, DESIGN_SPEED
        ),
        fuel_mainline=dict(
            length=0.80, hydraulic_diameter=0.025,
            cross_sectional_area=np.pi / 4 * 0.025**2,
            friction_factor=0.085,
        ),
        main_fuel_valve=dict(discharge_coefficient=0.80, cross_sectional_area=3.05e-4),
        regen_jacket=dict(
            length=1.20, hydraulic_diameter=0.0028,
            cross_sectional_area=64 * 0.0028 * 0.002032,
            friction_factor=0.031,
        ),
        regen_heat_load=1.2e6,      # W   source: combustion_chamber_cooling.py
        tca_fuel_injector=dict(discharge_coefficient=0.75, cross_sectional_area=1.62e-4),
        gg_fuel_injector=dict(discharge_coefficient=0.75, cross_sectional_area=1.20e-5),
        film_coolant_fraction=0.10,
        # ---- thrust chamber ----
        throat_area=5.0e-3,          # m^2  source: combustion_chamber_geometry.py
        expansion_ratio=4.0,
        c_star_efficiency=0.85,
        thrust_coefficient_efficiency=0.95,
        # ---- gas generator ----
        gg_c_star_efficiency=0.75,
        # ---- turbine ----
        turbine_nozzle_throat_area=1.435e-4,
        turbine_nozzle_discharge_coefficient=0.95,
        turbine_efficiency=0.40,
        turbine_exit_pressure=1.05e5,
        # ---- synchronous link ----
        pole_pair_ratio=1.0,
        alternator_efficiency=0.94,
        cable_efficiency=0.98,
        motor_efficiency=0.93,
        # ---- rotordynamics ----
        fuel_shaft_inertia=2.0e-3,
        oxidizer_shaft_inertia=1.5e-3,
        # ---- initial guesses ----
        guess=dict(
            fuel_shaft_speed=DESIGN_SPEED,
            oxidizer_mass_flow=OXIDIZER_FLOW,
            fuel_mass_flow=FUEL_FLOW,
            regen_mass_flow=FUEL_FLOW * 0.88,
            ox_pump_inlet_pressure=2.3e5,
            ox_pump_discharge_pressure=6.2e6,
            ox_manifold_pressure=6.0e6,
            fuel_pump_inlet_pressure=0.9e5,
            fuel_pump_discharge_pressure=8.0e6,
            fuel_mainline_pressure=7.8e6,
            fuel_manifold_pressure=7.5e6,
            jacket_outlet_pressure=6.0e6,
            chamber_pressure=CHAMBER_PRESSURE,
            gas_generator_pressure=4.0e6,
            turbine_torque=60.0,
        ),
    )
