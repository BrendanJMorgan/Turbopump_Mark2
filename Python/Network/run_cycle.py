"""
Run the Mark 2 cycle balance.

    python run_cycle.py sizing      solve orifice areas to hit the design point
    python run_cycle.py offdesign   fix areas, solve the operating point
    python run_cycle.py smoke       synthetic design point, no CEA, checks wiring

The smoke path exists because the wiring should be debuggable without a CEA
install and without waiting on the forward solve.  It uses invented numbers and
proves nothing physical.
"""

import sys

from mark2_cycle import (
    AffinityPumpCurve,
    AnchoredRegenJacket,
    ChamberDesign,
    DesignPoint,
    GasGeneratorDesign,
    OrificeAreas,
    PumpDesign,
    RegenDesign,
    TurbineDesign,
    build_network,
    frozen_chamber,
    frozen_gas_generator,
)


# Placeholder until turbines.py stores its design velocity ratio; see
# ImpulseCurveEfficiency._design_velocity_ratio.
class ConstantTurbineEfficiency:
    validated = False

    def __init__(self, efficiency: float):
        self.efficiency = efficiency

    def __call__(self, velocity_ratio: float) -> float:
        return self.efficiency


def build(design: DesignPoint, *, mode: str, areas: OrificeAreas,
          chamber_temperature: float, chamber_gas_constant: float,
          coolant_specific_heat: float, link_efficiency: float,
          link_pullout_power: float):
    return build_network(
        design,
        fuel_pump_curve=AffinityPumpCurve(design.fuel_pump),
        ox_pump_curve=AffinityPumpCurve(design.ox_pump),
        regen_jacket=AnchoredRegenJacket(design.regen, specific_heat=coolant_specific_heat),
        turbine_efficiency=ConstantTurbineEfficiency(design.turbine.power_efficiency),
        chamber_gas=frozen_chamber(
            design.chamber,
            gas_constant=chamber_gas_constant,
            temperature=chamber_temperature,
        ),
        gas_generator_gas=frozen_gas_generator(design.gas_generator),
        link_efficiency=link_efficiency,
        link_pullout_power=link_pullout_power,
        areas=areas,
        mode=mode,
    )


def report(model):
    if model.unvalidated:
        print("\nPLACEHOLDER MODELS IN THIS RUN:")
        for entry in model.unvalidated:
            print(f"  - {entry}")
        print("  Results are a topology and solver check, not a cycle balance.\n")

    S = model.states
    C = model.components
    print(f"  chamber pressure      {S['pc'].value / 1e5:8.2f} bar")
    print(f"  gg pressure           {S['p_gg'].value / 1e5:8.2f} bar")
    print(f"  fuel pump discharge   {S['p_A'].value / 1e5:8.2f} bar")
    print(f"  ox pump discharge     {S['p_D'].value / 1e5:8.2f} bar")
    print(f"  fuel shaft speed      {S['N_fuel'].value:8.0f} rpm")
    print(f"  ox shaft speed        {S['N_ox'].value:8.0f} rpm")
    print(f"  fuel pump flow        {S['mdot_fuel_pump'].value:8.3f} kg/s")
    print(f"  ox pump flow          {S['mdot_ox_pump'].value:8.3f} kg/s")
    print(f"  chamber mixture ratio {S['mr_chamber'].value:8.3f}")
    print(f"  gg mixture ratio      {S['mr_gg'].value:8.3f}")
    print(f"  turbine power         {C['turbine'].shaft_power.value / 1e3:8.1f} kW")
    print(f"  link power            {C['link'].transmitted_power.value / 1e3:8.1f} kW")
    print(f"  pull-out margin       {C['link'].pullout_margin.value:8.2f}")

    if model.mode == "sizing":
        for key, label in (
            ("A_mfv", "main fuel valve"),
            ("A_tca_fuel", "TCA fuel injector"),
            ("A_tca_ox", "TCA ox injector"),
            ("A_gg_fuel", "GG fuel injector"),
            ("A_gg_ox", "GG ox injector"),
        ):
            print(f"  {label:22s}{S[key].value * 1e6:8.2f} mm2")


def synthetic_design_point() -> DesignPoint:
    """Invented but self-consistent numbers, for wiring checks only."""
    return DesignPoint(
        ambient_pressure=93900.0,
        ambient_temperature=293.0,
        fuel_tank_pressure=1.0e5,
        ox_tank_pressure=2.5e5,
        fuel_tank_temperature=293.0,
        ox_tank_temperature=90.2,
        fuel_name="JetA",
        ox_name="LOX",
        chamber=ChamberDesign(
            pressure=50.0e5,
            mixture_ratio=2.0,
            mass_flow=14.0,
            fuel_mass_flow=4.667,
            ox_mass_flow=9.333,
            throat_area=0.0045,
            expansion_ratio=4.0,
            c_star=1550.0,
            c_star_efficiency=0.85,
            gamma=1.20,
            injector_manifold_pressure=62.5e5,
        ),
        gas_generator=GasGeneratorDesign(
            pressure=40.0e5,
            mixture_ratio=0.2,
            mass_flow=0.75,
            fuel_mass_flow=0.625,
            ox_mass_flow=0.125,
            temperature=950.0,
            gas_constant=380.0,
            gamma=1.25,
            throat_area=2.625e-4,
            c_star=1400.0,
            injector_pressure=48.0e5,
        ),
        turbine=TurbineDesign(
            nozzle_plate_throat_area=2.625e-4,
            pitchline_radius=0.1143,
            shaft_speed=3142.0,
            power_efficiency=0.35,
            admission_fraction=0.0556,
            blade_length=0.00635,
            exit_pressure=93_900.0,
        ),
        regen=RegenDesign(
            mass_flow=4.667,
            pressure_drop=8.0e5,
            temperature_rise=110.0,
            heat_load=1.0e6,
            flow_direction=-1,
        ),
        fuel_pump=PumpDesign(
            fluid="fuel",
            inlet_pressure=1.0e5,
            discharge_pressure=72.0e5,
            mass_flow=5.292,
            density=800.0,
            head_rise=904.0,
            shaft_speed=3142.0,
            hydraulic_efficiency=0.60,
            shaft_power=97_700.0,
        ),
        ox_pump=PumpDesign(
            fluid="ox",
            inlet_pressure=2.5e5,
            discharge_pressure=63.2e5,
            mass_flow=9.458,
            density=1140.0,
            head_rise=542.0,
            shaft_speed=3142.0,
            hydraulic_efficiency=0.65,
            shaft_power=77_300.0,
        ),
        pole_pairs_fuel=1,
        pole_pairs_ox=1,
    )


def main(mode: str):
    if mode == "smoke":
        design = synthetic_design_point()
        chamber_temperature, chamber_gas_constant = 3500.0, 350.0
        coolant_specific_heat, link_efficiency = 2200.0, 0.87
        link_pullout_power = 4.0 * design.ox_pump.shaft_power
        solve_mode = "sizing"
    else:
        from mark2_cycle import capture

        design = capture()
        solve_mode = mode
        raise NotImplementedError(
            "Supply chamber stagnation temperature, chamber gas constant, coolant "
            "specific heat, link efficiency and link pull-out power from the forward "
            "solve.  See notes: these are the five numbers the cycle model needs that "
            "DesignPoint does not currently carry."
        )

    areas = OrificeAreas(
        main_fuel_valve=3.0e-4,
        tca_fuel_injector=1.2e-4,
        tca_ox_injector=1.8e-4,
        gg_fuel_injector=6.0e-5,
        gg_ox_injector=1.5e-5,
    )

    model = build(
        design,
        mode=solve_mode,
        areas=areas,
        chamber_temperature=chamber_temperature,
        chamber_gas_constant=chamber_gas_constant,
        coolant_specific_heat=coolant_specific_heat,
        link_efficiency=link_efficiency,
        link_pullout_power=link_pullout_power,
    )

    model.solve(verbose=True)
    report(model)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "smoke")
