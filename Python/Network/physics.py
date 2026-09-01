"""
Named physics hooks for the Mark 2 cycle assembler.

``cycle.yaml`` owns the wiring; these functions own the equations that do not
map cleanly onto a single component - pump torque from the map, pump dissipation
heating, mixture ratios, the CEA equilibrium chain, the turbine efficiency
closure, the synchronous link, and engine performance.

Each hook takes the shared ``model.BuildContext`` and mutates it in place,
publishing results into ``ctx.derived`` (referenced from the YAML as ``d:<name>``),
``ctx.temperature``, ``ctx.performance`` and ``ctx.components``.  They run in the
order given by the YAML ``build`` manifest.
"""

from __future__ import annotations

import math

from fullflow import Balance, Lookup
from thermoprop import Equilibrium, Fluid, Propellant, Reactants

_RPM_TO_RAD = math.pi / 30.0


# ----------------------------------------------------------------------------
# shafts
# ----------------------------------------------------------------------------
def shaft_speeds(ctx) -> None:
    """One shaft-speed unknown; the pole lock fixes the other."""
    fuel_speed = ctx.scalar["fuel_shaft_speed_rpm"]
    ox_speed = fuel_speed * ctx.params["pole_pair_ratio"]

    ctx.derived["oxidizer_shaft_speed_rpm"] = ox_speed
    ctx.derived["fuel_omega"] = fuel_speed * _RPM_TO_RAD
    ctx.derived["oxidizer_omega"] = ox_speed * _RPM_TO_RAD
    ctx.track("shaft/oxidizer_speed", ox_speed)


# ----------------------------------------------------------------------------
# oxidizer side
# ----------------------------------------------------------------------------
def oxidizer_inlet(ctx) -> None:
    """Real-fluid LOX state at the pump inlet (density for the whole ox side)."""
    inlet = Lookup(
        "LOX Pump Inlet",
        ctx.net,
        Fluid,
        ctx.propellants.oxidizer_fluid,
        pressure=ctx.pressure["ox_pump_in"],
        temperature=ctx.temperature["ox_tank"],
    )
    ctx.derived["oxidizer_density"] = inlet.density


def oxidizer_pump(ctx) -> None:
    """Ox pump map + shaft torque formed as g*H*mdot/(eta*omega)."""
    density = ctx.derived["oxidizer_density"]
    mass_flow = ctx.mass_flow["ox_pump"]

    ox_map = ctx.load_map(
        "oxidizer_pump",
        inputs={
            "volumetric_flow": mass_flow / density,
            "rotor_speed": ctx.derived["oxidizer_shaft_speed_rpm"],
        },
        outputs=["head_rise", "efficiency"],
    )
    ctx.derived["ox_map"] = ox_map

    torque = (ctx.params["gravity"] * ox_map.head_rise * mass_flow) / (
        ox_map.efficiency * ctx.derived["oxidizer_omega"]
    )
    ctx.derived["ox_pump_torque"] = torque
    ctx.components["ox_pump_torque"] = torque


# ----------------------------------------------------------------------------
# fuel side
# ----------------------------------------------------------------------------
def fuel_inlet(ctx) -> None:
    inlet = Lookup(
        "Fuel Pump Inlet",
        ctx.net,
        Propellant,
        ctx.propellants.fuel_propellant,
        pressure=ctx.pressure["fuel_pump_in"],
        temperature=ctx.temperature["fuel_tank"],
    )
    ctx.derived["fuel_density"] = inlet.density
    ctx.derived["fuel_specific_heat"] = inlet.specific_heat_cp


def fuel_pump(ctx) -> None:
    """Fuel pump map, shaft torque, and the pump-dissipation temperature rise."""
    gravity = ctx.params["gravity"]
    density = ctx.derived["fuel_density"]
    specific_heat = ctx.derived["fuel_specific_heat"]
    mass_flow = ctx.mass_flow["fuel_pump"]

    fuel_map = ctx.load_map(
        "fuel_pump",
        inputs={
            "volumetric_flow": mass_flow / density,
            "rotor_speed": ctx.scalar["fuel_shaft_speed_rpm"],
        },
        outputs=["head_rise", "efficiency"],
    )
    ctx.derived["fuel_map"] = fuel_map

    torque = (gravity * fuel_map.head_rise * mass_flow) / (
        fuel_map.efficiency * ctx.derived["fuel_omega"]
    )
    ctx.derived["fuel_pump_torque"] = torque
    ctx.components["fuel_pump_torque"] = torque

    # dT = (1/eta - 1) * g * H / cp
    temperature_rise = (
        (1.0 / fuel_map.efficiency - 1.0) * gravity * fuel_map.head_rise / specific_heat
    )
    ctx.set_temperature(
        "fuel_pump_out", ctx.temperature["fuel_tank"] + temperature_rise
    )


def regen(ctx) -> None:
    """Lumped jacket heat pickup: one heat load, one dT."""
    outlet_temperature = ctx.temperature["fuel_pump_out"] + ctx.params[
        "regen_heat_load"
    ] / (ctx.mass_flow["regen"] * ctx.derived["fuel_specific_heat"])
    ctx.set_temperature("jacket_out", outlet_temperature)


# ----------------------------------------------------------------------------
# combustion (CEA equilibrium chain for TCA and GG)
# ----------------------------------------------------------------------------
def combustion(ctx) -> None:
    pressure = ctx.pressure
    mass_flow = ctx.mass_flow
    temperature = ctx.temperature
    net = ctx.net
    propellants = ctx.propellants

    # Film coolant is carried through the nozzle mass balance but excluded from
    # the core mixture ratio.  No mixing model: chamber gas properties are the
    # core values while the throat passes the total flow.
    core_fuel_flow = mass_flow["tca_fuel"] * (1.0 - ctx.params["film_coolant_fraction"])
    tca_mixture_ratio = mass_flow["tca_ox"] / core_fuel_flow
    gg_mixture_ratio = mass_flow["gg_ox"] / mass_flow["gg_fuel"]

    ctx.set_performance("tca_mixture_ratio", tca_mixture_ratio)
    ctx.set_performance("gg_mixture_ratio", gg_mixture_ratio)

    ox_reactant = Lookup(
        "Ox Reactant", net, Propellant,
        propellants.oxidizer_reactant,
        temperature=temperature["ox_tank"],
    )
    tca_fuel_reactant = Lookup(
        "TCA Fuel Reactant", net, Propellant,
        propellants.fuel_reactant,
        temperature=temperature["jacket_out"],
    )
    gg_fuel_reactant = Lookup(
        "GG Fuel Reactant", net, Propellant,
        propellants.fuel_reactant,
        temperature=temperature["fuel_pump_out"],
    )

    tca_reactants = Lookup(
        "TCA Reactants", net, Reactants,
        fuels=tca_fuel_reactant,
        oxidizers=ox_reactant,
        mixture_ratio=tca_mixture_ratio,
    )
    gg_reactants = Lookup(
        "GG Reactants", net, Reactants,
        fuels=gg_fuel_reactant,
        oxidizers=ox_reactant,
        mixture_ratio=gg_mixture_ratio,
    )

    tca_combustion = Lookup(
        "TCA Combustion", net, Equilibrium,
        reactants=tca_reactants, mode="hp", pressure=pressure["chamber"],
    )
    gg_combustion = Lookup(
        "GG Combustion", net, Equilibrium,
        reactants=gg_reactants, mode="hp", pressure=pressure["gas_generator"],
    )
    ctx.components["tca_combustion"] = tca_combustion
    ctx.components["gg_combustion"] = gg_combustion
    ctx.derived["tca_combustion"] = tca_combustion
    ctx.derived["gg_combustion"] = gg_combustion

    # c* efficiency applied as eta^2 on total temperature, since c* ~ sqrt(R*T0).
    tca_total_temperature = tca_combustion.temperature * ctx.params["c_star_efficiency"] ** 2
    gg_total_temperature = gg_combustion.temperature * ctx.params["gg_c_star_efficiency"] ** 2
    ctx.derived["tca_total_temperature"] = tca_total_temperature
    ctx.derived["gg_total_temperature"] = gg_total_temperature

    ctx.set_temperature("chamber", tca_total_temperature)
    ctx.set_temperature("gas_generator", gg_total_temperature)


# ----------------------------------------------------------------------------
# turbine closure and shaft balance
# ----------------------------------------------------------------------------
def turbine_balance(ctx) -> None:
    """GasTurbine takes torque as an input and reports efficiency; close the loop."""
    turbine = ctx.components["turbine"]
    Balance(
        "Turbine Efficiency",
        ctx.net,
        variable=ctx.scalar["turbine_torque"],
        function=turbine.efficiency - ctx.params["turbine_efficiency"],
    )


def shaft_balance(ctx) -> None:
    """Fuel-shaft net torque for the Rotor: turbine drive minus the two loads.

    The ``link`` element (SynchronousLink) is built from the YAML just before
    this hook and has already published ``d:link_reaction_torque`` and
    ``d:electrical_power`` through its ``exports``.
    """
    net_torque = (
        ctx.scalar["turbine_torque"]
        - ctx.derived["fuel_pump_torque"]
        - ctx.derived["link_reaction_torque"]
    )
    ctx.derived["net_torque"] = net_torque

    ctx.components["turbine_torque"] = ctx.scalar["turbine_torque"]
    ctx.components["link_reaction_torque"] = ctx.derived["link_reaction_torque"]
    ctx.components["electrical_power"] = ctx.derived["electrical_power"]


# ----------------------------------------------------------------------------
# performance
# ----------------------------------------------------------------------------
def performance(ctx) -> None:
    params = ctx.params
    mass_flow = ctx.mass_flow
    nozzle = ctx.components["nozzle"]
    ambient_pressure = ctx.derived["ambient_pressure"]

    exit_area = params["throat_area"] * params["expansion_ratio"]
    momentum_thrust = nozzle.mass_flow * nozzle.exit_velocity
    pressure_thrust = (nozzle.exit_static_pressure - ambient_pressure) * exit_area
    chamber_thrust = (momentum_thrust + pressure_thrust) * params["thrust_coefficient_efficiency"]

    total_mass_flow = mass_flow["nozzle"] + mass_flow["turbine"]

    ctx.set_performance("chamber_thrust", chamber_thrust)
    ctx.set_performance("total_mass_flow", total_mass_flow)
    ctx.set_performance("specific_impulse", chamber_thrust / (total_mass_flow * params["gravity"]))
    ctx.set_performance(
        "engine_mixture_ratio",
        (mass_flow["tca_ox"] + mass_flow["gg_ox"])
        / (mass_flow["tca_fuel"] + mass_flow["gg_fuel"]),
    )
    # NOTE: GG exhaust thrust is not included - the dump nozzle is not modelled.
