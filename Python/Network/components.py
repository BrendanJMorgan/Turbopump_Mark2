"""
Custom FullFlow components for Mark 2.

FullFlow's catalog covers the cycle except for three small pieces: the
synchronous electrical link between the two shafts, the choked-nozzle-plate flow
parameter that feeds the turbine, and a partial-admission impulse stage.
"""

import math

from fullflow import Component, Network, State


class SynchronousLink(Component):
    """
    Pole-locked electrical coupling between the fuel-side alternator and the
    LOX-side motor.

    Steady state, both machines locked:

        omega_ox = omega_fuel * pole_pair_ratio

    so the two rotating assemblies are one shaft-speed degree of freedom.  The
    LOX pump load is reflected back onto the fuel shaft as an alternator
    reaction torque through the round-trip efficiency:

        P_transmitted     = torque_ox * omega_ox
        drive_shaft_torque = P_transmitted / (eta_roundtrip * omega_fuel)

    with ``eta_roundtrip = alternator * cable * motor``.  The fuel-side Rotor
    then closes ``torque_turbine - torque_fuel_pump - drive_shaft_torque = 0``
    on rotor speed.

    Power-balance model only: no torque-angle relation, so it cannot predict
    pull-out.  Check the margin separately against ``synchronous_links.py``.
    """

    def __init__(
        self,
        name: str,
        network: Network,
        drive_shaft_speed_rpm: State,        # fuel shaft, rpm
        driven_shaft_speed_rpm: State,       # LOX shaft, rpm (kinematically = drive * ratio)
        driven_shaft_torque: State,          # LOX pump load, N*m
        alternator_efficiency: State | float,
        cable_efficiency: State | float,
        motor_efficiency: State | float,
        drive_shaft_torque: State | None = None,
        transmitted_power: State | None = None,
        electrical_power: State | None = None,
        electrical_frequency: State | None = None,
    ):
        self.setup()

    def evaluate_states(self):
        drive_speed = self.drive_shaft_speed_rpm.value
        if drive_speed <= 0.0:
            raise ValueError(
                f"{self.name}: drive shaft speed {drive_speed} rpm is non-positive; "
                "the synchronous link has no meaning at or below zero speed"
            )

        efficiency = (
            self.alternator_efficiency.value
            * self.cable_efficiency.value
            * self.motor_efficiency.value
        )

        omega_drive = drive_speed * math.pi / 30.0
        omega_driven = self.driven_shaft_speed_rpm.value * math.pi / 30.0

        transmitted_power = self.driven_shaft_torque.value * omega_driven

        self.transmitted_power.value = transmitted_power
        self.drive_shaft_torque.value = transmitted_power / (efficiency * omega_drive)
        # power drawn from the alternator terminals (before cable + motor losses)
        self.electrical_power.value = transmitted_power / (
            self.cable_efficiency.value * self.motor_efficiency.value
        )
        self.electrical_frequency.value = drive_speed / 60.0


class ChokedFlowParameter(Component):
    """
    Turbine-nozzle-plate flow parameter for a choked orifice.

    ``GasTurbine`` computes ``mdot = FP * Pt / sqrt(R * Tt)``.  For a choked
    plate of total area ``A`` and discharge coefficient ``Cd``:

        FP = Cd * A * sqrt(gamma) * (2 / (gamma + 1)) ** ((gamma + 1) / (2 (gamma - 1)))

    (the gas constant cancels against the ``sqrt(R)`` inside ``GasTurbine``).
    """

    def __init__(
        self,
        name: str,
        network: Network,
        throat_area: State | float,
        discharge_coefficient: State | float,
        specific_heat_ratio: State | float,
        flow_parameter: State | None = None,
    ):
        self.setup()

    def evaluate_states(self):
        gamma = self.specific_heat_ratio.value
        exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
        self.flow_parameter.value = (
            self.discharge_coefficient.value
            * self.throat_area.value
            * math.sqrt(gamma)
            * (2.0 / (gamma + 1.0)) ** exponent
        )


class ImpulseTurbineStage(Component):
    """
    Single-row partial-admission impulse stage, fed by an upstream choked
    nozzle plate that sets the mass flow.

    Explicit component, no closure equation of its own.  Mass flow comes from
    the nozzle-plate orifice; this converts it to torque:

        dh_ideal = cp * Tt * (1 - (p_exit / Pt)^((g-1)/g))
        P_shaft  = efficiency * mdot * dh_ideal
        torque   = P_shaft / omega

    ``efficiency`` is the total-to-static power efficiency including partial
    admission, i.e. turbines.py's ``power_efficiency``.
    """

    def __init__(
        self,
        name: str,
        network: Network,
        mass_flow: State,
        rotor_speed: State,                     # rpm
        upstream_total_pressure: State,
        upstream_total_temperature: State,
        downstream_pressure: State,
        gas_constant: State | float,
        specific_heat_ratio: State | float,
        efficiency: State,
        torque: State | None = None,
        shaft_power: State | None = None,
        ideal_enthalpy_change: State | None = None,
        spouting_velocity: State | None = None,
        discharge_total_enthalpy: State | None = None,
        upstream_total_enthalpy: State | None = None,
    ):
        self.setup()

    def evaluate_states(self):
        gamma = self.specific_heat_ratio.value
        gas_constant = self.gas_constant.value
        total_pressure = self.upstream_total_pressure.value
        total_temperature = self.upstream_total_temperature.value
        exit_pressure = self.downstream_pressure.value
        mass_flow = self.mass_flow.value
        speed = self.rotor_speed.value

        pressure_ratio = exit_pressure / total_pressure
        if pressure_ratio >= 1.0:
            raise ValueError(
                f"{self.name}: turbine exit pressure {exit_pressure:.0f} Pa is at or above "
                f"inlet total pressure {total_pressure:.0f} Pa"
            )

        specific_heat = gamma * gas_constant / (gamma - 1.0)
        ideal_enthalpy_change = specific_heat * total_temperature * (
            1.0 - pressure_ratio ** ((gamma - 1.0) / gamma)
        )
        shaft_power = self.efficiency.value * mass_flow * ideal_enthalpy_change
        omega = speed * math.pi / 30.0

        self.ideal_enthalpy_change.value = ideal_enthalpy_change
        self.spouting_velocity.value = math.sqrt(2.0 * ideal_enthalpy_change)
        self.shaft_power.value = shaft_power
        self.torque.value = shaft_power / omega

        if self.upstream_total_enthalpy.is_assigned and abs(mass_flow) > 1e-12:
            self.discharge_total_enthalpy.value = (
                self.upstream_total_enthalpy.value - shaft_power / mass_flow
            )
