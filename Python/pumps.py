from __future__ import annotations
from dataclasses import dataclass, field
from engine_state import engine, tca, gg, pump, turbine
from typing import List
import math

from impellers import impellers
from blades import blades
from volute import volute
from inducer import inducer

def pumps(p: pump):

    p.head_pump = (p.p_out - p.p_in) / (p.density_pump * math.g) # m

    p.vdot_pump = p.mdot_pump / p.density_pump # m3/s - Volumetric Flow Rate

    # Sizing the Pump --------------------------------------------------------
    p.specific_speed = abs(p.pump_shaft_speed) * math.sqrt(p.vdot_pump) / (p.head_pump ** 0.75) # unitless 

    p.specific_speed_imperial = p.specific_speed * 229.6        # Ns pseudo unitless
    p.outlet_flow_coeff = 0.1715 * math.sqrt(p.specific_speed)  # unitless - phi_i in pump handbook

    ft = 1.0; # unitless - could range to 1.1 to be more aggressive on head, on low head pumps
    p.shutoff_coeff_op = 0.605 * ft * math.exp(-0.408 * p.specific_speed)  # unitless - psi_op; eqn 3.26 in Gulich
    p.shutoff_coeff = 1.25; # 1.25 for volute designs and 1.31 for diffuser designs (kept as input)
    p.shutoff_head_coeff = 0.5 * p.shutoff_coeff * 0.625 * math.exp(-4.15 * p.specific_speed) # unitless - psi_0 in the limit pf Q=0

    # specific speed check - should be in the centrifugal range
    if p.specific_speed < 0.1 or p.specific_speed > 0.6:
        print(f"[{p.name}] Pump Specific Speed is {p.specific_speed:.4g}. "
              f"Recommended is 0.1 to 0.6 for centrifugal geometries.")
        
    p.blockage = 0.85       # unitless - iterative initial guess - 1 is completely open, pump handbook says 0.85 is typical?
    p.old_blockage = 0      # placeholder
    p.blade_count = 16      # unitless - iterative initial guess
    p.old_blade_count = 0   # placeholder
    p.slip_factor = 0.15;   # unitless - iterative initial guess

    # Subfunction Calls --------------------------------------------------------
    while abs(p.blockage - old_blockage) > 0.01:
        old_slip_factor = p.slip_factor
        old_blockage = p.blockage
        old_blade_count = p.blade_count

        # TODO: add functionality for a multistage pump (duplicate impellers, blades; add crossover module)
        impellers(p)
        blades(p)

    volute(p)
    inducer(p)

    # Requirements for Turbine ----------------------------------------------
    term1 = (p.vdot_pump / (p.surface_roughness**2 * abs(p.pump_shaft_speed))) ** -0.2133
    term2 = (math.log10(0.83644 / p.specific_speed)) ** 2
    p.pump_efficiency = 0.94 - 6.983e-4 * term1 - 0.29 * term2 # unitless - Anderson; pump handbook pg 2.25

    p.shaft_power = p.vdot_pump * p.head_pump * inputs.g / p.pump_efficiency # W
    p.shaft_power_theory = p.vdot_pump * p.head_pump * p.density_pump * math.g # W

    # Other
    p.v_inlet_pump = p.vdot_pump / ( (math.pi/4.0) * (p.d_inlet_pump ** 2) ) # m/s - inlet velocity
    p.cavitation_inlet = (p.p_in - p.vapor_pressure_pump) / (0.5 * p.density_pump * (p.v_inlet_pump ** 2)) # unitless
    p.v_outlet_pump = p.vdot_pump / ( (math.pi/4.0) * (p.d_outlet_pump ** 2) ) # m/s - outlet velocity
    p.cavitation_outlet = (p.p_out - p.vapor_pressure_pump) / (0.5 * p.density_pump * (p.v_outlet_pump ** 2)) # unitless
    p.v_volute_pump = math.sqrt(2.0 * p.power_pump / p.mdot_pump) # m/s - average tangential velocity along the contour of the volute

