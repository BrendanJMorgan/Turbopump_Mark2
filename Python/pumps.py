from __future__ import annotations
from dataclasses import dataclass, field
from engine_state import engine, tca, gg, pump, turbine
from typing import List
import numpy as np

from impellers import impellers
from blades import blades
from volute import volute
from inducer import inducer

def pumps(p: pump):

    p.head_rise = (p.p_out - p.p_in) / (p.density * engine.g) # m
    p.vdot = p.mdot / p.density # m3/s - Volumetric Flow Rate

    # Sizing the Pump --------------------------------------------------------
    p.specific_speed = abs(p.shaft_speed) * np.sqrt(p.vdot) / (engine.g * p.head_rise ** 0.75) # unitless 

    p.specific_speed_imperial = p.specific_speed * 229.6        # Ns pseudo unitless
    p.outlet_flow_coeff = 0.1715 * np.sqrt(p.specific_speed)  # unitless - phi_i in pump handbook

    ft = 1.0; # unitless - could range to 1.1 to be more aggressive on head, on low head pumps
    p.head_coeff = 0.5 * 1.21 * ft * np.exp(-0.408 * p.specific_speed)  # unitless - psi_op; eqn 3.26 in Gulich
    p.shutoff_coeff = 1.25; # 1.25 for volute designs and 1.31 for diffuser designs (kept as input)
    p.shutoff_head_coeff = 0.5 * p.shutoff_coeff * 0.625 * np.exp(-4.15 * p.specific_speed) # unitless - psi_0 in the limit pf Q=0

    # specific speed check - should be in the centrifugal range
    if p.specific_speed < 0.1 or p.specific_speed > 0.6:
        print(f"{p.fluid} pump specific speed is {round(p.specific_speed,2)}. Recommended is 0.1 to 0.6 for centrifugal geometries.")
        
    p.blockage = 0.85       # unitless - iterative initial guess - 1 is completely open, pump handbook says 0.85 is typical?
    p.old_blockage = 0      # placeholder
    p.blade_count = 16      # unitless - iterative initial guess
    p.old_blade_count = 0   # placeholder
    p.slip_factor = 0.15;   # unitless - iterative initial guess

    old_slip_factor = 0     # placeholder
    old_blockage = 0        # placeholder
    old_blade_count = 0     # placeholder

    # Subfunction Calls --------------------------------------------------------
    while abs(p.blockage - old_blockage) > 0.01:
        old_slip_factor = p.slip_factor
        old_blockage = p.blockage
        old_blade_count = p.blade_count

        # TODO: add functionality for a multistage pump (duplicate impellers, blades; add crossover module)
        impellers(p)
        blades(p)

    #volute(p)
    inducer(p)

    # Requirements for Turbine ----------------------------------------------
    term1 = (p.vdot / (p.surface_roughness**2 * abs(p.shaft_speed))) ** -0.2133
    term2 = (np.log10(0.83644 / p.specific_speed)) ** 2
    p.efficiency = 0.94 - 6.983e-4 * term1 - 0.29 * term2 # unitless - Anderson; pump handbook pg 2.25

    p.shaft_power = p.vdot * p.head_rise * engine.g / p.efficiency # W
    p.shaft_power_theory = p.vdot * p.head_rise * p.density * engine.g # W

    # Other
    # p.v_inlet = p.vdot / ( (np.pi/4.0) * (p.d_inlet_pump ** 2) ) # m/s - inlet velocity
    # p.cavitation_inlet = (p.p_in - p.pvap_inlet) / (0.5 * p.density * (p.v_inlet ** 2)) # unitless
    # p.v_outlet_pump = p.vdot / ( (np.pi/4.0) * (p.d_outlet_pump ** 2) ) # m/s - outlet velocity
    # p.cavitation_outlet = (p.p_out - p.pvap_inlet) / (0.5 * p.density * (p.v_outlet_pump ** 2)) # unitless
    # p.v_volute = np.sqrt(2.0 * p.power / p.mdot) # m/s - average tangential velocity along the contour of the volute

