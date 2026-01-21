from __future__ import annotations
from dataclasses import dataclass, field
from engine_state import engine, tca, gg, pump, turbine
from typing import List
import numpy as np

from impellers import impellers
from blades import blades
from volute import volute
from inducers import inducers

def pumps(p: pump):

    p.head_rise = (p.p_out - p.p_in) / (p.density * engine.g) # m
    p.impeller[0].head_rise = p.head_rise  # assign head rise to impeller object
    p.vdot = p.mdot / p.density # m3/s - Volumetric Flow Rate

    p.n_pts_meridional = 1000
    inlet_idx = int(0.1*p.n_pts_meridional)
    p.impeller[0].blockage = 0.85*np.ones(p.n_pts_meridional)   # unitless - initial guess - proportional of meridional flow area physically blocked by blade material
    p.impeller[0].blockage[0:inlet_idx,] = 1 # update blockage vector to reflect blades only along the bladed portion of the meanline
    old_blockage = 0.0
    
    p.impeller[0].slip_factor = 0.15; # unitless - initial guess - proportional reduction in tangential velocity from ideal velocity triangle

    p.impeller[0].hydraulic_efficiency = 1.0 - 0.071 * (p.vdot ** 0.25) # unitless - initial guess, Jekat's Empirical Formula; CONVERT TO ANDERSON?

    # the initial run is just to find shaft speed
    if not p.shaft_speed:
        inducers(p)
        return

    # the real runs
    inducers(p)
    while abs(np.mean(p.impeller[0].blockage) - old_blockage) > 0.01:
        old_blockage = np.mean(p.impeller[0].blockage)
        # TODO: add functionality for a multistage pump (duplicate impellers, blades; add crossover module)
        impellers(p)
        blades(p)
    
    # volute(p)
    # inducer(p)

    # Requirements for Turbine ----------------------------------------------
    term1 = (p.vdot / (p.surface_roughness**2 * abs(p.shaft_speed))) ** -0.2133
    term2 = (np.log10(0.83644 / p.impeller[0].specific_speed)) ** 2
    p.efficiency = 0.94 - 6.983e-4 * term1 - 0.29 * term2 # unitless - Anderson; pump handbook pg 2.25

    p.shaft_power = p.vdot * p.head_rise * engine.g / p.efficiency # W
    p.shaft_power_theory = p.vdot * p.head_rise * p.density * engine.g # W

    # Other
    # p.v_inlet = p.vdot / ( (np.pi/4.0) * (p.d_inlet_pump ** 2) ) # m/s - inlet velocity
    # p.cavitation_inlet = (p.p_in - p.pvap_inlet) / (0.5 * p.density * (p.v_inlet ** 2)) # unitless
    # p.v_outlet_pump = p.vdot / ( (np.pi/4.0) * (p.d_outlet_pump ** 2) ) # m/s - outlet velocity
    # p.cavitation_outlet = (p.p_out - p.pvap_inlet) / (0.5 * p.density * (p.v_outlet_pump ** 2)) # unitless
    # p.v_volute = np.sqrt(2.0 * p.power / p.mdot) # m/s - average tangential velocity along the contour of the volute

