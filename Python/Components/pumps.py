from __future__ import annotations
from dataclasses import dataclass, field
from engine_state import engine, tca, gg, pump, turbine
from typing import List
import numpy as np

from impellers import impellers
from blades import blades
from volute import volute
from inducers import inducer

def pumps(p: pump):

    p.head_rise = (p.p_out - p.p_in) / (p.density * engine.g) # m
    p.impeller[0].head_rise = p.head_rise  # assign head rise to impeller object - SHOULD SHARE WITH INDUCER
    p.vdot = p.mdot / p.density # m3/s - Volumetric Flow Rate

    p.n_pts_meridional = 100
    inlet_idx = int(0.1*p.n_pts_meridional) # THIS IS ARBITRARY
    p.impeller[0].blockage = 1.0*np.ones(p.n_pts_meridional)   # unitless - initial guess - proportional of meridional flow area physically blocked by blade material
    p.impeller[0].blockage[0:inlet_idx,] = 1 # update blockage vector to reflect blades only along the bladed portion of the meanline
    old_blockage = 0.0
    
    p.impeller[0].slip_factor = 0*np.ones_like(p.impeller[0].blockage) # unitless - initial guess - proportional reduction in tangential velocity from ideal velocity triangle

    p.impeller[0].hydraulic_efficiency = 1.0 - 0.071 / (p.vdot ** 0.25) # unitless - initial guess, Jekat's Empirical Formula; CONVERT TO ANDERSON?

    # the initial run is just to find shaft speed
    if not p.shaft_speed:
        inducer(p)
        return

    # the real runs
    inducer(p)
    # while abs(np.mean(p.impeller[0].blockage) - old_blockage) > 0.01:
    for iteration in range(1):
        # TODO: add functionality for a multistage pump (duplicate impellers, blades; add crossover module)
        impellers(p)
        blades(p)
    
    volute(p)

    #######################################################################################################
    ### Turbine Requirements
    #######################################################################################################
    
    p.head_coeff = p.impeller[0].head_coeff # unitless
    p.inlet_loss = 0.0
    p.outlet_loss = 0.0
    p.hydraulic_efficiency = p.head_coeff / (p.impeller[0].total_loss + p.volute.total_loss+ p.inlet_loss + p.outlet_loss) # Gulich 3.8.1 (the subscripts are German)
    
    p.shaft_power = p.vdot * p.head_rise * p.density * engine.g / p.hydraulic_efficiency # W
    p.shaft_power_theory = p.vdot * p.head_rise * p.density * engine.g # W
