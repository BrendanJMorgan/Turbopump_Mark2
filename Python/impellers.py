from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from engine_state import engine, pump
from typing import List
import math

def impellers(p: pump):
    # Impeller Contours - pg 2.29-2.32 in pump handbook
    r_eye_margin = 1.15 + 0.2116*(p.specific_speed-0.2836); 
        # unitless - margin correction factor, Gulich 7.1.4; interpolated for a normal impeller (not a suction impeller)
    swirl_number = 1; # NEEDS TO BE A FUNCTION OF INDUCER OUTLET CONDITIONS (ITERATIVE)
    r_inlet_impeller = r_eye_margin * np.sqrt(p.r_hub**2 + 0.0026418 * p.outlet_flow_coeff* p.specific_speed**1.33 / swirl_number**0.67) 
        # m - impeller inlet radius, Gulich eqn 7.1.4

    r_outlet_impeller = 1 / abs(p.shaft_speed) * np.sqrt (engine.g*p.head_rise / p.head_coeff); # m - impeller exit radius, Gulich eqn 7.1.3

    # The exit width b2 ------------------------------------------------------------------------------
    w_exit = 2*r_outlet_impeller * (0.017 + 0.1386*p.specific_speed - 0.022387*p.specific_speed**2 + 0.0013767*p.specific_speed**3) 
        # m - outlet width; empirical; Gulich eqn 7.1

    du = 0.001 # unitless - differential arc length, normalized
    u = np.linspace(0, 1, np.round(1/ds)) # m - normalized arc length along meanline from inlet (0) to outlet (1)
    weight = 0.1 # unitless - weighting factor for area distribution shape

    # this needs boundary layer and leakage effect - pg 2.60?
    # could parameterize a width at each station with a blockage vector, linear ramp of area?

    # Hub and shroud profiles ------------------------------------------------------------------
    # specify a smoothly varying curvature and derive a curve out of that?
    cavitation_curvature = 2/r_inlet_impeller # pump handbook recommends a minimum turning radius of hlaf the eye radius
    kappa = np.linspace(cavitation_curvature, 0, len(u)) # 1/m - curvature of the shroud curve; flattens out at the end
        # this should be modified to have a short inlet section with low curvature, to achieve a flatter pressure distribution impinging on blades

    inlet_azimuth_shroud = -np.pi/2 # rad
    outlet_azimuth_shroud = 0 # rad - replace with epsilon_ds from Gulich

    length_shroud = (outlet_azimuth_shroud - inlet_azimuth_shroud) / np.trapezoid(kappa, u) # m - total arc length of the shroud

    theta = np.empty(len(u)) # rad
    theta[0] = inlet_azimuth_shroud # rad

    for i in range(1, len(u)):
        k_avg = 0.5 * (kappa[i-1] + kappa[i])
        theta[i] = theta[i-1] + length_shroud*k_avg*ds

    r_shroud = np.empty(len(u))
    z_shroud = np.empty(len(u))
    r_shroud[0] = r_inlet_impeller
    z_shroud[0] = 0 # puts the curve into quadrant IV of the coordinate plane

    ds_shroud = length_shroud * du # m - differential arc length
    s_shroud = length_shroud * s # m - arc length
    for i in range(1, len(u)):
        cos_avg = 0.5 * (np.cos(theta[i - 1]) + np.cos(theta[i]))
        sin_avg = 0.5 * (np.sin(theta[i - 1]) + np.sin(theta[i]))
        r_shroud[i] = r_shroud[i - 1] + cos_avg * ds_shroud
        z_shroud[i] = z_shroud - 1 + sin_avg * ds_shroud

    A_inlet = np.pi*(r_outlet_impeller**2 - r_inlet_impeller**2) # m2 - annular flow area at impeller inlet
    A_outlet = 2*np.pi*r_outlet_impeller*w_exit # m2 - annular flow area at impeller outlet
    A = A_inlet-A_outlet * (2*(weight-1)*s_shroud**3 - 3*(weight-1)*s_shroud**2 + weight*s_shroud) + A_outlet # m2 - area distribution along meanline








    