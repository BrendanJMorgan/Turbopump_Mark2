from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from scipy.special import fresnel
from engine_state import engine, pump
from typing import List

def impellers(p: pump):
    # Impeller Contours - pg 2.29-2.32 in pump handbook
    r_eye_margin = 1.15 + 0.2116*(p.specific_speed-0.2836); 
        # unitless - margin correction factor, Gulich 7.1.4; interpolated for a normal impeller (not a suction impeller)
    swirl_number = 1; # NEEDS TO BE A FUNCTION OF INDUCER OUTLET CONDITIONS (ITERATIVE)
    p.r_inlet_impeller = r_eye_margin * np.sqrt(p.r_hub_impeller**2 + 0.0026418 * p.outlet_flow_coeff* p.specific_speed**1.33 / swirl_number**0.67) 
        # m - impeller inlet radius, Gulich eqn 7.1.4

    p.r_outlet_impeller = 1 / abs(p.shaft_speed) * np.sqrt (engine.g*p.head_rise / p.head_coeff); # m - impeller exit radius, Gulich eqn 7.1.3

    # The exit width b2 ------------------------------------------------------------------------------
    w_outlet = 2*p.r_outlet_impeller * (0.017 + 0.1386*p.specific_speed - 0.022387*p.specific_speed**2 + 0.0013767*p.specific_speed**3) 
        # m - outlet width; empirical; Gulich eqn 7.1

    # this needs boundary layer and leakage effect - pg 2.60?
    # could parameterize a width at each station with a blockage vector, linear ramp of area?

    # Hub and shroud profiles ------------------------------------------------------------------
    n_pts_meridional = 1000

    eps_ds = 10*np.pi/180 # rad - outlet shroud angle from Gulich table 7.1
    X = lambda t: np.sqrt(np.pi) * fresnel(t/np.sqrt(np.pi))[0] # clothoid / fresnel integral
    Y = lambda t: np.sqrt(np.pi) * fresnel(t/np.sqrt(np.pi))[1] # clothoid / fresnel integral
    T = np.sqrt(np.pi - 2*eps_ds) # upper limit of clothoid to achieve desired outlet angle

    clothoid = lambda t: ((p.r_inlet_impeller-p.r_outlet_impeller)/(np.sin(eps_ds)*Y(T) - np.cos(eps_ds)*X(T)))*np.vstack([X(t)*np.cos(eps_ds) - Y(t)*np.sin(eps_ds),
                                                                                                                X(t)*np.sin(eps_ds) + Y(t)*np.cos(eps_ds)]).T
    shroud_curve_raw = clothoid(np.linspace(T, 0, n_pts_meridional))    # [m,m] - (r,z) coordinates of shroud curve from outlet to inlet
    #shroud_curve_raw = np.flip(shroud_curve_raw, axis=1)                # [m,m] - reverse to go from inlet to outlet
    p.shroud_curve = (shroud_curve_raw + np.ones_like(shroud_curve_raw)*[-p.r_outlet_impeller, shroud_curve_raw[-1,1]])*[-1,1]  # [m,m] - translate the origin to the inlet central point, and flip from quadrant IV to I
    
    ds = np.linalg.norm(np.diff(p.shroud_curve, axis=0), axis=1)    # m - differential arc length along shroud curve
    s  = np.concatenate(([0.0], np.cumsum(ds)))                     # m - arc length along shroud curve
    s_shroud = s / s[-1]                                            # unitless - normalized arc length along shroud curve
    s_shroud = s_shroud.T                                           # make into column vector

    A_inlet = np.pi*(p.r_inlet_impeller**2 - p.r_hub_impeller**2)      # m2 - annular flow area at impeller inlet
    A_outlet = 2*np.pi*p.r_outlet_impeller*w_outlet                       # m2 - annular flow area at impeller outlet
    weight = 2*np.pi*p.r_outlet_impeller/np.tan(eps_ds) - np.pi*w_outlet  # unitless - weighting factor for area distribution along meanline
    p.A_meanline = (A_outlet-A_inlet) * (2*(weight-1)*s_shroud**3 - 3*(weight-1)*s_shroud**2 + weight*s_shroud) + A_inlet # m2 - area distribution along meanline

    tangents_shroud = np.empty_like(p.shroud_curve)
    tangents_shroud[1:-1] = p.shroud_curve[2:] - p.shroud_curve[:-2]
    tangents_shroud[0]    = p.shroud_curve[1] - p.shroud_curve[0]
    tangents_shroud[-1]   = p.shroud_curve[-1] - p.shroud_curve[-2]
    normals_shroud = np.column_stack([-tangents_shroud[:,1], tangents_shroud[:,0]]) # rotate tangent +90° in (r,z) plane -> normal
    normals_shroud /= np.linalg.norm(normals_shroud, axis=1, keepdims=True)

    determinant = p.shroud_curve[:,0]**2 - p.A_meanline*p.shroud_curve[:,0]/np.pi
    crosswise_meanline = (-p.shroud_curve[:,0] - np.sqrt(p.shroud_curve[:,0]**2 - p.A_meanline*normals_shroud[:,0]/np.pi)) / normals_shroud[:,0] # m - crosswise distance between shroud and hub curves
    p.hub_curve = p.shroud_curve + crosswise_meanline[:, None] * normals_shroud # m - (r,z) coordinates of hub curve from outlet to inlet
    
    # dummy = p.shroud_curve[:,0]**2 - p.A_meanline*normals_shroud[:,0]/np.pi
    # d1 = p.shroud_curve[999,0]
    # d2 = p.A_meanline[999]
    # d3 = normals_shroud[999,0]
    # d4 = shroud_curve_raw[999,0]
    p.meanline_curve = 0.5 * (p.shroud_curve + p.hub_curve); # [m,m] - halfway between shroud and impeller 





    