from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from scipy.special import fresnel
import matplotlib.pyplot as plt

from engine_state import engine, pump
from typing import List

def impellers(p: pump):

    p.r_outlet_impeller = 1 / abs(p.shaft_speed) * np.sqrt (engine.g*p.head_rise / p.head_coeff); # m - impeller exit radius, Gulich eqn 7.1.3

    r_eye_margin = 1.15 + 0.2116*(p.specific_speed-0.2836); 
        # unitless - margin correction factor, Gulich 7.1.4; interpolated for a normal impeller (not a suction impeller)
    swirl_number = 0.5; # NEEDS TO BE ITERATIVE 
    p.r_inlet_impeller = p.r_outlet_impeller * r_eye_margin * np.sqrt((p.r_hub_impeller/p.r_outlet_impeller)**2 + 
                                                                    0.00148*2*p.head_coeff*(p.specific_speed*52.9186)**1.33/swirl_number**0.67) 
        # m - impeller inlet radius, Gulich eqn 7.1.4
    p.eye_flow_coeff = p.vdot / (np.pi*(p.r_inlet_impeller**2-p.r_hub_impeller**2)*p.r_inlet_impeller * p.shaft_speed) # unitless

    # The exit width b2 ------------------------------------------------------------------------------
    w_outlet = 2*p.r_outlet_impeller * (0.017 + 0.1386*p.specific_speed - 0.022387*p.specific_speed**2 + 0.0013767*p.specific_speed**3) 
        # m - outlet width; empirical; Gulich eqn 7.1

    # this needs boundary layer and leakage effect - pg 2.60?
    # could parameterize a width at each station with a blockage vector, linear ramp of area?

    # Hub and shroud profiles ------------------------------------------------------------------
    n_pts_meridional = 1000

    eps_ds = 10*np.pi/180 # rad - outlet shroud angle from Gulich table 7.1
    X = lambda t: fresnel(t/np.sqrt(np.pi))[1] # clothoid / fresnel integral
    Y = lambda t: fresnel(t/np.sqrt(np.pi))[0] # clothoid / fresnel integral
    T = np.sqrt(np.pi - 2*eps_ds) # upper limit of clothoid to achieve desired outlet angle

    clothoid = lambda t: (p.r_inlet_impeller-p.r_outlet_impeller)/(np.sin(eps_ds)*Y(T) - np.cos(eps_ds)*X(T)) * np.vstack([X(t), Y(t)]).T  # [m,m] - parametric clothoid function scaled to desired inlet and outlet radii
    c = clothoid(np.linspace(T,0,n_pts_meridional))    # [m,m] - (r,z) coordinates of shroud curve from outlet to inlet
    shroud_curve_raw = np.vstack([c[:,0]*np.cos(eps_ds) - c[:,1]*np.sin(eps_ds), 
                                  c[:,0]*np.sin(eps_ds) + c[:,1]*np.cos(eps_ds)]).T  # [m,m] - rotate clothoid to desired outlet angle
    p.shroud_curve = (shroud_curve_raw + np.ones_like(shroud_curve_raw)*[-p.r_outlet_impeller, -shroud_curve_raw[0,1]])*[-1,1]  # [m,m] - translate the origin to the inlet central point, and flip from quadrant III to IV
    


    ds = np.linalg.norm(np.diff(p.shroud_curve, axis=0), axis=1)    # m - differential arc length along shroud curve
    s  = np.concatenate(([0.0], np.cumsum(ds)))                     # m - arc length along shroud curve
    s_shroud = s / s[-1]                                            # unitless - normalized arc length along shroud curve
    s_shroud = s_shroud.T                                           # make into column vector

    A_inlet = np.pi*(p.r_inlet_impeller**2 - p.r_hub_impeller**2)      # m2 - annular flow area at impeller inlet
    A_outlet = 2*np.pi*p.r_outlet_impeller*w_outlet                       # m2 - annular flow area at impeller outlet
    weight = 2*np.pi*p.r_outlet_impeller/np.tan(eps_ds) - np.pi*w_outlet  # unitless - weighting factor for area distribution along meanline
    p.A_meanline = (A_outlet-A_inlet) * (2*(weight-1)*s_shroud**3 - 3*(weight-1)*s_shroud**2 + weight*s_shroud) + A_inlet # m2 - area distribution along meanline

    # tangents_shroud = np.empty_like(p.shroud_curve)
    # tangents_shroud[1:-1] = p.shroud_curve[2:] - p.shroud_curve[:-2]
    # tangents_shroud[0]    = p.shroud_curve[1] - p.shroud_curve[0]
    # tangents_shroud[-1]   = p.shroud_curve[-1] - p.shroud_curve[-2]
    tangents_shroud = np.arctan2(np.gradient(p.shroud_curve[:, 1]), np.gradient(p.shroud_curve[:, 0])) # rad - approximate tangent vectors along shroud curve
    normals_shroud = np.column_stack([np.sin(tangents_shroud), -np.cos(tangents_shroud)]) # [m,m] - outward normal vectors along shroud curve

    p.crosswise_meanline = (p.shroud_curve[:,0] - np.sqrt(p.shroud_curve[:,0]**2 - p.A_meanline*normals_shroud[:,0]/np.pi)) / normals_shroud[:,0] # m - crosswise distance between shroud and hub curves
    p.hub_curve = p.shroud_curve + p.crosswise_meanline[:, None] * normals_shroud # m - (r,z) coordinates of hub curve from outlet to inlet

    p.meanline_curve = 0.5 * (p.shroud_curve + p.hub_curve); # [m,m] - halfway between shroud and impeller

    plt.plot(p.shroud_curve[:,0],   p.shroud_curve[:,1],   label="shroud")
    plt.plot(p.hub_curve[:,0],      p.hub_curve[:,1],      label="hub")
    plt.plot(p.meanline_curve[:,0], p.meanline_curve[:,1], label="meanline")

    plt.grid(True)
    plt.axis('equal')
    plt.legend()
    plt.show()






    