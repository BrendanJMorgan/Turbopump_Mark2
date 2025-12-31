from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from scipy.special import fresnel
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize
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
    p.w_outlet = 2*p.r_outlet_impeller * (0.017 + 0.1386*p.specific_speed - 0.022387*p.specific_speed**2 + 0.0013767*p.specific_speed**3) 
        # m - outlet width; empirical; Gulich eqn 7.1

    # this needs boundary layer and leakage effect - pg 2.60?
    # could parameterize a width at each station with a blockage vector, linear ramp of area?

    #######################################################################################################
    ### Shroud curve
    #######################################################################################################

    n_pts_meridional = 1000

    shroud_out_angle = 10*np.pi/180 # rad
    X = lambda t: fresnel(t/np.sqrt(np.pi))[1] # clothoid / fresnel integral
    Y = lambda t: fresnel(t/np.sqrt(np.pi))[0] # clothoid / fresnel integral
    T = np.sqrt(np.pi - 2*shroud_out_angle) # upper limit of clothoid to achieve desired outlet angle

    clothoid = lambda t: (p.r_inlet_impeller-p.r_outlet_impeller)/(np.sin(shroud_out_angle)*Y(T) - np.cos(shroud_out_angle)*X(T)) * np.vstack([X(t), Y(t)]).T  # [m,m] - parametric clothoid function scaled to desired inlet and outlet radii
    c = clothoid(np.linspace(T,0,n_pts_meridional))    # [m,m] - (r,z) coordinates of shroud curve from outlet to inlet
    shroud_curve_raw = np.vstack([c[:,0]*np.cos(shroud_out_angle) - c[:,1]*np.sin(shroud_out_angle), 
                                  c[:,0]*np.sin(shroud_out_angle) + c[:,1]*np.cos(shroud_out_angle)]).T  # [m,m] - rotate clothoid to desired outlet angle
    p.shroud_curve = (shroud_curve_raw + np.ones_like(shroud_curve_raw)*[-p.r_outlet_impeller, -shroud_curve_raw[0,1]])*[-1,1]  # [m,m] - translate the origin to the inlet central point, and flip from quadrant III to IV
    
    #######################################################################################################
    ### Hub curve
    #######################################################################################################

    hub_in_angle = -90*np.pi/180 # rad
    hub_out_angle = 0*np.pi/180 # rad

    p.hub_curve = np.zeros_like(p.shroud_curve)  # [m,m]
    p.hub_curve[0,:] = [p.r_hub_impeller, 0.0]  # m - starting point of hub curve at inlet
    p.hub_curve[-1,:] = [p.r_outlet_impeller, p.shroud_curve[-1,1] - p.w_outlet]  # m - end point of hub curve at outlet

    points = np.array([[p.r_hub_impeller, 0.0],[p.r_outlet_impeller, p.shroud_curve[-1,1] - p.w_outlet]])
    angles = np.array([-np.pi/2, 0.0])

    scales_opt = minimize(area_linearity_error, 
                    x0=[0.1,0.1],  # Initial guess
                    args=(points, angles, p.shroud_curve, n_pts_meridional),
                    method='L-BFGS-B').x

    p.hub_curve = generate_hub_curve(points, angles, scales_opt, n_pts_meridional)
    p.meanline_curve = 0.5 * (p.shroud_curve + p.hub_curve); # [m,m] - halfway between shroud and impeller hub curves
    p.crosswise_meanline = np.linalg.norm(p.shroud_curve - p.hub_curve, axis=1); # m - crosswise distance between shroud and hub along meanline
    p.A_meanline = np.pi * (p.shroud_curve[:,0] + p.hub_curve[:,0]) * p.crosswise_meanline # m2 - annular flow area along meanline

    ds = np.linalg.norm(np.diff(p.shroud_curve, axis=0), axis=1)    # m - differential arc length along shroud curve
    p.arc_meanline = np.concatenate(([0.0], np.cumsum(ds))).T       # m - arc length along shroud curve

def area_linearity_error(scales, points, angles, shroud_curve, n_pts):
    """RMS error from linear area distribution."""
    hub_curve = generate_hub_curve(points, angles, scales, n_pts)
    A = np.pi * (shroud_curve[:,0] + hub_curve[:,0]) * np.sqrt( (shroud_curve[:,0] - hub_curve[:,0])**2 + (shroud_curve[:,1] - hub_curve[:,1])**2 )
    A_linear = np.linspace(A[0], A[-1], len(A))
    error = np.trapezoid((A - A_linear)**2)
    return error

def generate_hub_curve(points, angles, scales, n_pts):
    """Generate hub curve using cubic Hermite spline interpolation."""
    tangents = np.array([[scales[0]*np.cos(angles[0]), scales[0]*np.sin(angles[0])],
                         [scales[1]*np.cos(angles[1]), scales[1]*np.sin(angles[1])]])
    t_eval = np.linspace(0, 1, n_pts)
    return np.column_stack([CubicHermiteSpline([0.0, 1.0], points[:,i], tangents[:,i])(t_eval) for i in range(2)])