from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from scipy.special import fresnel
from scipy.interpolate import CubicHermiteSpline
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
from scipy.optimize import brute, minimize
import matplotlib.pyplot as plt

from engine_state import engine, pump
from typing import List

def impellers(p: pump):
    imp = p.impeller[0]  # Currently only single stage pumps are supported
    ind = p.inducer
    imp.head_rise = p.head_rise  # assign head rise to impeller object

    #######################################################################################################
    ### Basic Sizing
    #######################################################################################################
    imp.specific_speed = abs(p.shaft_speed) * np.sqrt(p.vdot) / (engine.g * imp.head_rise ** 0.75) # unitless 

    imp.specific_speed_imperial = imp.specific_speed * 229.6        # Ns pseudo unitless
    imp.outlet_flow_coeff = 0.1715 * np.sqrt(imp.specific_speed)  # unitless - phi_i in pump handbook

    ft = 1.0; # unitless - could range to 1.1 to be more aggressive on head, on low head pumps
    imp.head_coeff = 0.5 * 1.21 * ft * np.exp(-0.408 * imp.specific_speed)  # unitless - psi_op; eqn 3.26 in Gulich THIS IS TECHNICALLY ONLY FOR SINGLE STAGE PUMPS
    imp.shutoff_coeff = 1.25; # 1.25 for volute designs and 1.31 for diffuser designs (kept as input)
    imp.shutoff_head_coeff = 0.5 * imp.shutoff_coeff * 0.625 * np.exp(-4.15 * imp.specific_speed) # unitless - psi_0 in the limit pf Q=0

    # specific speed check - should be in the centrifugal range
    if imp.specific_speed < 0.1 or imp.specific_speed > 0.6:
        print(f"{p.fluid} impeller specific speed is {round(imp.specific_speed,2)}. Recommended is 0.1 to 0.6 for centrifugal geometries.")

    imp.r_outlet = 1 / abs(p.shaft_speed) * np.sqrt (engine.g*imp.head_rise / imp.head_coeff); # m - impeller exit radius, Gulich eqn 7.1.3

    if imp.r_outlet < ind.r_tip:
        raise ValueError(f"Impeller outlet radius {round(imp.r_outlet,3)} m is smaller than inducer tip radius {round(ind.r_tip,3)} m")

    r_eye_margin = 1.15 + 0.2116*(imp.specific_speed-0.2836); # unitless - margin correction factor, Gulich 7.1.4; interpolated for a normal impeller (not a suction impeller)
    swirl_number = 0.5; # NEEDS TO BE ITERATIVE 
    
    imp.r_inlet = ind.r_tip

    imp.v_inlet = imp.flow_coeff * p.shaft_speed * imp.r_inlet # m/s
    imp.r_hub = np.sqrt(imp.r_inlet**2 - p.vdot / imp.v_inlet / np.pi) # m - impeller hub radius, from flow continuity and inducer outlet velocity; this is just a starting point for the hub radius, which will be iteratively updated to achieve the desired blade blockage

    # imp.r_inlet = imp.r_outlet * r_eye_margin * np.sqrt((imp.r_hub/imp.r_outlet)**2 + 
    #                                                                 0.00148*2*imp.head_coeff*(imp.specific_speed*52.9186)**1.33/swirl_number**0.67) 
    #     # m - impeller inlet radius, Gulich eqn 7.1.4
    # imp.eye_flow_coeff = p.vdot / (np.pi*(imp.r_inlet**2-imp.r_hub**2)*imp.r_inlet * p.shaft_speed) # unitless

    # The exit width b2
    imp.w_outlet = 2*imp.r_outlet * (0.017 + 0.1386*imp.specific_speed - 0.022387*imp.specific_speed**2 + 0.0013767*imp.specific_speed**3) 
        # m - outlet width; empirical; Gulich eqn 7.1

    # this needs boundary layer and leakage effect - pg 2.60?
    # could parameterize a width at each station with a blockage vector

    #######################################################################################################
    ### Shroud curve
    #######################################################################################################

    X = lambda t: fresnel(t/np.sqrt(np.pi))[1] # clothoid / fresnel integral
    Y = lambda t: fresnel(t/np.sqrt(np.pi))[0] # clothoid / fresnel integral
    T = np.sqrt(np.pi + 2*imp.shroud_out_angle) # upper limit of clothoid to achieve desired outlet angle

    clothoid = lambda t: (imp.r_inlet-imp.r_outlet)/(np.sin(imp.shroud_out_angle)*Y(T) - np.cos(imp.shroud_out_angle)*X(T)) * np.vstack([X(t), Y(t)]).T  # [m,m] - parametric clothoid function scaled to desired inlet and outlet radii
    c = clothoid(np.linspace(T,0,p.n_pts_meridional))    # [m,m] - (r,z) coordinates of shroud curve from outlet to inlet
    imp.shroud_curve_raw = np.vstack([c[:,0]*np.cos(imp.shroud_out_angle) + c[:,1]*np.sin(imp.shroud_out_angle), 
                                  -c[:,0]*np.sin(imp.shroud_out_angle) + c[:,1]*np.cos(imp.shroud_out_angle)]).T  # [m,m] - rotate clothoid to desired outlet angle
    imp.shroud_curve = (imp.shroud_curve_raw + np.ones_like(imp.shroud_curve_raw)*[-imp.r_outlet, -imp.shroud_curve_raw[0,1]])*[-1,1]  # [m,m] - translate the origin to the inlet central point, and flip from quadrant III to IV
    
    #######################################################################################################
    ### Hub curve
    #######################################################################################################

    imp.hub_curve = np.nan * np.ones_like(imp.shroud_curve)  # [m,m] - initialize hub curve array with NaNs
    imp.hub_curve[0,:] = [imp.r_hub, 0.0]  # m - starting point of hub curve at inlet
    imp.hub_curve[-1,:] = [imp.r_outlet, imp.shroud_curve[-1,1] - imp.w_outlet]  # m - end point of hub curve at outlet

    tangents_meanline = np.gradient(imp.shroud_curve, axis=0) # rad
    normals_meanline = np.column_stack([tangents_meanline[:, 1], -tangents_meanline[:, 0]]) / np.linalg.norm(tangents_meanline, axis=1)[:, np.newaxis] # rad

    hub_in_tangent = np.array([np.cos(imp.hub_in_angle), np.sin(imp.hub_in_angle)]) # rad
    hub_out_tangent = np.array([np.cos(imp.hub_out_angle), np.sin(imp.hub_out_angle)]) # rad

    def line_from_point_tangent(p, t): # construct a line as a vector
        return np.array([-t[1], t[0], t[1]*p[0] - t[0]*p[1]])

    L1 = line_from_point_tangent(imp.shroud_curve[0,:], normals_meanline[1,:])
    L2 = line_from_point_tangent(imp.hub_curve[0,:], hub_in_tangent)
    P_homogeneous = np.cross(L1, L2)
    imp.hub_curve[1,:] = P_homogeneous[:2] / P_homogeneous[2] # find intersection point

    L1 = line_from_point_tangent(imp.shroud_curve[-1,:], normals_meanline[-1,:])
    L2 = line_from_point_tangent(imp.hub_curve[-1,:], hub_out_tangent)
    P_homogeneous = np.cross(L1, L2)
    imp.hub_curve[-2,:] = P_homogeneous[:2] / P_homogeneous[2] # find intersection point

    gap = np.linalg.norm(imp.shroud_curve - imp.hub_curve, axis=1); # m - crosswise distance between shroud and hub along meanline    
    area = np.pi * (imp.shroud_curve[:,0] + imp.hub_curve[:,0]) * gap * imp.blockage # m2 - annular flow area along meanline
    # Interpolate areas between inlet and outlet to be a smooth curve in area space; this minimizes abrupt meridional velocity changes:
    area[np.isnan(area)] = PchipInterpolator(np.where(~np.isnan(area))[0], area[~np.isnan(area)])(np.where(np.isnan(area))[0]) # m2
    imp.gap_meanline = (imp.shroud_curve[:,0] / normals_meanline[:,0]) * (np.sqrt(1 + (area*normals_meanline[:,0]) / (np.pi*imp.shroud_curve[:,0]**2)) - 1) # m - crosswise gap along meanline
    imp.area_meanline = area # m2 - crosswise area along meanline

    imp.hub_curve = imp.shroud_curve + imp.gap_meanline[:, np.newaxis] * normals_meanline # [m,m] - hub curve generated by offsetting shroud curve by crosswise gap along meanline normals
    imp.meanline_curve = 0.5 * (imp.shroud_curve + imp.hub_curve); # [m,m] - halfway between shroud and impeller hub curves

    ds = np.linalg.norm(np.diff(imp.shroud_curve, axis=0), axis=1)    # m - differential arc length along shroud curve
    imp.arc_meanline = np.concatenate(([0.0], np.cumsum(ds))).T       # m - arc length along shroud curve
