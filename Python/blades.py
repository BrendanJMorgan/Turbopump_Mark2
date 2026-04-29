import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_trapezoid
import math
from dataclasses import dataclass
from engine_state import engine, pump

def blades(p: pump):
    imp = p.impeller[0]  # Currently only single stage pumps are supported

    # Meanline Definition
    inlet_idx = int(0.1*imp.meanline_curve.shape[0]) # index where blades start
    imp.meanline_curve_bladed = imp.meanline_curve[inlet_idx:, :] # [m,m] - the portion of the meanline with a blade
    
    #######################################################################################################
    ### Velocity triangles
    #######################################################################################################
    imp.u_meanline = p.shaft_speed * imp.meanline_curve[:, 0] # m/s - tangential velocity along the meanline
    imp.u_bladed = p.shaft_speed * imp.meanline_curve_bladed[:, 0] # m/s - blade tangential velocity
    imp.v_merid = p.vdot / imp.area_meanline # m/s - meridional velocity
    imp.v_merid_bladed = imp.v_merid[inlet_idx:] # m/s - meridional velocity with blockage correction

    imp.v_tangential_inlet = imp.u_meanline[0] # m/s - assuming ideal swirl off inducer UPDATE WITH ACTUAL VALUE FROM INDUCER

    imp.v_tangential_outlet = engine.g * p.head_rise / (imp.hydraulic_efficiency * imp.u_bladed[-1]) + imp.v_tangential_inlet # m/s
    imp.v_tangential = np.linspace(imp.v_tangential_inlet, imp.v_tangential_outlet, imp.meanline_curve_bladed.shape[0]) # m/s - assuming linear ramp of fluid tangential velocity

    #######################################################################################################
    ### Camber Line Integration - Blade Curve Construction
    #######################################################################################################
    blade_angle = np.arctan2(imp.v_merid_bladed, (imp.u_bladed * (1.0 - imp.slip_factor) - imp.v_tangential)) # rad - angle between blade tip azimuth and local tangential azimuth (beta_2 in pump handbook)
    dm = np.sqrt(np.diff(imp.meanline_curve_bladed[:,0])**2 + np.diff(imp.meanline_curve_bladed[:,1])**2) # m - differential arc length along meanline
    imp.meanline_arc_length = np.concatenate([[0.0], np.cumsum(dm)]) # m - meanline cumulative arc length

    theta_blade = cumulative_trapezoid(1.0/np.tan(blade_angle)/imp.meanline_curve_bladed[:,0], imp.meanline_curve_bladed[:,0], initial=0)  # rad - amgular step to advance blade curve

    imp.blade_curve = np.column_stack([imp.meanline_curve_bladed[:,0]*np.cos(theta_blade),imp.meanline_curve_bladed[:,0]*-p.clocking*np.sin(theta_blade),imp.meanline_curve_bladed[:,1],])

    #######################################################################################################
    ### Number of Blades + Blockage
    #######################################################################################################
    imp.blade_arc_length = float(np.linalg.norm(np.diff(imp.blade_curve, axis=0), axis=1).sum())

    solidity_ideal = float(CubicSpline([0.0, 0.4, 3.0], [1.8, 1.8, 1.0])(imp.specific_speed))
        # unitless - solidity is the optimal ratio of blade chord to blade spacing. Pump handbook page 2.36 (sigma)
    test = solidity_ideal * (2.0 * np.pi * imp.r_outlet) / imp.blade_arc_length
    imp.blade_count = int(np.round(solidity_ideal * (2.0 * np.pi * imp.r_outlet) / imp.blade_arc_length)) # number of blades

    imp.solidity = imp.blade_count * imp.blade_arc_length / (2.0 * np.pi * imp.r_outlet) # unitless

    imp.blade_thickness = 0.04 * imp.r_outlet # m - empirical assumption: 4% of outlet radius as blade thickness
    imp.boundary_layer_thickness = 0.002 * imp.blade_arc_length # m

    imp.meridional_length = np.sum(np.sqrt(np.diff(imp.blade_curve[:, 0])**2 + np.diff(imp.blade_curve[:, 2])**2)) # m - arc length of blade projected onto the meridional plane

    imp.blockage[inlet_idx:] = 1 - imp.blade_count * (2*imp.boundary_layer_thickness + imp.blade_thickness) / (2*np.pi*imp.meanline_curve_bladed[:,0]) / np.sin(blade_angle) # unitless

    #######################################################################################################
    ### Velocity Slip
    #######################################################################################################
    
    # Pfleiderer’s Method - pump handbook pg 2.35 - Extending this to the full arc length of the meanline, instead of just the tip
    a = 0.75 # unitless - 0.65 to 0.85 for volute 
    r_ref = imp.meanline_curve_bladed[0,0] # m - scalar reference radius
    psi_prime = a * (1.0 + np.sin(blade_angle)) # unitless
    denom = imp.blade_count * 0.5 * imp.meanline_curve_bladed[:,0]**2 - r_ref**2 
    imp.slip_factor = psi_prime * (imp.meanline_curve_bladed[:,0]**2) / denom # unitless

    # Cavitation Check
    # eqns 25 and 26 in pump handbook
    # add cavitation criteria at blade inlet/outlet

    #######################################################################################################
    ### Hydraulic Efficiency
    #######################################################################################################

    v_merid_avg = (imp.v_merid_bladed[0] + imp.v_merid_bladed[-1]) / 2 # m/s - average meridional velocity along the blade
    passage_width_bladed = 2*np.pi*imp.meanline_curve_bladed[:,0]*imp.blockage[inlet_idx:] / imp.blade_count # m - circumferential distance between blades at impeller outlet
    hydraulic_diameter = 2*(imp.area_meanline[inlet_idx]+imp.area_meanline[-1]) / (passage_width_bladed[inlet_idx]+passage_width_bladed[-1]+imp.gap_meanline[inlet_idx]+imp.gap_meanline[-1]) # m - hydraulic diameter of flow passage - Gulich 3.8.3
    reynolds = v_merid_avg*imp.blade_arc_length*p.density / p.viscosity # unitless - Gulich 3.8.3
    friction_coeff = 0.136/(-np.log10(0.2*imp.surface_roughness/imp.blade_arc_length+12.5/reynolds))**2.15 # unitless - friction coefficient - Gulich 3.8.2
    dissipation_coeff = (friction_coeff+0.0015)*(1.1+4*passage_width_bladed[-1]/(2*np.pi*imp.r_outlet)) # unitless - dissipation coefficient - Gulich 3.8.5
    friction_loss = 4*dissipation_coeff*imp.blade_arc_length/hydraulic_diameter*(v_merid_avg/imp.u_bladed[-1])**2 # unitless - friction loss coefficient - Gulich 3.8.6

    w_inlet = np.sqrt(imp.v_merid[0]**2 + (imp.u_meanline[0] - imp.v_tangential_inlet)**2) # m/s - inlet relative velocity - Gulich 3.8.7
    w_throat = imp.v_merid_bladed[0] # m/s
    if w_throat/w_inlet > 0.65:
        shock_loss = 0.3*((w_inlet-w_throat)/imp.u_bladed[-1])**2 # unitless - shock loss coefficient - Gulich 3.8.8
    else:
        raise ValueError("Excessive shock loss")
    
    imp.hydraulic_efficiency = imp.head_coeff / (imp.head_coeff + friction_loss + shock_loss) # unitless - Gulich 3.8.1

    

