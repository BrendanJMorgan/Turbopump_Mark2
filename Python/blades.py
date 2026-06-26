import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_trapezoid
import math
from dataclasses import dataclass
from engine_state import engine, pump

def blades(p: pump):
    imp = p.impeller[0]  # Currently only single stage pumps are supported
    b = imp.blade

    # Meanline Definition
    inlet_idx = int(0.1*imp.meanline_curve.shape[0]) # index where blades start
    b.meanline_curve = imp.meanline_curve[inlet_idx:, :] # [m,m] - the portion of the meanline with a blade
    
    #######################################################################################################
    ### Velocity triangles
    #######################################################################################################
    imp.u_meanline = p.shaft_speed * imp.meanline_curve[:, 0] # m/s - tangential velocity along the meanline
    b.u = p.shaft_speed * b.meanline_curve[:, 0] # m/s - blade tangential velocity
    imp.v_merid = p.vdot / imp.area_meanline # m/s - meridional velocity
    b.v_merid = imp.v_merid[inlet_idx:] # m/s - meridional velocity with blockage correction

    # imp.v_tangential_inlet = ind.v_meanline[-1]
    imp.v_tangential_inlet = imp.u_meanline[0] # m/s - assuming ideal swirl off inducer UPDATE WITH ACTUAL VALUE FROM INDUCER

    imp.v_tangential_outlet = engine.g * p.head_rise / (imp.hydraulic_efficiency * b.u[-1]) + imp.v_tangential_inlet # m/s
    b.v_tangential = np.linspace(imp.v_tangential_inlet, imp.v_tangential_outlet, b.meanline_curve.shape[0]) # m/s - assuming linear ramp of fluid tangential velocity
    imp.c_bladed = np.sqrt(b.v_merid**2 + b.v_tangential**2) # m/s - absolute velocity along the bladed portion of the meanline

    #######################################################################################################
    ### Camber Line Integration - Blade Curve Construction
    #######################################################################################################
    b.angle = np.arctan2(b.v_merid, (b.u * (1.0 - imp.slip_factor) - b.v_tangential)) # rad - angle between blade tip azimuth and local tangential azimuth (beta_2 in pump handbook)
    dm = np.sqrt(np.diff(b.meanline_curve[:,0])**2 + np.diff(b.meanline_curve[:,1])**2) # m - differential arc length along meanline
    imp.meanline_arc_length = np.concatenate([[0.0], np.cumsum(dm)]) # m - meanline cumulative arc length

    theta_blade = cumulative_trapezoid(1.0/np.tan(b.angle)/b.meanline_curve[:,0], b.meanline_curve[:,0], initial=0)  # rad - amgular step to advance blade curve

    imp.blade_curve = np.column_stack([b.meanline_curve[:,0]*np.cos(theta_blade),b.meanline_curve[:,0]*-p.clocking*np.sin(theta_blade),b.meanline_curve[:,1],])

    #######################################################################################################
    ### Number of Blades + Blockage
    #######################################################################################################
    b.ds = np.linalg.norm(np.diff(imp.blade_curve, axis=0), axis=1) # m
    b.arc_length = np.concatenate([[0.0], np.cumsum(b.ds)]) # m - cumulative arc length along the blade curve (length N)
    b.arc_length_total = float(b.ds.sum()) # m

    solidity_ideal = float(CubicSpline([0.0, 0.4, 3.0], [1.8, 1.8, 1.0])(imp.specific_speed))
        # unitless - solidity is the optimal ratio of blade chord to blade spacing. Pump handbook page 2.36 (sigma)
    test = solidity_ideal * (2.0 * np.pi * imp.r_outlet) / b.ds
    imp.blade_count = int(np.round(solidity_ideal * (2.0 * np.pi * imp.r_outlet) / b.arc_length_total)) # number of blades

    imp.solidity = imp.blade_count * b.arc_length_total / (2.0 * np.pi * imp.r_outlet) # unitless

    imp.blade_thickness = 0.04 * imp.r_outlet # m - empirical assumption: 4% of outlet radius as blade thickness
    imp.boundary_layer_thickness = 0.002 * b.arc_length_total # m

    imp.meridional_length = np.sum(np.sqrt(np.diff(imp.blade_curve[:, 0])**2 + np.diff(imp.blade_curve[:, 2])**2)) # m - arc length of blade projected onto the meridional plane

    imp.blockage[inlet_idx:] = 1 - imp.blade_count * (2*imp.boundary_layer_thickness + imp.blade_thickness) / (2*np.pi*b.meanline_curve[:,0]) / np.sin(b.angle) # unitless

    b.area_meanline = imp.area_meanline[inlet_idx:] # m2 - crosswise area along meanline

    #######################################################################################################
    ### Velocity Slip
    #######################################################################################################
    
    # Pfleiderer’s Method - pump handbook pg 2.35 - Extending this to the full arc length of the meanline, instead of just the tip
    a = 0.75 # unitless - 0.65 to 0.85 for volute 
    r_ref = b.meanline_curve[0,0] # m - scalar reference radius
    psi_prime = a * (1.0 + np.sin(b.angle)) # unitless
    denom = imp.blade_count * 0.5 * b.meanline_curve[:,0]**2 - r_ref**2 
    imp.slip_factor = psi_prime * (b.meanline_curve[:,0]**2) / denom # unitless

    # Cavitation Check
    # eqns 25 and 26 in pump handbook
    # add cavitation criteria at blade inlet/outlet

    #######################################################################################################
    ### Hydraulic Efficiency
    #######################################################################################################

    b.v_merid_avg = np.trapz(b.v_merid, b.arc_length) / b.arc_length_total # m/s - average meridional velocity along the blade
    b.passage_width = 2*np.pi*b.meanline_curve[:,0]*imp.blockage[inlet_idx:] / imp.blade_count # m - circumferential distance between blades at impeller outlet
    b.area_meanline_avg = np.trapezoid(b.area_meanline, b.arc_length) / b.arc_length_total # m2 - average crosswise area along the bladed portion of the meanline
    b.hydraulic_diameter_avg = np.sqrt(4/np.pi*b.area_meanline_avg) # m - hydraulic diameter of flow passage - Gulich 3.8.4
    
    reynolds = b.v_merid_avg*b.arc_length_total*p.density/p.viscosity # unitless - Gulich 3.8.3
    friction_coeff = 0.136/(-np.log10(0.2*imp.surface_roughness/b.arc_length_total+12.5/reynolds))**2.15 # unitless - friction coefficient - Gulich 3.8.3
    dissipation_coeff = (friction_coeff+0.0015)*(1.1+4*b.passage_width[-1]/(2*np.pi*imp.r_outlet)) # unitless - dissipation coefficient - Gulich 3.8.5
    friction_loss = 4*dissipation_coeff*b.arc_length_total/b.hydraulic_diameter_avg*(b.v_merid_avg/b.u[-1])**2 # unitless - friction loss coefficient - Gulich 3.8.6

    b.u_meanline = p.shaft_speed * b.meanline_curve[:,0] # m/s - tangential velocity along the meanline
    b.w_inlet = np.sqrt(b.v_merid[0]**2 + (b.u_meanline[0] - b.v_tangential[0])**2) # m/s - inlet relative velocity - Gulich 3.8.7
    if  b.v_merid[0]/b.w_inlet > 0.65:
        shock_loss = 0.3*((b.w_inlet-b.v_merid[0])/b.u[-1])**2 # unitless - shock loss coefficient - Gulich 3.8.8
    else:
        shock_loss = 0.0
    
    imp.total_loss = friction_loss + shock_loss # unitless - total loss coefficient - Gulich 3.8.10

    

