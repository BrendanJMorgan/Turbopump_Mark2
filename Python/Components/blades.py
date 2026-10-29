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
    ind = p.inducer

    # Meanline Definition
    inlet_idx = int(0.1*imp.meanline_curve.shape[0]) # index where blades start SHOULD LESS ARBITRARY
    b.meanline_curve = imp.meanline_curve[inlet_idx:, :] # [m,m] - the portion of the meanline with a blade
    
    #######################################################################################################
    ### Velocity triangles - Inlet and Outlet
    #######################################################################################################
    imp.u_meanline = p.shaft_speed * imp.meanline_curve[:, 0] # m/s - tangential velocity along the meanline
    b.u = p.shaft_speed * b.meanline_curve[:, 0] # m/s - blade tangential velocity
    imp.v_merid = p.vdot / imp.area_meanline # m/s - meridional velocity
    b.v_merid = imp.v_merid[inlet_idx:] # m/s - meridional velocity with blockage correction

    imp.v_tangential_inlet = ind.v_tang_meanline[-1] * ind.r_meanline[-1] / b.meanline_curve[0,0] # m/s - r*cu conserved from inducer meanline to impeller blade start
    # imp.v_tangential_inlet = imp.u_meanline[0] # m/s - assuming ideal swirl off inducer UPDATE WITH ACTUAL VALUE FROM INDUCER

    imp.v_tangential_outlet = (engine.g * imp.head_rise / imp.hydraulic_efficiency + b.u[0]*imp.v_tangential_inlet)/b.u[-1] # m/s

    b.v_tangential = np.linspace(imp.v_tangential_inlet, imp.v_tangential_outlet, b.meanline_curve.shape[0]) # m/s - assuming linear ramp of fluid tangential velocity
    imp.c_bladed = np.sqrt(b.v_merid**2 + b.v_tangential**2) # m/s - absolute velocity along the bladed portion of the meanline

    #######################################################################################################
    ### Camber Line Integration - Blade Curve Construction
    #######################################################################################################
    b.slip_factor = imp.slip_factor[inlet_idx:] # unitless - slip factor along the bladed portion of the meanline
    b.angle = np.arctan2(b.v_merid, (b.u * (1.0 - b.slip_factor) - b.v_tangential)) # rad - angle between blade tip azimuth and local tangential azimuth (beta_2 in pump handbook)
    dm = np.sqrt(np.diff(b.meanline_curve[:,0])**2 + np.diff(b.meanline_curve[:,1])**2) # m - differential arc length along meanline
    imp.meanline_arc_length = np.concatenate([[0.0], np.cumsum(dm)]) # m - meanline cumulative arc length

    theta_blade = cumulative_trapezoid(1.0/np.tan(b.angle)/b.meanline_curve[:,0], imp.meanline_arc_length, initial=0)  # rad - dtheta = dm / (r tan(beta)), integrated over meridional length

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
    b.count = int(np.round(solidity_ideal * (2.0 * np.pi * imp.r_outlet) / b.arc_length_total)) # number of blades

    imp.solidity = b.count * b.arc_length_total / (2.0 * np.pi * imp.r_outlet) # unitless

    imp.blade_thickness = imp.blade_thickness_ratio * 2*imp.r_outlet # m
    imp.boundary_layer_thickness = 0.002 * b.arc_length_total # m

    imp.meridional_length = np.sum(np.sqrt(np.diff(imp.blade_curve[:, 0])**2 + np.diff(imp.blade_curve[:, 2])**2)) # m - arc length of blade projected onto the meridional plane

    imp.blockage[inlet_idx:] = 1 - b.count * (2*imp.boundary_layer_thickness + imp.blade_thickness) / (2*np.pi*b.meanline_curve[:,0]) / np.sin(b.angle) # unitless

    b.area_meanline = imp.area_meanline[inlet_idx:] # m2 - crosswise area along meanline

    #######################################################################################################
    ### Velocity Slip
    #######################################################################################################
    
    # Pfleiderer’s Method - pump handbook pg 2.35 - Extending this to the full arc length of the meanline, instead of just the tip
    # a = 0.75 # unitless - 0.65 to 0.85 for volute 
    # r_ref = b.meanline_curve[0,0] # m - scalar reference radius
    # psi_prime = a * (1.0 + np.sin(b.angle)) # unitless
    # denom = b.count * 0.5 * b.meanline_curve[:,0]**2 - r_ref**2 
    # imp.slip_factor = psi_prime * (b.meanline_curve[:,0]**2) / denom # unitless

    # Gulich Slip
    epsilon_lim = np.exp(-8.16*np.sin(b.angle)/b.count) # pseudo  unitless - Gulich 3.2.4
    kw = 1 - ((b.meanline_curve[0,0]/b.meanline_curve[:,0]-epsilon_lim)/(1-epsilon_lim))**3 # unitless - Gulich 3.2.5
    imp.slip_factor[-1] = 1 - 0.98 * (1 - np.sqrt(np.sin(b.angle[-1])) / b.count**0.7) # unitless - 1 - gamma; Gulich 3.2.6 gives gamma (c_2u ~ gamma*u_2), stored here as the fractional reduction used in u*(1-slip)

    # Cavitation Check
    # eqns 25 and 26 in pump handbook
    # add cavitation criteria at blade inlet/outlet

    #######################################################################################################
    ### Hydraulic Efficiency
    #######################################################################################################

    b.v_merid_avg = np.trapz(b.v_merid, b.arc_length) / b.arc_length_total # m/s - average meridional velocity along the blade
    b.passage_width = 2*np.pi*b.meanline_curve[:,0]*imp.blockage[inlet_idx:] / b.count # m - circumferential distance between blades at impeller outlet
    b.area_meanline_avg = np.trapezoid(b.area_meanline, b.arc_length) / b.arc_length_total # m2 - average crosswise area along the bladed portion of the meanline
    b.throat_width_inlet  = 2*np.pi*b.meanline_curve[0,0]*np.sin(b.angle[0])/b.count - imp.blade_thickness   # m - a1, normal to the blades
    b.throat_width_outlet = 2*np.pi*b.meanline_curve[-1,0]*np.sin(b.angle[-1])/b.count - imp.blade_thickness # m - a2
    b.height_inlet = imp.gap_meanline[inlet_idx] # m - b1
    b.throat_area_inlet = b.throat_width_inlet*b.height_inlet # m2 - A1q
    b.hydraulic_diameter_avg = 2*(b.throat_width_outlet*imp.width_outlet + b.throat_area_inlet)/(b.throat_width_inlet + b.height_inlet + b.throat_width_outlet + imp.width_outlet) # m - Gulich 3.8.4
    
    b.w = np.sqrt(b.v_merid**2 + (b.u - b.v_tangential)**2) # m/s - relative velocity along the blade
    b.w_avg = np.trapezoid(b.w, b.arc_length) / b.arc_length_total # m/s - average relative velocity along the blade

    reynolds = b.w_avg*b.arc_length_total*p.density/p.viscosity # unitless - Gulich 3.8.3
    friction_coeff = 0.136/(-np.log10(0.2*imp.surface_roughness/b.arc_length_total+12.5/reynolds))**2.15 # unitless - friction coefficient - Gulich 3.8.3
    dissipation_coeff = (friction_coeff+0.0015)*(1.1+4*imp.width_outlet/(2*imp.r_outlet)) # unitless - dissipation coefficient - Gulich 3.8.5 (b2/d2)
    friction_loss = 0.5*4*dissipation_coeff*b.arc_length_total/b.hydraulic_diameter_avg*(b.w_avg/b.u[-1])**2 # unitless - friction loss coefficient - Gulich 3.8.6

    b.u_meanline = p.shaft_speed * b.meanline_curve[:,0] # m/s - tangential velocity along the meanline
    b.w_inlet = np.sqrt(b.v_merid[0]**2 + (b.u_meanline[0] - b.v_tangential[0])**2) # m/s - inlet relative velocity - Gulich 3.8.7
    b.w_throat_inlet = imp.vdot/(b.count*b.throat_area_inlet) # m/s - w1q, single entry (f_q = 1) - Gulich 3.8.8
    if  b.w_throat_inlet/b.w_inlet > 0.65:
        shock_loss = 0.5*0.3*((b.w_inlet-b.w_throat_inlet)/b.u[-1])**2 # unitless - shock loss coefficient - Gulich 3.8.9
    else:
        shock_loss = 0.5*0.0
    
    imp.total_loss = friction_loss + shock_loss # unitless - total loss coefficient - Gulich 3.8.10

    

