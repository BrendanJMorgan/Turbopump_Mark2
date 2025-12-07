import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_trapezoid
import math
from dataclasses import dataclass
from engine_state import engine, pump

def blades(p: pump):
    # Meanline Definition
    r_inlet_blade = p.r_inlet_impeller # m - radial location of blade index - typically vertically aligned with eye radius - this assumption can be changed

    inlet_idx = int(np.argmax(p.meanline_curve[:, 0] > r_inlet_blade))

    p.meanline_curve_bladed = p.meanline_curve[inlet_idx:, :] # [m,m] - the portion of the meanline with a blade

    # Velocity Triangles
    p.u_blade = p.shaft_speed * p.meanline_curve_bladed[:, 0] # m/s - blade tangential velocity
    p.v_merid = p.vdot / (2*math.pi*p.meanline_curve_bladed[:,0]*np.linalg.norm(p.shroud_curve-p.hub_curve, axis=1)) # m/s - meridional velocity

    p.hydraulic_efficiency = 1.0 - 0.071 * (p.vdot ** 0.25) # unitless - Jekat's Empirical Formula - valid for all specific speeds CONVERT TO ANDERSON?

    p.v_tangential_inlet = 0.0 # m/s - assuming no pre-swirl from the inducer

    p.v_tangential_outlet = engine.g * p.head_rise / (p.hydraulic_efficiency * p.u_blade[-1]) + p.v_tangential_inlet # m/s
    v_tangential = np.linspace(p.v_tangential_inlet, p.v_tangential_outlet, p.meanline_curve_bladed.shape[0]) # m/s - assuming linear ramp of fluid tangential velocity
    fluid_azimuth = np.arctan(p.v_merid /(p.u_blade - v_tangential)) # rad - angle between fluid flow azimuth and local tangential azimuth (beta_f2 in pump handbook)

    # Camber Line Integration - Blade Curve Construction ---------------------------------------------
    blade_azimuth = np.arctan2(p.v_merid, (p.u_blade * (1.0 - p.slip_factor) - v_tangential)) # rad - angle between blade tip azimuth and local tangential azimuth (beta_2 in pump handbook)
    dm = np.sqrt(np.diff(p.meanline_curve_bladed[:,0])**2 + np.diff(p.meanline_curve_bladed[:,1])**2) # m - differential arc length along meanline
    p.meanline_arc_length = np.concatenate([[0.0], np.cumsum(dm)]) # m - meanline cumulative arc length

    theta_blade = np.concatenate([[0.0], cumulative_trapezoid(1.0 / np.tan(blade_azimuth), p.meanline_curve_bladed[:,0])]) / p.meanline_curve_bladed[:,0] # rad

    blade_curve = np.column_stack([
        p.meanline_curve_bladed[:,0] * np.cos(theta_blade),-p.clocking * np.sin(theta_blade),p.meanline_curve_bladed[:,1],
    ])

    # Number of Blades -------------------------------------------------------------
    db1 = np.diff(blade_curve[:, 0])
    db2 = np.diff(blade_curve[:, 1])
    db3 = np.diff(blade_curve[:, 2])
    blade_arc_length = float(np.sum(np.sqrt(db1**2 + db2**2 + db3**2)))

    solidity_ideal = float(CubicSpline([0.0, 0.4, 3.0], [1.8, 1.8, 1.0])(p.specific_speed))
        # unitless - solidity is the optimal ratio of blade chord to blade spacing. Pump handbook page 2.36 (sigma)

    blade_count = int(np.round(solidity_ideal * (2.0 * np.pi * p.r_outlet_impeller) / blade_arc_length)) # number of blades

    solidity = blade_count * blade_arc_length / (2.0 * np.pi * p.r_outlet_impeller) # unitless

    p.blade_arc_length = blade_arc_length
    p.solidity_ideal = float(solidity_ideal)
    p.blade_count = blade_count
    p.solidity = float(solidity)

    # Compute Blockage (pump handbook page 2.63) ---------------------------------------------
    p.blade_thickness = 0.04 * p.r_outlet_impeller # m - empirical assumption: 4% of outlet radius as blade thickness
    p.boundary_layer_thickness = 0.002 * blade_arc_length # m

    p.meridional_length = np.sum(np.sqrt(np.diff(blade_curve[:, 0])**2 + np.diff(blade_curve[:, 2])**2)) # m - arc length of blade projected onto the meridional plane

    p.blockage = 1.0 - solidity / p.meanline_arc_length[-1] * (p.boundary_layer_thickness + p.blade_thickness)

    # Velocity Slip - Pfleiderer’s Method - pump handbook pg 2.35
    # Extending this to the full arc length of the meanline, instead of just the tip
    a = 0.75 # unitless - 0.65 to 0.85 for volute 
    r_ref = p.meanline_curve_bladed[0,0] # m - scalar reference radius
    psi_prime = a * (1.0 + np.sin(blade_azimuth)) # unitless
    denom = blade_count * 0.5 * p.meanline_curve_bladed[:,0]**2 - r_ref**2 
    slip_factor = psi_prime * (p.meanline_curve_bladed[:,0]**2) / denom # unitless

    p.slip_factor = slip_factor  # update for next outer iteration (volute loop)

    # Cavitation Check
    # eqns 25 and 26 in pump handbook
    # add cavitation criteria at blade inlet/outlet
