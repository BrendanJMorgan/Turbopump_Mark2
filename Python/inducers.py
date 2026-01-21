import numpy as np
from scipy.optimize import minimize

from engine_state import engine, pump

def inducers(p) -> None:
    ind = p.inducer
    imp = p.impeller[0]

    ind.NPSH_available = (p.p_in - p.pvap_inlet) / (p.density * engine.g)  # m - Net positive suction head. Assuming zero plumbing water column and inlet velocity
    ind.NPSH_required = ind.NPSH_available / ind.suction_margin  # m - NPSH required with margin

    # Brumfield Criterion  ---------------------------------------------
    ind.tip_speed = np.sqrt((1-2*ind.flow_coeff**2)/(3*ind.flow_coeff**2)*2*engine.g*ind.NPSH_required) # m/s - inducer tip speed
    flow_corrected = p.vdot * ind.flow_margin  # m^3/s - flow rate with margin
    flow_area = p.vdot / (ind.flow_coeff * ind.tip_speed)  # m^2 - flow area at inducer inlet

    if not p.shaft_speed:
        ind.r_tip = np.sqrt(flow_area / (np.pi * (1 - ind.hub_tip_ratio**2)))  # m - inducer tip radius
        ind.r_hub_in = ind.r_tip * ind.hub_tip_ratio  # m - inducer hub radius
        p.shaft_speed = ind.tip_speed / ind.r_tip  # rad/s - shaft speed
        return 0

    ind.r_tip = ind.tip_speed / p.shaft_speed  # m - inducer tip radius
    ind.r_hub_in = ind.r_tip * ind.hub_tip_ratio  # m - inducer hub radius

    ind.head = (imp.NPSH_required - ind.NPSH_available) * imp.suction_margin  # m
    ind.head_coeff = ind.head * engine.g / (ind.tip_speed**2)  # unitless
    ind.suction_specific_speed = p.shaft_speed * np.sqrt(flow_corrected) / (engine.g * ind.NPSH_required)**0.75  # unitless - suction specific speed

    #######################################################################################################
    ### Blade Construction
    #######################################################################################################
    ind.r_hub_out = ind.r_hub_in  # m - inducer hub radius at outlet, set to impeller inlet radius for smooth flow transition
    ind.blade_angle_in = np.arctan(ind.flow_coeff)  # rad - angle between blade tip azimuth and local tangential azimuth (beta_2 in pump handbook)
    ind.blade_pitch = 2 * np.pi * ind.r_tip / ind.blade_number  # m - blade pitch - circumferential distance between blades at inducer tip

    num_streamlines = 14
    ind.r0_streamline = np.linspace(ind.r_hub_out, ind.r_tip, num_streamlines)  # m - radial locations along blade center streamline
    ind.beta_streamline0 = np.arctan(ind.flow_coeff*ind.r_hub_out/ind.r0_streamline)  # rad - streamline angle at blade inlet

    [ind.blade_angle, ind.blade_length] = minimize(psi_error, 
                x0=[0.05,ind.blade_pitch],  # Initial guess
                args=(ind.blade_pitch, ind.flow_coeff, ind.r_tip, ind.r_hub_in, ind.r_hub_out, ind.beta_streamline0[int(num_streamlines/2)], ind.head_coeff),
                method='L-BFGS-B').x

    p.inducer = ind

def psi_error(params, blade_pitch, flow_coeff, r_tip, r_hub_in, r_hub_out, beta0, head_coeff):
    return np.abs(find_psi(params, blade_pitch, flow_coeff, r_tip, r_hub_in, r_hub_out, beta0) - head_coeff)

def find_psi(params, blade_pitch, flow_coeff, r_tip, r_hub_in, r_hub_out, beta0):
    blade_angle, length = params
    head_efficiency = 1-0.11*length/blade_pitch # unitless - Gulich 7.6.13
    deviation_angle = (2+(blade_angle-beta0)) * (blade_pitch/length)**(1/3) # rad - flow deviation from blade, Gulich 7.6.14
    psi = head_efficiency * (1-flow_coeff**2)*((r_tip-r_hub_in)/(r_tip-r_hub_out)**2 / np.sin(blade_angle-deviation_angle)**2 - 1) # unitless - head coefficient, Gulich 7.6.15
    return psi
    
    # v_in = np.array([0.0, 0.0])
    # v_in[0] = p.shaft_speed * p.r_tip  # m/s - tangential velocity at inducer pitchline tip
    # v_in[1] = p.vdot_corrected / (np.pi * p.r_tip**2)  # m/s - axial velocity at inducer pitchline
    # p.v_inlet = v_in

    # p.flow_blade_ratio = 0.35  # unitless
    # p.flow_incidence_angle = np.arctan(v_in[1] / v_in[0])  # rad - alpha
    # p.blade_angle = p.flow_incidence_angle / (1 - p.flow_blade_ratio)  # rad - beta

    # p.blade_lead = 2 * np.pi * p.r_tip * np.tan(p.blade_angle)  # m - blade lead - how much the blade advances per turn
    # p.solidity = 2.5
    # p.h_min = (p.blade_lead * p.solidity / p.blade_number * np.sin(p.blade_angle))  # m - height of inducer

    # # Axial Clearance
    # p.clearance_axial = p.blade_lead / p.blade_number  # m

    # # Head
    # p.v_tip = p.shaft_speed * p.r_tip  # m/s
    # p.head_coeff = p.head * engine.g / (p.v_tip**2)  # unitless - 0.15 is the cutoff between low and high head - camber needed if >0.075

    # # Hub Construction
    # p.r_hub_base = p.r_hub_impeller  # m
    # p.r_hub_top = 2 * p.hub_tip_ratio * p.r_tip - p.r_hub_base  # m

    # # Clearance Losses - not relevant with a proper shroud
    # p.r_cavity = p.r_tip + p.clearance_radial  # m
    # p.ss_speed_cleared = p.ss_speed * (1 - 0.575 * np.sqrt((p.r_cavity - p.r_tip) / (p.r_tip - p.r_hub)))  # unitless
    # p.head_coeff_cleared = p.head_coeff * (1 - 1.0 * np.sqrt((p.r_cavity - p.r_tip) / (p.r_tip - p.r_hub)))  # unitless

    # # Other
    # p.inducer_torque = p.mdot * p.shaft_speed * p.r_tip**2  # N*m - this is a conservative maximum
