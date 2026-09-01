import numpy as np
from scipy.integrate import cumulative_trapezoid

from engine_state import engine, pump

def inducer(p) -> None:
    ind = p.inducer
    imp = p.impeller[0]

    ind.n_blade = 1000
    ind.n_bladeline = 20
    n_meanline = ind.n_bladeline // 2

    #######################################################################################################
    ### High Level Parameters
    #######################################################################################################
    ind.NPSH_available = (p.p_in - p.pvap_inlet) / (p.density * engine.g)  # m - Net positive suction head. Assuming zero plumbing water column and inlet velocity
    ind.flow_required = p.vdot * ind.flow_margin_factor  # m^3/s - flow rate with margin

    if not p.shaft_speed: # find shaft speed assuming the pump is suction limited
        ind.tip_speed = np.sqrt((1-2*ind.flow_coeff**2)/(3*ind.flow_coeff**2)*2*engine.g*ind.NPSH_available) # m/s - inducer tip speed - Brumfield criterion
        flow_area = ind.flow_required / (ind.flow_coeff * ind.tip_speed)  # m2 - flow area at inducer inlet
        p.shaft_speed = ind.tip_speed / np.sqrt(flow_area / (np.pi * (1 - ind.hub_tip_ratio**2)))  # rad/s - shaft speed
        return 0

    ind.r_tip = np.zeros(ind.n_blade)
    ind.r_tip[0] = (ind.flow_required / (np.pi * p.shaft_speed * ind.flow_coeff * (1 - ind.hub_tip_ratio**2)))**(1/3)
    ind.r_tip[-1] = ind.r_tip[0] # m - straight
    ind.r_tip = np.linspace(ind.r_tip[0], ind.r_tip[-1], ind.n_blade) # [m] - tip contour
    ind.r_hub = np.linspace(ind.r_tip[0]*ind.hub_tip_ratio, ind.r_tip[0]*ind.hub_tip_ratio_out, ind.n_blade) # [m] - hub contour

    ind.tip_speed = p.shaft_speed * ind.r_tip[0]   # m/s - actual tip speed at the clamped shaft speed

    ind.v_axial = np.zeros(ind.n_blade) # m/s - axial component, assumed constant
    ind.v_axial = ind.flow_required / (np.pi * (ind.r_tip**2-ind.r_hub**2))   # m/s - axial component, assumed constant THIS NEEDS A BLOCKAGE CORRECTION

    
    ind.head = (imp.NPSH_required - ind.NPSH_available)  # m
    ind.head_coeff = ind.head * engine.g / (ind.tip_speed**2)  # unitless
    suction_specific_speed_required = p.shaft_speed * np.sqrt(ind.flow_required) / (engine.g * ind.NPSH_available)**0.75  # unitless - suction specific speed

    #######################################################################################################
    ### Blade Geometry
    #######################################################################################################

    ind.blade_pitch = 2*np.pi*ind.r_tip[0]/ind.blade_number  # m - blade pitch - circumferential distance between blades at inducer tip

    # Pre-allocate  ---------------------------------------------
    for attr in ['r_bladeline', 'theta_bladeline', 'x_bladeline', 'y_bladeline', 'z_bladeline', 'blade_angle', 'flow_angle']:
        setattr(ind, attr, np.zeros((ind.n_bladeline+1, ind.n_blade))) 
    ind.blade_length =     np.zeros(ind.n_bladeline+1)
    ind.blade_length[-1] = ind.blade_pitch * ind.length_pitch_ratio  # tip

    # Flow streamlines  ---------------------------------------------
    incidence = np.deg2rad(3.0)  # 2-4 deg per Gulich
    fraction = np.arange(ind.n_bladeline+1)[:, None] / ind.n_bladeline
    ind.r_bladeline = fraction * (ind.r_tip - ind.r_hub) + ind.r_hub
    ind.flow_angle[:, 0] = np.arctan(ind.r_tip[0] / ind.r_bladeline[:, 0] * ind.flow_coeff) # rad - Gulich 7.6.8 - free vortex condition
    ind.blade_angle[:, 0] = ind.flow_angle[:, 0] + incidence

    # Outlet Conditions  ---------------------------------------------
    ind.head_efficiency = 1-0.11*ind.length_pitch_ratio # unitless - Gulich 7.6.13 - long blades = thick boundary layer at TE → more deviation, more blockage at outlet

    ind.flow_coeff_realized = ind.v_axial / (p.shaft_speed * 0.5 * (ind.r_tip[0] + ind.r_hub[0])) # unitless - average flow coefficient (not tip), for performance estimation

    # Outlet on meanline (Eq. 7.6.15 inverted)
    A_ratio = (ind.r_tip[0]**2 - ind.r_hub[0]**2) / (ind.r_tip[-1]**2 - ind.r_hub[-1]**2)
    flow_angle_meanline_out = np.arcsin(A_ratio / np.sqrt(1 + (1 - ind.head_coeff/ind.head_efficiency)/ind.flow_coeff_realized**2))

    # Meanline velocity triangle (at trailing edge)
    ind.r_meanline = ind.r_bladeline[n_meanline]
    ind.u_meanline = p.shaft_speed * ind.r_meanline
    ind.v_tang_meanline = ind.u_meanline - ind.v_axial / np.tan(flow_angle_meanline_out)

    # Eq. 7.6.14 inverted on meanline → blade angle on meanline
    blade_pitch_meanline = 2 * np.pi * ind.r_meanline[n_meanline] / ind.blade_number
    k_meanline = (blade_pitch_meanline / ind.blade_length[-1]) ** (1/3)
    ind.blade_angle[n_meanline, -1] = (
        3*ind.flow_angle[n_meanline, -1] + 6*k_meanline*np.pi/180 - k_meanline*ind.blade_angle[n_meanline, 0]) / (3 - k_meanline)
    deviation_meanline = ind.blade_angle[n_meanline, -1] - ind.flow_angle[n_meanline, -1]

    # Free-vortex propagation to all streamlines (including tip)
    ind.v_tang_local = ind.v_tang_meanline * ind.r_bladeline / ind.r_bladeline[n_meanline]
    ind.v_tang_capped = np.minimum(ind.v_tang_local, 0.7*p.shaft_speed*ind.r_bladeline) # m/s - cap tangential velocity to keep blade geometry sane

    ind.flow_angle = np.arctan(ind.v_axial / (p.shaft_speed*ind.r_bladeline - ind.v_tang_capped))
    ind.blade_angle[:, -1] = ind.flow_angle[:, -1] + deviation_meanline

    clipped = ind.blade_angle[:, -1] > np.deg2rad(40)

    #######################################################################################################
    ### Blade Construction / Interpolation
    #######################################################################################################

    weight = 3
    shape = np.linspace(0, 1, ind.n_blade)**weight # 0 at LE, 1 at TE, slow ramp early, fast ramp late

    a0 = ind.blade_angle[:, 0, None]
    a1 = ind.blade_angle[:, -1, None]
    ind.blade_angle = a0 + (a1 - a0) * shape # rad/s - user-defined distribution of blade angles

    ds = ind.blade_length[-1]/ind.n_blade # m - arc length differential
    ind.axial_height = np.trapezoid(np.sqrt(np.sin(ind.blade_angle[-1])**2 - np.gradient(ind.r_bladeline[-1], ds)**2), dx=ds)
    dz = ind.axial_height/(ind.n_blade-1) # m - axial position differential

    r_leading_edge = ind.r_bladeline[:, 0] # m - leading edge radius as a function of blade number, for calculating sweepback offset
    theta_leading_edge_offset = ind.sweepback * (r_leading_edge-r_leading_edge[0]) / (r_leading_edge[-1] - r_leading_edge[0]) # rad - sweepback offset angle at leading edge, scaled

    for n in range(ind.n_bladeline+1):
        dr_ds = np.gradient(ind.r_bladeline[n], ds) # m/m - radial change over arc length
        dz_ds = np.sin(ind.blade_angle[n]) * np.sqrt(np.maximum(1 - dr_ds**2, 0)) # m/m - axial change over arc length

        z_of_s = cumulative_trapezoid(dz_ds, dx=ds, initial=0) # [m] - axial z as a function of arc length s
        z_of_s_scaled = z_of_s * (ind.axial_height / z_of_s[-1]) # [m] - stretch bladeline to terminate at the tip exit plane

        # Resample beta and r onto uniform z
        ind.z_bladeline[n] = np.linspace(0, -ind.axial_height, ind.n_blade)   # [m] - negative projects it downwards from the inlet plane
        ind.r_bladeline[n] = np.interp(ind.z_bladeline[n], z_of_s_scaled, ind.r_bladeline[n]) # [m] - radial positions scaled and uniformly spaced in z

        ind.blade_angle[n] = np.interp(ind.z_bladeline[n], z_of_s_scaled, ind.blade_angle[n]) # [rad] - radial positions scaled and uniformly spaced in z
        dtheta_dz = np.sqrt(1 + np.gradient(ind.r_bladeline[n], dz)**2) / (ind.r_bladeline[n] * np.tan(ind.blade_angle[n]))
        ind.theta_bladeline[n] = -theta_leading_edge_offset[n] - cumulative_trapezoid(dtheta_dz, dx=dz, initial=0) # [m]

    ind.x_bladeline = ind.r_bladeline*np.cos(ind.theta_bladeline) # [m]
    ind.y_bladeline = ind.r_bladeline*np.sin(ind.theta_bladeline) # [m]

    #######################################################################################################
    ### Performance Prediction (approximation before CFD)
    #######################################################################################################

    ind.blockage = 0.98 # DUMMY

    # Gulich empirical correlation - equation 7.23a
    x_sigma = (ind.flow_coeff**0.43 * (np.tan(ind.blade_angle[-1, 0])/ind.blockage)**0.75 * min(ind.length_pitch_ratio, 2.0)**-0.2 *
            ind.hub_tip_ratio**0.3 * (ind.r_hub[-1]/ind.r_tip[-1])**0.1 * (1 + ind.sweepback*180/np.pi/100)**(-0.47) * (ind.blade_number/3)**0.24)

    sigma_3 = 1.26 * x_sigma
    NPSH_3_predicted = sigma_3 * ind.tip_speed**2 / (2 * engine.g)
    raw_Nss = p.shaft_speed * np.sqrt(p.vdot) / (engine.g * NPSH_3_predicted)**0.75
    k4 = (0.5+0.65)/2 # unitless - empirical coefficient
    ind.suction_specific_speed_available = raw_Nss * (1 - k4*np.sqrt(ind.tip_clearance/(ind.r_tip[0]-ind.r_hub[0])))

    if ind.suction_specific_speed_available < suction_specific_speed_required:
        raise ValueError(f"{ind.fluid} inducer ωss of {ind.suction_specific_speed_available:0.2f} is too low for the required value of {suction_specific_speed_required:0.2f}\nshaft speed of {30/np.pi*p.shaft_speed:0.2f} RPM")

