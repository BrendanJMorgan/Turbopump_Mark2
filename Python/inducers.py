import numpy as np
from scipy.optimize import minimize
from scipy.integrate import cumulative_trapezoid

from engine_state import engine, pump

def inducer(p) -> None:
    ind = p.inducer
    imp = p.impeller[0]

    #######################################################################################################
    ### High Level Parameters
    #######################################################################################################
    ind.NPSH_available = (p.p_in - p.pvap_inlet) / (p.density * engine.g)  # m - Net positive suction head. Assuming zero plumbing water column and inlet velocity
    ind.NPSH_required = ind.NPSH_available / ind.suction_margin  # m - NPSH required with margin

    # Brumfield Criterion  ---------------------------------------------
    ind.tip_speed = np.sqrt((1-2*ind.flow_coeff**2)/(3*ind.flow_coeff**2)*2*engine.g*ind.NPSH_required) # m/s - inducer tip speed
    flow_corrected = p.vdot * ind.flow_margin  # m^3/s - flow rate with margin
    flow_area = p.vdot / (ind.flow_coeff * ind.tip_speed)  # m^2 - flow area at inducer inlet

    if not p.shaft_speed:
        r_tip_in = np.sqrt(flow_area / (np.pi * (1 - ind.hub_tip_ratio**2)))  # m - inducer tip radius
        r_hub_in = r_tip_in * ind.hub_tip_ratio  # m - inducer hub radius
        p.shaft_speed = ind.tip_speed / r_tip_in  # rad/s - shaft speed
        return 0

    r_tip_in = ind.tip_speed / p.shaft_speed  # m - inducer tip radius
    r_hub_in = r_tip_in * ind.hub_tip_ratio  # m - inducer hub radius
    r_hub_out = r_tip_in * ind.hub_tip_ratio_out  # m - inducer hub radius at outlet, set to impeller inlet radius for smooth flow transition

    ind.head = (imp.NPSH_required - ind.NPSH_available)  # m
    ind.head_coeff = ind.head * engine.g / (ind.tip_speed**2)  # unitless
    suction_specific_speed_required = p.shaft_speed * np.sqrt(flow_corrected) / (engine.g * ind.NPSH_required)**0.75  # unitless - suction specific speed
    k4 = (0.5+0.65)/2 # unitless - empirical coefficient
    ind.suction_specific_speed_corrected = suction_specific_speed_required * (1 - k4*np.sqrt(ind.tip_clearance/(r_tip_in-r_hub_in)))

    #######################################################################################################
    ### Blade Geometry
    #######################################################################################################
    ind.blade_pitch = 2 * np.pi * r_tip_in / ind.blade_number  # m - blade pitch - circumferential distance between blades at inducer tip

    n_blade = 1000
    num_bladeline = 21
    n_bladeline = num_bladeline-1
    n_meanline = int(n_bladeline/2)

    # Pre-allocate  ---------------------------------------------
    ind.r_bladeline =      np.zeros((n_bladeline+1, n_blade))
    ind.theta_bladeline =  np.zeros((n_bladeline+1, n_blade))
    ind.z_bladeline =      np.zeros((n_bladeline+1, n_blade))
    ind.x_bladeline =      np.zeros((n_bladeline+1, n_blade))
    ind.y_bladeline =      np.zeros((n_bladeline+1, n_blade))
    ind.blade_angle  =      np.zeros((n_bladeline+1, n_blade))
    ind.flow_angle  =       np.zeros((n_bladeline+1, n_blade))
    ind.blade_length =      np.zeros(n_bladeline+1)
    
    ind.blade_length[-1] = ind.blade_pitch * ind.length_pitch_ratio  # tip

    # Starting conditions - inlet plane and meridional plane  ---------------------------------------------
    ind.r_hub = np.linspace(r_hub_in, r_hub_out, n_blade) # [m] - hub contour
    ind.r_tip = np.linspace(r_tip_in, r_tip_in, n_blade)  # [m] - tip contour
    incidence = np.deg2rad(3.0)  # 2-4 deg per Gulich
    for n in range(n_bladeline+1):
        ind.r_bladeline[n] = n/n_bladeline * (ind.r_tip - ind.r_hub) + ind.r_hub # m
        ind.flow_angle[n,0] = np.arctan(ind.r_tip[0]/ind.r_bladeline[n,0]*ind.flow_coeff) # rad - Gulich 7.6.8 - free vortex condition
        ind.blade_angle[n,0] = ind.flow_angle[n,0] + incidence
   
    # Outlet Conditions  ---------------------------------------------
    ind.head_efficiency = 1-0.11*ind.length_pitch_ratio # unitless - Gulich 7.6.13 - long blades = thick boundary layer at TE → more deviation, more blockage at outlet

    c_axial = ind.flow_coeff * ind.tip_speed   # m/s - axial component, assumed constant   
    ind.flow_coeff_meanline = c_axial / (p.shaft_speed * 0.5 * (ind.r_tip[0] + r_hub_in)) # unitless - average flow coefficient (not tip), for performance estimation

    # Outlet on meanline (Eq. 7.6.15 inverted)
    A_ratio = (ind.r_tip[-1]**2 - r_hub_in**2) / (ind.r_tip[-1]**2 - r_hub_out**2)
    ind.flow_angle[n_meanline, -1] = np.arcsin(
        A_ratio / np.sqrt(1 + (1 - ind.head_coeff/ind.head_efficiency)/ind.flow_coeff_meanline**2)
    )

    # Meanline velocity triangle (at TE)
    r_meanline_TE = ind.r_bladeline[n_meanline, -1]
    u_meanline = p.shaft_speed * r_meanline_TE
    c_tang_meanline = u_meanline - c_axial / np.tan(ind.flow_angle[n_meanline, -1])

    # Eq. 7.6.14 inverted on meanline → blade angle on meanline
    blade_pitch_meanline = 2 * np.pi * r_meanline_TE / ind.blade_number
    k_meanline = (blade_pitch_meanline / ind.blade_length[-1]) ** (1/3)
    ind.blade_angle[n_meanline, -1] = (
        3*ind.flow_angle[n_meanline, -1] + 6*k_meanline*np.pi/180 - k_meanline*ind.blade_angle[n_meanline, 0]
    ) / (3 - k_meanline)
    deviation_meanline = ind.blade_angle[n_meanline, -1] - ind.flow_angle[n_meanline, -1]

    # Free-vortex propagation to all streamlines (including tip)
    for n in range(n_bladeline + 1):
        r_TE = ind.r_bladeline[n, -1]
        u_local = p.shaft_speed * r_TE
        c_tang_local = c_tang_meanline * r_meanline_TE / r_TE
        c_tang_capped = min(c_tang_local, 0.7 * u_local)
        
        ind.flow_angle[n, -1] = np.arctan(c_axial / (u_local - c_tang_capped))
        ind.blade_angle[n, -1] = ind.flow_angle[n, -1] + deviation_meanline
        
        if ind.blade_angle[n, -1] > np.deg2rad(40):
            ind.blade_angle[n, -1] = np.deg2rad(40)
            print(f"Clipping bladeline {n}")

    #######################################################################################################
    ### Blade Construction / Interpolation
    #######################################################################################################

    weight = 3
    shape = np.linspace(0, 1, n_blade)**weight # 0 at LE, 1 at TE, slow ramp early, fast ramp late
    for n in range(n_bladeline+1):
        ind.blade_angle[n] = ind.blade_angle[n,0] + (ind.blade_angle[n,-1] - ind.blade_angle[n,0]) * shape # rad/s - user-defined distribution of blade angles

    ds = ind.blade_length[-1]/n_blade # m - arc length differential
    ind.axial_height = np.trapezoid( np.sqrt(np.sin(ind.blade_angle[-1])**2 - np.gradient(ind.r_bladeline[-1], ds)**2 ), dx=ds)
    dz = ind.axial_height/(n_blade-1) # m - axial position differential

    ind.z = np.linspace(0, ind.axial_height, n_blade) # re-parameterize by uniform z spacing

    for n in range(n_bladeline+1):

        dr_ds = np.gradient(ind.r_bladeline[n], ds) # m/m - radial change over arc length
        dz_ds = np.sin(ind.blade_angle[n]) * np.sqrt(np.maximum(1 - dr_ds**2, 0)) # m/m - axial change over arc length
        
        z_of_s = cumulative_trapezoid(dz_ds, dx=ds, initial=0) # [m] - axial z as a function of arc length s
        z_of_s_scaled = z_of_s * (ind.axial_height / z_of_s[-1]) # [m] - stretch bladeline to terminate at the tip exit plane
        
        # Resample beta and r onto uniform z
        ind.blade_angle[n] = np.interp(ind.z, z_of_s_scaled, ind.blade_angle[n]) # [rad] - radial positions scaled and uniformly spaced in z
        r_z    = np.interp(ind.z, z_of_s_scaled, ind.r_bladeline[n]) # [m] - radial positions scaled and uniformly spaced in z
    
        ind.r_bladeline[n] = r_z    
        ind.z_bladeline[n] = -ind.z   # [m] - negative projects it downwards from the inlet plane
        ind.theta_bladeline[n] = -cumulative_trapezoid(np.sqrt(1 + np.gradient(r_z, dz)**2) / (r_z * np.tan(ind.blade_angle[n])), dx=dz, initial=0) # [m]

    ind.x_bladeline = ind.r_bladeline*np.cos(ind.theta_bladeline) # [m]
    ind.y_bladeline = ind.r_bladeline*np.sin(ind.theta_bladeline) # [m]

    


    #######################################################################################################
    ### Performance Prediction (approximation before CFD)
    #######################################################################################################

    ind.blockage = 0.98 # DUMMY
    ind.sweep_back_angle = 0 # DUMMY

    # Gulich empirical correlation - equation 7.23a
    x_sigma = (ind.flow_coeff**0.43 * (np.tan(ind.blade_angle[-1, 0])/ind.blockage)**0.75 * min(ind.length_pitch_ratio, 2.0)**-0.2 * 
            ind.hub_tip_ratio**0.3 * (ind.r_hub[-1]/ind.r_tip[-1])**0.1 * (1 + ind.sweep_back_angle*180/np.pi/100)**(-0.47) * (ind.blade_number/3)**0.24)

    sigma_3 = 1.26 * x_sigma
    NPSH_3_predicted = sigma_3 * ind.tip_speed**2 / (2 * engine.g)
    ind.suction_specific_speed = p.shaft_speed * np.sqrt(p.vdot) / (engine.g * NPSH_3_predicted)**0.75


def reparameterize_s_to_z(s_array, *quantities, n_z=None, z_max=None):
    dz_ds = quantities[0]
    other_quantities = quantities[1:]
    
    # Integrate dz/ds to get z(s)
    ds = s_array[1] - s_array[0]
    z_of_s = cumulative_trapezoid(dz_ds, dx=ds, initial=0)
    
    # Build target uniform-z grid
    if z_max is None:
        z_max = z_of_s[-1]
    if n_z is None:
        n_z = len(s_array)
    z = np.linspace(0, z_max, n_z)
    
    # Resample each quantity (which is a function of s) onto ind.z
    resampled = []
    for q in other_quantities:
        # q is sampled at s_array, which corresponds to z_of_s
        # Interpolate q vs z_of_s, then evaluate at ind.z
        q_resampled = np.interp(z, z_of_s, q)
        resampled.append(q_resampled)
    
    return z, *resampled
