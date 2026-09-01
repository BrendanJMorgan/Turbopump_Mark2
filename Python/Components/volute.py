import numpy as np
from scipy.interpolate import LinearNDInterpolator
from pathlib import Path

from engine_state import engine, pump

# Load and build interpolator once at module load
diffuser_data = np.load(Path(__file__).parent / 'diffuser_loss.npy')
diffuser_interp = LinearNDInterpolator(np.log(diffuser_data[:, :2]), diffuser_data[:, 2])

def volute(p) -> None:

    #######################################################################################################
    ### High Level Parameters
    #######################################################################################################
    vol = p.volute
    imp = p.impeller[0]  # Currently only single stage pumps are supported

    vol.r_tongue = imp.r_outlet * (1.03 + 0.0025 * imp.specific_speed + 7e-8 * p.density * p.head_rise)  # m - Gulich table 10.2, pg 729
    vol.t_tongue = 0.04 * imp.r_outlet # m - Gulich pg 493
    vol.r_tongue_stag = vol.r_tongue + 0.5 * vol.t_tongue  # m - stagnation radius

    p.vdot_seals = 0 # PLACEHOLDER
    vol.vdot = p.vdot - p.vdot_seals

    # Eq. T7.7.7 - circular throat (or equivalent circle for arbitrary section)
    vol.r_throat_inner = vol.r_tongue_stag # m - this assumption is mostly for manufacturability (meanline on constant z). 
        # Gulich uses stagnation point instead of r_tongue_outer because the outer tongue flow doesn't full stagnate (approximation())
    X_Sp = (vol.vdot / (np.pi * imp.v_tangential_outlet * imp.r_outlet)) * vol.wrap_angle / (2.0 * np.pi)  # m
    vol.d_throat = X_Sp + 2*np.sqrt(vol.r_throat_inner * X_Sp)  # m - hydraulic diameter of throat; Gulich 7.7.7
    vol.A_throat = np.pi/4 * vol.d_throat ** 2  # m2

    vol.r_throat_outer = vol.r_throat_inner + vol.d_throat  # m
    vol.c_throat = vol.vdot / vol.A_throat  # m/s - throat velocity (= c_3q)

    #######################################################################################################
    ### Contour Construction
    #######################################################################################################
    n_theta = 128
    vol.theta = np.linspace(0, vol.wrap_angle, n_theta) # rad - azimuth
    dtheta = 2*np.pi/n_theta # rad - azimuth increment
    vol.A_section = np.zeros_like(vol.theta) # m2 - cross-sectional area along the contour

    vol.r_inner = np.linspace(vol.r_tongue, vol.r_throat_inner, n_theta) # m - inner radius along the contour
    
    # Design volutes for nq < 25 to 35 for constant velocity as per Eqs. (7.26 to 7.27)
    # Design volutes for nq > 25 to 35 according to the conservation of angular
    # momentum as per Table 7.7.
    if imp.specific_speed < 30.0/52.9186:  # threshold in Gulich's 25-35 band
        vol.design_method = "constant_velocity"
        vol.A_section = vol.A_throat * vol.theta / vol.wrap_angle  # m2 - gulich 7.27
 
    else:
        vol.design_method = "angular_momentum"
        
        X = (vol.vdot / (np.pi * imp.v_tangential_outlet * imp.r_outlet)) * (vol.theta / (2.0 * np.pi))  # m
        d_section = X + np.sqrt(2.0 * vol.r_throat_inner * X)  # m - Gulich 7.7.7 generalized to partial wrap angle theta
        vol.A_section = np.pi/4 * d_section ** 2  # m^2

    vol.d_hydraulic_section = np.sqrt(4.0 * vol.A_section / np.pi)  # m - hydraulic diameter of section

    vol.r_meanline = vol.r_inner + 0.5 * vol.d_hydraulic_section  # m - meanline radius of section

    n_phi = 32
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi)        # rad - angle around the circular cross section
    vol.z_meanline = imp.meanline_curve[-1, 1]              # m - axial center (impeller discharge meanline)

    # Cross-section point in the meridional (r, z) plane, broadcast over (theta, phi)
    r_wall = vol.r_meanline[:, None] + 0.5 * vol.d_hydraulic_section[:, None] * np.cos(phi)[None, :]  # m
    z_wall = vol.z_meanline + 0.5 * vol.d_hydraulic_section[:, None] * np.sin(phi)[None, :]                 # m

    # Sweep around the axis with clocking sign
    vol.contour = np.stack([
        r_wall * np.cos(p.clocking * vol.theta[:, None]),   # x
        r_wall * np.sin(p.clocking * vol.theta[:, None]),   # y
        z_wall,                                             # z
    ], axis=-1)  # [m,m,m] - shape (n_theta, n_phi, 3)

    #######################################################################################################
    ### Losses
    #######################################################################################################

    vol.c_section = imp.v_tangential_outlet * imp.r_outlet / vol.r_meanline  # m/s
    vol.deceleration_ratio = vol.c_throat / imp.v_merid[-1] # unitless

    # Friction loss
    reynolds = vol.c_throat * vol.d_throat / (p.viscosity/p.density) # unitless
    roughness_ratio = p.surface_roughness / vol.d_throat
    cf_friction = 0.0625 / (np.log10(roughness_ratio / 3.7 + 5.74 / reynolds ** 0.9)) ** 2 # COMES FROM IMPELLER TABLE? Friction coefficient from flat-plate model (Eq. T3.8.3)
 
    # Wetted perimeter and arc-length increments per section
    perimeter = np.pi*vol.d_hydraulic_section  # m - circular
    dA_wetted = perimeter * vol.r_meanline * dtheta  # m2 per section
 
    power_friction_section = 0.5 * p.density * np.sum((cf_friction + 0.0015) * vol.c_section**3 * dA_wetted)  # W - Gulich 3.8.21 - power dissipated via friction
    vol.friction_loss = 2.0 * np.sum(power_friction_section) / (p.density * imp.vdot * (p.shaft_speed * imp.r_outlet) ** 2) # friction loss - Gulich 3.8.21

    # Discharge nozzle deceleration loss
    p.d_mainline = 0.5*0.0254 # HARDCODED PLACEHOLDER
    p.A_mainline = np.pi/4*p.d_mainline**2 # m2 - cross-sectional area of mainline
    vol.discharge_length = 2*0.0254 # m
    vol.discharge_cp = discharge_cp(vol.discharge_length*2/vol.d_hydraulic_section[-1], p.A_mainline/vol.A_section[-1]) # pressure recovery coefficient for discharge nozzle - Gulich Fig. 1.19
    vol.nozzle_loss = vol.c_throat**2 / (p.shaft_speed*imp.r_outlet)**2 * (1.0 - vol.discharge_cp - 1.0 / (p.A_mainline/vol.A_throat)**2) # diffuser nozzle loss - Gulich 3.8.22

    vol.w_inlet = 1.5*imp.w_outlet
    vol.shock_loss = imp.flow_coeff[-1]**2 * (imp.blockage[-1] - imp.w_outlet/vol.w_inlet)**2 # vaneless diffuser shock loss - Gulich 3.8.24
    c = 0.63 # 0.65 for diffusers
    x0 = 0.87 # 0.93 for diffusers
    vol.off_bep_loss = c*x0**2 - 2*x0*c*vol.c_throat/imp.c_bladed[-1] + c*(vol.c_throat/imp.c_bladed[-1])**2 # off-BEP loss - Gulich 3.8.25 - empirical curvefit
        
    vol.total_loss = vol.friction_loss + vol.nozzle_loss + vol.shock_loss + (vol.off_bep_loss) # total loss - Gulich 3.8.22a

    #######################################################################################################
    ### Curve Generation
    #######################################################################################################
 
    # "...considerable freedom to configure the cross section without risking major losses ... flat cross sections result in less intense secondary flow than circular cross
    # sections, thereby generating fewer losses ... a ratio of width B to height H in the range of B/H = 2 to 3 may be considered the optimum." - Gulich pg 492

def discharge_cp(L_R1: float, A_R: float) -> float:
    """
    Pressure recovery coefficient c_p for a conical diffuser.
    Interpolated from digitized Gulich Fig. 1.19 (Sovran-Klomp data).
    """
    query = np.log([L_R1, A_R])
    pressure = float(diffuser_interp(query))
    if np.isnan(pressure):
        raise ValueError(f"(L/R1={L_R1}, A_R={A_R}) outside Fig. 1.19 data range")
    return pressure