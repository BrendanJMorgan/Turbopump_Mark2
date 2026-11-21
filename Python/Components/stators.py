import numpy as np
from engine_state import engine, gg, turbine

def stator_flow(t: turbine):
    """Nozzle manifold to nozzle exit (spouting) state. Independent of flow rate, so it runs once before the triangles."""
    s = t.stator

    #######################################################################################################
    ### Nozzle Spouting Properties
    #######################################################################################################
    s.p_manifold = gg.pc        # Pa - pressure at the nozzle manifold - assuming negligible dp from chamber to manifold
    s.p_spouting = engine.p_amb # Pa - assuming perfectly expanded flow - PROBABLY A BAD ASSUMPTION
    s.T_manifold = gg.Tt        # K - stagnation temperature at the nozzle manifold - assumes no losses in the hot gas circuit
    s.T_spouting = gg.Tt*(s.p_spouting/s.p_manifold)**((gg.gamma_avg-1)/gg.gamma_avg)   # K - static temperature for perfectly expanded flow
    s.density_spouting = s.p_spouting/(gg.R_gas*s.T_spouting)      # kg/m3 - static density for perfectly expanded flow (ideal gas at exit static state)

    lambda_cone = 0.5*(1+np.cos(s.diverge_angle))         # Divergence correction factor for conical nozzle
    s.nozzle_thrust_coefficient_ideal = np.sqrt( (2*s.gamma_avg_nozzle**2/(s.gamma_avg_nozzle-1) * (2/(s.gamma_avg_nozzle+1))**((s.gamma_avg_nozzle+1)/(s.gamma_avg_nozzle-1)) * \
                    (1-(t.p_exit/s.pc)**((s.gamma_avg_nozzle-1)/s.gamma_avg_nozzle) ) ) ) + s.Ae_At*(t.p_exit-engine.p_amb)/s.pc # Ideal Thrust Coefficient
    s.nozzle_thrust_coefficient = s.nozzle_thrust_coefficient_ideal*lambda_cone*s.nozzle_thrust_coefficient_efficiency       # Thrust Coefficient
    s.v_spouting_norm = gg.c_star * s.nozzle_thrust_coefficient     # m/s - this is equaivalent to the enthalpy approach (NASA eqn 1, 5, 6)
    s.mach_spouting = s.v_spouting_norm / np.sqrt(gg.gamma_avg * gg.R_gas * s.T_spouting)  # unitless

def stator_geometry(t: turbine):
    """Nozzle plate sized for the current GG flow. Runs inside the efficiency loop because it depends on gg.mdot"""
    s = t.stator

    #######################################################################################################
    ### Nozzle Geometry
    #######################################################################################################
    s.A_throat_cum = gg.mdot * gg.c_star / s.p_manifold  # m2 - cumulative throat area on the nozzle plate
    s.nozzle_number = int(np.ceil(s.A_throat_cum / (np.pi/4*s.d_throat_nozzle**2))) # unitless - whole nozzles, rounded up

    # drilled round nozzles - each exit meets the plate face as an ellipse, long axis along the arc (Aghaei & Tousi 2015, Fig. 5-6)
    s.d_exit_nozzle = s.d_throat_nozzle*np.sqrt(gg.Ae_At) # m - same expansion ratio as the GG
    s.nozzle_footprint = nozzle_footprint(s) # m - exit length along the arc, where the diverging cone meets the plate face
    s.nozzle_spacing = s.nozzle_spacing_ratio*s.nozzle_footprint # m - center to center along the pitchline
    s.nozzle_arc = (s.nozzle_number - 1)*s.nozzle_spacing + s.nozzle_footprint # m - active arc, first exit edge to last (Stenning sector loss)
    if s.nozzle_arc > 2*np.pi*t.r_pitchline:
        raise ValueError(f"{s.nozzle_number} turbine nozzles need a {s.nozzle_arc*1000:.0f} mm arc, more than the {2*np.pi*t.r_pitchline*1000:.0f} mm circumference")

    # Active fraction of stator-exit area (SP-290 eq. 8-23), overlap counted once (Aghaei & Tousi 2015 eqs. 1-7)
    area_projected = np.pi/4*s.d_exit_nozzle*s.nozzle_footprint # m2 - one exit ellipse on the plate face
    R, c = s.d_exit_nozzle/2, min(s.nozzle_spacing*s.d_exit_nozzle/s.nozzle_footprint, s.d_exit_nozzle) # m - ellipses scaled to circles along the arc
    area_lens = (2*R**2*np.arccos(c/(2*R)) - c/2*np.sqrt(4*R**2 - c**2))*s.nozzle_footprint/s.d_exit_nozzle # m2 - overlap of neighbouring exits
    s.admission_fraction = (s.nozzle_number*area_projected - (s.nozzle_number - 1)*area_lens)/(2*np.pi*t.r_pitchline*s.d_exit_nozzle) # unitless - epsilon

    A_choked = gg.mdot*np.sqrt(gg.Tt)/s.p_manifold * np.sqrt(gg.R_gas/gg.gamma_avg) * ((gg.gamma_avg+1)/2)**(0.5*(gg.gamma_avg+1)/(gg.gamma_avg-1)) # m2 - minimum throughput area needed to prevent gas backup
    A_exit = s.A_throat_cum * ( (2 + (gg.gamma_avg - 1)*s.mach_spouting**2 ) / (gg.gamma_avg + 1) ) ** ((gg.gamma_avg + 1)/(2*(gg.gamma_avg - 1))) / s.mach_spouting # m2 - cumulative area of all blade gaps in rotor
    if A_exit < A_choked:
        raise ValueError(f"Cumulative blade gap area {A_exit:.4f} m2 is less than the minimum required area {A_choked:.4f} m2")

def nozzle_footprint(s: turbine.stator) -> float:
    """Tangential length of one tilted conical nozzle's cut through the plate face (cascade plane).
    Exit circle sits on the axis at the plate face; the two walls lean at tilt -/+ half-angle, so the cut is lopsided."""
    r_throat, r_exit = s.d_throat_nozzle/2, s.d_exit_nozzle/2 # m
    length_diverging = (r_exit - r_throat)/np.tan(s.diverge_angle) # m
    along = np.array([np.sin(s.angle_nozzle), -np.cos(s.angle_nozzle)]) # unitless - down the jet, plate face at y = 0
    across = np.array([along[1], -along[0]]) # unitless
    throat = -length_diverging*along # m
    x = []
    for side in (1.0, -1.0):
        wall_start = throat + side*r_throat*across
        tilt = -side*s.diverge_angle # rad - wall opens away from the axis
        wall = np.array([along[0]*np.cos(tilt) - along[1]*np.sin(tilt), along[0]*np.sin(tilt) + along[1]*np.cos(tilt)])
        x.append(wall_start[0] - wall_start[1]/wall[1]*wall[0]) # m - where this wall meets y = 0
    return abs(x[0] - x[1])
