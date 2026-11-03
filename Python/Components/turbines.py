import numpy as np
import numpy.linalg as l
from scipy.interpolate import CubicSpline
from scipy.interpolate import RegularGridInterpolator
from engine_state import engine, tca, gg, turbine
from scipy.interpolate import LinearNDInterpolator
from pathlib import Path

import os

# Load and build interpolators once at module load
incidence_loss = np.load(Path(__file__).parent / 'incidence_loss.npy')                    # deg vs loss coeff
mach_loss = np.load(Path(__file__).parent / 'mach_loss.npy')                              # mach vs loss coeff
impulse_turbine_efficiency = np.load(Path(__file__).parent / 'impulse_turbine_efficiency.npy')    # v_ratio vs efficiency
pitch_chord_zweifel = np.load(Path(__file__).parent / 'pitch_chord_zweifel.npy')                    # blade angle deg vs pitch/chord ratio 
expansion_loss = np.load(Path(__file__).parent / 'expansion_loss.npy')                    # blade angle deg and blade width in vs loss coeff

def turbines(t: turbine):

    # Nozzle Manifold --------------------------------------------------------
    t.p_manifold = gg.pc        # Pa - pressure at the nozzle manifold - assuming negligible dp from chamber to manifold
    t.p_spouting = engine.p_amb # Pa - assuming perfectly expanded flow - PROBABLY A BAD ASSUMPTION
    t.p_exit = engine.p_amb     # Pa - assuming perfectly expanded flow - PROBABLY A BAD ASSUMPTION
    t.T_spouting = gg.Tt*(t.p_spouting/t.p_manifold)**((gg.gamma_avg-1)/gg.gamma_avg)   # K - static temperature for perfectly expanded flow
    t.dens_spouting = t.p_spouting/(gg.R_gas*t.T_spouting)      # kg/m3 - static density for perfectly expanded flow (ideal gas at exit static state)

    t.c_tau = 1.50                              # unitless - initial guess  NEEDS UPDATING
    t.v_spouting_norm = gg.c_star * t.c_tau     # m/s
    t.mach_spouting = t.v_spouting_norm / np.sqrt(gg.gamma_avg * gg.R_gas * t.T_spouting)  # unitless
    t.r_pitchline = min(t.velocity_ratio*t.v_spouting_norm/t.shaft_speed, t.r_pitchline_max) # m - optimum velocity ratio, capped by manufacturing
    t.v_pitchline = t.r_pitchline*t.shaft_speed # m/s - tangential rotor speed at pitchline
    isentropic_v_ratio = t.v_pitchline/t.v_spouting_norm            # unitless - "A single-row impulse stage delivers best performance at velocity ratios between 0.30 and 0.40"
    if t.v_spouting_norm**2 < 4*t.v_pitchline**2:
        raise ValueError(f"Spouting velocity of {t.v_spouting_norm:0.2f} m/s is too low for the given shaft speed of {t.shaft_speed*30/np.pi:0.2f} rpm and radius of {t.r_pitchline*1000:0.1f} mm\n")
    v_axial = np.sqrt(t.v_spouting_norm**2 - 4*t.v_pitchline**2)    # m/s

    ## Velocities -------------------------------------------------
    # One stage impulse rotor - does not work for multistage and/or velocity driven turbines
    t.v_spouting = np.array([2*t.v_pitchline, v_axial])   # [m/s, m/s] - [tangential, axial]
    t.angle_nozzle = np.atan(t.v_spouting[0]/t.v_spouting[1]) # rad
    t.w_in = np.array([t.v_pitchline, v_axial])           # m/s
    t.incidence_angle = np.atan(t.w_in[0]/t.w_in[1])        # rad - this is measured FROM THE VERTICAL
    t.blade_angle = t.incidence_angle                   # rad - measured from the vertical

    t.w_mach_norm = l.norm(t.w_in) / np.sqrt(gg.gamma_avg * gg.R_gas * t.T_spouting)  # unitless - static temp different because not full pressure drop yet?

    t.torque = t.shaft_power / t.shaft_speed      # N*m
    power_efficiency_old = 0
    t.power_efficiency = 1

    while abs(power_efficiency_old - t.power_efficiency) > 0.01:
        power_efficiency_old = t.power_efficiency

        t.power = t.shaft_power / t.power_efficiency    # W - ideal power needed to meet pump requirements
        gg.mdot = 0.5 * t.power / (t.shaft_speed*t.r_pitchline)**2 # kg/s - mass flow rate through the gg / manifold / rotor - this is actually equation 10b simplified with impulse rotor assumptions

        t.A_throat_cum = gg.mdot * gg.c_star / t.p_manifold  # m2 - cumulative throat area on the nozzle plate

        # Blade Spacing -------------------------------------------------
        A_choked = gg.mdot*np.sqrt(gg.Tt)/t.p_manifold * np.sqrt(gg.R_gas/gg.gamma_avg) * ((gg.gamma_avg+1)/2)**(0.5*(gg.gamma_avg+1)/(gg.gamma_avg-1)) # m2 - minimum throughput area needed to prevent gas backup
        A_exit = t.A_throat_cum * ( (2 + (gg.gamma_avg - 1)*t.mach_spouting**2 ) / (gg.gamma_avg + 1) ) ** ((gg.gamma_avg + 1)/(2*(gg.gamma_avg - 1))) / t.mach_spouting # m2 - cumulative area of all blade gaps in rotor

        t.nozzle_number = int(np.ceil(t.A_throat_cum / (np.pi/4*t.d_throat_nozzle**2))) # unitless - whole nozzles, rounded up

        # Admission Geometry - drilled round nozzles -------------------------------------------------
        # Each exit meets the plate face as an ellipse, long axis along the arc (Aghaei & Tousi 2015, Fig. 5-6)
        t.d_exit_nozzle = t.d_throat_nozzle*np.sqrt(gg.Ae_At) # m - same expansion ratio as the GG nozzle
        t.nozzle_footprint = nozzle_footprint(t) # m - exit length along the arc, where the diverging cone meets the plate face
        t.nozzle_spacing = t.nozzle_spacing_ratio*t.nozzle_footprint # m - center to center along the pitchline
        t.nozzle_arc = (t.nozzle_number - 1)*t.nozzle_spacing + t.nozzle_footprint # m - active arc, first exit edge to last (Stenning sector loss)
        if t.nozzle_arc > 2*np.pi*t.r_pitchline:
            raise ValueError(f"{t.nozzle_number} turbine nozzles need a {t.nozzle_arc*1000:.0f} mm arc, more than the {2*np.pi*t.r_pitchline*1000:.0f} mm circumference")

        # Active fraction of stator-exit area (SP-290 eq. 8-23), overlap counted once (Aghaei & Tousi 2015 eqs. 1-7)
        area_projected = np.pi/4*t.d_exit_nozzle*t.nozzle_footprint # m2 - one exit ellipse on the plate face
        R, c = t.d_exit_nozzle/2, min(t.nozzle_spacing*t.d_exit_nozzle/t.nozzle_footprint, t.d_exit_nozzle) # m - ellipses scaled to circles along the arc
        area_lens = (2*R**2*np.arccos(c/(2*R)) - c/2*np.sqrt(4*R**2 - c**2))*t.nozzle_footprint/t.d_exit_nozzle # m2 - overlap of neighbouring exits
        t.admission_fraction = (t.nozzle_number*area_projected - (t.nozzle_number - 1)*area_lens)/(2*np.pi*t.r_pitchline*t.d_exit_nozzle) # unitless - epsilon

        if np.pi/2-t.blade_angle < 30*np.pi/180:
            raise ValueError(f"Rotor blade angle of {t.blade_angle*180/np.pi} deg is too aggressive (Zweifel has no correlation)\n")
        elif np.pi/2-t.blade_angle > 80*np.pi/180:
            raise ValueError(f"Rotor blade angle of {t.blade_angle*180/np.pi} deg is too shallow (Zweifel has no correlation)\nShaft Speed: {round(t.shaft_speed*30/np.pi,2)} rpm")
        else:
            pitch_chord_zweifel_ascending = pitch_chord_zweifel[::-1,:]  # reverse to be ascending order
            t.pitch_chord_ratio = CubicSpline(pitch_chord_zweifel_ascending[:,0], pitch_chord_zweifel_ascending[:,1])(-90+t.blade_angle*180/np.pi) # unitless - NASA Turbines 1974, figure 26, Zweifel impulse blades curve

        t.blade_chord = t.blade_width                           # m - in an impulse rotor these are equal
        t.blade_pitch = t.blade_chord * t.pitch_chord_ratio     # m - the ideal circumferential distance between neighboring blades
        t.blade_opening = t.blade_pitch * np.sin(t.blade_angle) # m

        t.blade_count =  4*round(2*np.pi*t.r_pitchline/t.blade_pitch/4) # unitless - number of blades along the circumference of the rotor

        #######################################################################################################
        ### Efficiency
        #######################################################################################################
        
        blade_angle_interp = np.arange(0, 141, 1)                               # deg
        blade_width_interp = np.array([0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0])      # in
        expansion_loss_interp = RegularGridInterpolator((blade_angle_interp, blade_width_interp), expansion_loss)
        
        t.blade_angle_in = t.blade_angle # rad - measured from the vertical
        t.blade_angle_out = t.blade_angle # rad - measured from the vertical
        expansion_loss_coeff = float(expansion_loss_interp(np.array([(t.blade_angle_in+t.blade_angle_out)*180.0/np.pi, t.blade_width/0.0254]))[0]) # unitless
        kinetic_loss_coeff = 2*expansion_loss_coeff**2 - 1 # unitless, equation 3

        incidence_loss_coeff = CubicSpline(incidence_loss[:,0], incidence_loss[:,1])(t.incidence_angle*180/np.pi) # unitless - NASA Turbines 1976 fig. 17
        t.mach_loss_coeff = CubicSpline(mach_loss[:,0], mach_loss[:,1])(l.norm(t.w_mach_norm)) # unitless - NASA Turbines 1976 fig. 18
        t.stage_efficiency = CubicSpline(impulse_turbine_efficiency[:,0], impulse_turbine_efficiency[:,1])(isentropic_v_ratio) # unitless
        clearance_loss_coeff = -1.63*t.tip_clearance/t.blade_length+1

        isentropic_enthalpy = gg.R_gas*gg.Tt*(gg.gamma_avg/(gg.gamma_avg-1))*(1-(t.p_exit/t.p_manifold)**((gg.gamma_avg-1)/gg.gamma_avg)) # J/kg - equation 1b, equation 2
        kinetic_enthalpy = 0.5*t.v_spouting_norm**2 # J/kg
        specific_work = isentropic_enthalpy*expansion_loss_coeff + kinetic_enthalpy*kinetic_loss_coeff*incidence_loss_coeff*t.mach_loss_coeff # J/kg - equation 5b

        turbine_efficiency_full = specific_work / (isentropic_enthalpy + kinetic_enthalpy) # unitless
        t.w_out = t.w_in # this is a very bold assumption that NEEDS UPDATING
        wheel_ratio = t.w_out[0] / t.w_in[0] # unitless - ratio between outlet/inlet velocity of blades
        turbine_efficiency_partial = (1+wheel_ratio*(1-t.blade_pitch/(3*t.nozzle_arc)))/(1+wheel_ratio) * turbine_efficiency_full - 1.539E-9*isentropic_v_ratio*t.dens_spouting*t.shaft_speed*(1-t.admission_fraction)/t.admission_fraction # unitless - equation 14

        t.power_efficiency = turbine_efficiency_partial*clearance_loss_coeff

    t.r_tip = t.r_pitchline + t.blade_length / 2    # m
    t.r_base = t.r_pitchline - t.blade_length / 2   # m

    # Flow Adjustments
    engine.mdot = tca.mdot + gg.mdot               # kg/s
    gg.fraction = gg.mdot / engine.mdot            # unitless
    gg.mdot_fuel = gg.mdot*(1/(1+gg.mixture_ratio))           # kg/s - Fuel Mass Flow Rate
    gg.mdot_ox = gg.mdot*(gg.mixture_ratio/(1+gg.mixture_ratio))         # kg/s - Oxidizer Mass Flow Rate
    engine.mdot_fuel = tca.mdot_fuel+gg.mdot_fuel   # kg/s
    engine.mdot_ox = tca.mdot_ox+gg.mdot_ox         # kg/s

def nozzle_footprint(t: turbine) -> float:
    """Tangential length of one tilted conical nozzle's cut through the plate face (cascade plane).
    Exit circle sits on the axis at the plate face; the two walls lean at tilt -/+ half-angle, so the cut is lopsided."""
    r_throat, r_exit = t.d_throat_nozzle/2, t.d_exit_nozzle/2 # m
    length_diverging = (r_exit - r_throat)/np.tan(t.diverge_angle_gg) # m
    along = np.array([np.sin(t.angle_nozzle), -np.cos(t.angle_nozzle)]) # unitless - down the jet, plate face at y = 0
    across = np.array([along[1], -along[0]]) # unitless
    throat = -length_diverging*along # m
    x = []
    for side in (1.0, -1.0):
        wall_start = throat + side*r_throat*across
        tilt = -side*t.diverge_angle_gg # rad - wall opens away from the axis
        wall = np.array([along[0]*np.cos(tilt) - along[1]*np.sin(tilt), along[0]*np.sin(tilt) + along[1]*np.cos(tilt)])
        x.append(wall_start[0] - wall_start[1]/wall[1]*wall[0]) # m - where this wall meets y = 0
    return abs(x[0] - x[1])


