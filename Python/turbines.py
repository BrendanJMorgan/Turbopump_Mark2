import numpy as np
import numpy.linalg as l
from scipy.interpolate import CubicSpline
from scipy.interpolate import RegularGridInterpolator
from engine_state import engine, tca, gg, turbine

import os

def turbines(t: turbine):

    # Data Loading --------------------------------------------------------
    path = os.path.dirname(__file__)
    blade_angle_interp = np.arange(0, 141, 1)                               # deg
    blade_width_interp = np.array([0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0])      # in
    incidence_loss = np.load(os.path.join(path, 'incidence_loss.npy'))                    # deg vs loss coeff
    mach_loss = np.load(os.path.join(path, 'mach_loss.npy'))                              # mach vs loss coeff
    impulse_turbine_efficiency = np.load(os.path.join(path, 'impulse_turbine_efficiency.npy'))    # v_ratio vs efficiency
    pitch_chord_zweifel = np.load(os.path.join(path, 'pitch_chord_zweifel.npy'))                    # blade angle deg vs pitch/chord ratio 
    expansion_loss = np.load(os.path.join(path, 'expansion_loss.npy'))                    # blade angle deg and blade width in vs loss coeff

    pitch_chord_zweifel = pitch_chord_zweifel[::-1,:]  # reverse to be ascending order

    # Nozzle Manifold --------------------------------------------------------
    t.p_manifold = gg.pc  # Pa - pressure at the nozzle manifold - assuming negligible dp from chamber to manifold
    t.p_spouting = engine.p_amb # Pa - assuming perfectly expanded flow - PROBABLY A BAD ASSUMPTION
    t.p_exit = engine.p_amb # Pa - assuming perfectly expanded flow - PROBABLY A BAD ASSUMPTION
    t.T_spouting = gg.Tc*(t.p_spouting/t.p_manifold)**((gg.gamma-1)/gg.gamma) # K - static temperature for perfectly expanded flow
    t.dens_spouting = gg.density*(t.p_spouting/t.p_manifold)**(1/gg.gamma)  # kg/m3 - static density for perfectly expanded flow

    t.c_tau = 1.50                              # unitless - initial guess  NEEDS UPDATING
    t.v_spouting_norm = gg.c_star * t.c_tau     # m/s
    t.mach_spouting = t.v_spouting_norm / np.sqrt(gg.gamma * gg.R * t.T_spouting)  # unitless

    t.A_throat_cum = gg.mdot * gg.c_star / t.p_manifold  # m2 - cumulative throat area on the nozzle plate
    nozzle_number = t.A_throat_cum / (np.pi/4*t.d_throat_nozzle**2) # unitless

    ## Velocities -------------------------------------------------
    # One stage impulse rotor - does not work for multistage and/or velocity driven turbines
    t.v_pitchline = t.r_pitchline*t.shaft_speed # m/s - tangential rotor speed at pitchline


    isentropic_v_ratio = t.v_pitchline/t.v_spouting_norm                      # unitless - "A single-row impulse stage delivers best performance at velocity ratios between 0.30 and 0.40"
    v_axial = np.sqrt(t.v_spouting_norm**2 - 4*t.v_pitchline**2)          # m/s
    v_spouting = np.array([2*t.v_pitchline, v_axial])                 # [m/s, m/s] - [tangential, axial]
    angle_nozzle = np.atan(v_spouting[0]/v_spouting[1])                         # rad
    w_in = np.array([t.v_pitchline, v_axial])             # m/s
    t.incidence_angle = np.atan(w_in[0]/w_in[1])    # rad - this is measured FROM THE VERTICAL
    t.blade_angle = t.incidence_angle                                   # rad

    w_mach_norm = l.norm(w_in) / np.sqrt(gg.gamma * gg.R * t.T_spouting)  # unitless - static temp different because not full pressure drop yet?

    # Blade Spacing -------------------------------------------------
    A_choked = gg.mdot*np.sqrt(gg.Tc)/t.p_manifold * np.sqrt(gg.R/gg.gamma) * ((gg.gamma+1)/2)**(0.5*(gg.gamma+1)/(gg.gamma-1)) # m2 - minimum throughput area needed to prevent gas backup
    A_exit = t.A_throat_cum * ( (2 + (gg.gamma - 1)*t.mach_spouting**2 ) / (gg.gamma + 1) ) ** ((gg.gamma + 1)/(2*(gg.gamma - 1))) / t.mach_spouting # m2 - cumulative area of all blade gaps in rotor

    if np.pi/2-t.blade_angle < 30*np.pi/180:
        raise ValueError(f"Rotor blade angle of {t.blade_angle*180/np.pi} deg is too agressive (Zweifel has no correlation)\n")
    elif np.pi/2-t.blade_angle > 80*np.pi/180:
        raise ValueError(f"Rotor blade angle of {t.blade_angle*180/np.pi} deg is too shallow (Zweifel has no correlation)\n")
    else:
        t.pitch_chord_ratio = CubicSpline(pitch_chord_zweifel[:,0], pitch_chord_zweifel[:,1])(np.pi/2-t.blade_angle*180/np.pi) # unitless - NASA Turbines 1974, figure 26, Zweifel impulse blades curve

    # zweifel_coeff = (0.90+1.15)/2 # unitless
    t.blade_chord = t.blade_width           # m - in an impulse rotor these are equal
    t.blade_pitch = t.blade_chord * t.pitch_chord_ratio                 # m - the ideal circumferential distance between neighboring blades
    t.blade_opening = t.blade_pitch * np.sin(t.blade_angle)   # m

    t.blade_count =  4*round(2*np.pi*t.r_pitchline/t.blade_pitch/4) # unitless - number of blades along the circumference of the rotor

    ## Energies and Losses -------------------------------------------------
    expansion_loss_interp = RegularGridInterpolator((blade_angle_interp, blade_width_interp), expansion_loss, bounds_error=False)
    expansion_loss_coeff = expansion_loss_interp(np.array([(np.pi - t.blade_angle - t.blade_angle)*180.0/np.pi, t.blade_width/0.0254]))
    kinetic_loss_coeff = 2*expansion_loss_coeff**2 - 1 # unitless, equation 3

    incidence_loss_coeff = CubicSpline(incidence_loss[:,0], incidence_loss[:,1])(t.incidence_angle*180/np.pi) # unitless - NASA Turbines 1976 fig. 17
    mach_loss_coeff = CubicSpline(mach_loss[:,0], mach_loss[:,1])(l.norm(w_mach_norm)) # unitless - NASA Turbines 1976 fig. 18
    # if np.sqrt(mach_relative[0]**2 + mach_relative[1]**2) <= 1:
    #     mach_loss_coeff = 1 # assuming that mach loss does not apply to a subsonic case?
    stage_efficiency = CubicSpline(impulse_turbine_efficiency[:,0], impulse_turbine_efficiency[:,1])(isentropic_v_ratio) # unitless
    clearance_loss_coeff = -1.63*t.tip_clearance/t.blade_length+1

    isentropic_enthalpy = gg.R*t.T_spouting*(gg.gamma/(gg.gamma-1))*(1-(t.p_exit/t.p_manifold)**((gg.gamma-1)/gg.gamma)) # J/kg - equation 1b, equation 2
    kinetic_enthalpy = 0.5*t.v_spouting_norm**2 # J/kg
    specific_work = isentropic_enthalpy*expansion_loss_coeff + kinetic_enthalpy*kinetic_loss_coeff*incidence_loss_coeff*mach_loss_coeff # J/kg - equation 5b

    turbine_efficiency_full = specific_work / (isentropic_enthalpy + kinetic_enthalpy) # unitless
    nozzle_arc = 2*np.pi*t.admission_fraction*t.r_pitchline # m
    w_out = w_in # this is a very bold assumption that NEEDS UPDATING
    wheel_ratio = w_out / w_in # unitless - ratio between outlet/inlet velocity of blades
    turbine_efficiency_partial = (1+wheel_ratio*(1-t.blade_pitch/(3*nozzle_arc)))/(1+wheel_ratio) * turbine_efficiency_full 
    - 1.539E-9*isentropic_v_ratio*t.dens_spouting*t.shaft_speed*(1-t.admission_fraction)/t.admission_fraction # unitless - equation 14

    turbine_efficiency = turbine_efficiency_partial*clearance_loss_coeff

    power_turbine = t.shaft_power/turbine_efficiency    # W - ideal power needed to meet pump requirements
    mdot_gg = 0.5 * power_turbine / (t.shaft_speed*t.r_pitchline)**2 # kg/s - mass flow rate through the gg / manifold / rotor - this is actually equation 10b simplified with impulse rotor assumptions
    torque_turbine = t.shaft_power / t.shaft_speed      # N*m

    # Blade Depth
    # if t.blade_length < 1.5*A_choked/(t.blade_count*t.blade_opening):
    #     raise ValueError("Rotor blades are too short and choking the flow")

    t.r_tip = t.r_pitchline + t.blade_length / 2    # m
    t.r_base = t.r_pitchline - t.blade_length / 2   # m

    # Flow Adjustments
    engine.mdot = tca.mdot + gg.mdot               # kg/s
    gg.fraction = gg.mdot / engine.mdot            # unitless
    gg.mdot_fuel = gg.mdot*(1/(1+gg.OF))           # kg/s - Fuel Mass Flow Rate
    gg.mdot_ox = gg.mdot*(gg.OF/(1+gg.OF))         # kg/s - Oxidizer Mass Flow Rate
    engine.mdot_fuel = tca.mdot_fuel+gg.mdot_fuel   # kg/s
    engine.mdot_ox = tca.mdot_ox+gg.mdot_ox         # kg/s
