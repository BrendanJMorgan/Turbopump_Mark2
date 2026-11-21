import numpy as np
import numpy.linalg as l
from scipy.interpolate import CubicSpline
from pathlib import Path
from engine_state import engine, tca, gg, turbine

from stators import stator_flow, stator_geometry
from rotors import rotor

# Load and build interpolators once at module load
impulse_turbine_efficiency = np.load(Path(__file__).parent / 'impulse_turbine_efficiency.npy')    # v_ratio vs efficiency

def turbines(t: turbine):
    s, r = t.stator, t.rotor

    t.p_exit = engine.p_amb     # Pa - assuming perfectly expanded flow - PROBABLY A BAD ASSUMPTION

    stator_flow(t) # manifold to spouting state

    #######################################################################################################
    ### Velocity Triangles - the interface between stator exit and rotor inlet
    #######################################################################################################
    t.r_pitchline = min(t.velocity_ratio*s.v_spouting_norm/t.shaft_speed, t.r_pitchline_max) # m - optimum velocity ratio, capped by manufacturing
    t.v_pitchline = t.r_pitchline*t.shaft_speed # m/s - tangential rotor speed at pitchline
    isentropic_v_ratio = t.v_pitchline/s.v_spouting_norm            # unitless - "A single-row impulse stage delivers best performance at velocity ratios between 0.30 and 0.40"
    if s.v_spouting_norm**2 < 4*t.v_pitchline**2:
        raise ValueError(f"Spouting velocity of {s.v_spouting_norm:0.2f} m/s is too low for the given shaft speed of {t.shaft_speed*30/np.pi:0.2f} rpm and radius of {t.r_pitchline*1000:0.1f} mm\n")
    v_axial = np.sqrt(s.v_spouting_norm**2 - 4*t.v_pitchline**2)    # m/s

    # One stage impulse rotor - does not work for multistage and/or velocity driven turbines
    s.v_spouting = np.array([2*t.v_pitchline, v_axial])   # [m/s, m/s] - [tangential, axial]
    s.angle_nozzle = np.atan(s.v_spouting[0]/s.v_spouting[1]) # rad
    r.w_in = np.array([t.v_pitchline, v_axial])           # m/s
    r.incidence_angle = np.atan(r.w_in[0]/r.w_in[1])        # rad - this is measured FROM THE VERTICAL
    r.blade_angle = r.incidence_angle                   # rad - measured from the vertical

    r.w_mach_norm = l.norm(r.w_in) / np.sqrt(gg.gamma_avg * gg.R_gas * s.T_spouting)  # unitless - static temp different because not full pressure drop yet?

    rotor(t) # blade geometry and loss coefficients
    t.stage_efficiency = CubicSpline(impulse_turbine_efficiency[:,0], impulse_turbine_efficiency[:,1])(isentropic_v_ratio) # unitless

    #######################################################################################################
    ### Stage Efficiency and Flow - iterated because the nozzle plate depends on the GG flow
    #######################################################################################################
    t.torque = t.shaft_power / t.shaft_speed      # N*m
    power_efficiency_old = 0
    t.power_efficiency = 1

    while abs(power_efficiency_old - t.power_efficiency) > 0.01:
        power_efficiency_old = t.power_efficiency

        t.power = t.shaft_power / t.power_efficiency    # W - ideal power needed to meet pump requirements
        gg.mdot = 0.5 * t.power / (t.shaft_speed*t.r_pitchline)**2 # kg/s - mass flow rate through the gg / manifold / rotor - this is actually equation 10b simplified with impulse rotor assumptions

        stator_geometry(t) # nozzle plate for this flow

        isentropic_enthalpy = gg.R_gas*gg.Tt*(gg.gamma_avg/(gg.gamma_avg-1))*(1-(t.p_exit/s.p_manifold)**((gg.gamma_avg-1)/gg.gamma_avg)) # J/kg - equation 1b, equation 2
        kinetic_enthalpy = 0.5*s.v_spouting_norm**2 # J/kg

        turbine_efficiency_full = specific_work / (isentropic_enthalpy + kinetic_enthalpy) # unitless
        wheel_ratio = r.w_out[0] / r.w_in[0] # unitless - ratio between outlet/inlet velocity of blades
        turbine_efficiency_partial = (1+wheel_ratio*(1-r.blade_pitch/(3*s.nozzle_arc)))/(1+wheel_ratio) * turbine_efficiency_full - 1.539E-9*isentropic_v_ratio*s.density_spouting*t.shaft_speed*(1-s.admission_fraction)/s.admission_fraction # unitless - equation 14

        t.power_efficiency = turbine_efficiency_partial*r.clearance_loss_coeff

    # Flow Adjustments
    engine.mdot = tca.mdot + gg.mdot               # kg/s
    gg.fraction = gg.mdot / engine.mdot            # unitless
    gg.mdot_fuel = gg.mdot*(1/(1+gg.mixture_ratio))           # kg/s - Fuel Mass Flow Rate
    gg.mdot_ox = gg.mdot*(gg.mixture_ratio/(1+gg.mixture_ratio))         # kg/s - Oxidizer Mass Flow Rate
    engine.mdot_fuel = tca.mdot_fuel+gg.mdot_fuel   # kg/s
    engine.mdot_ox = tca.mdot_ox+gg.mdot_ox         # kg/s
