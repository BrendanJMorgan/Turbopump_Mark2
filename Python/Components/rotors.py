import numpy as np
import numpy.linalg as l
from scipy.interpolate import CubicSpline
from scipy.interpolate import RegularGridInterpolator
from pathlib import Path
from engine_state import turbine

# Load and build interpolators once at module load
incidence_loss = np.load(Path(__file__).parent / 'incidence_loss.npy')                    # deg vs loss coeff
mach_loss = np.load(Path(__file__).parent / 'mach_loss.npy')                              # mach vs loss coeff
pitch_chord_zweifel = np.load(Path(__file__).parent / 'pitch_chord_zweifel.npy')                    # blade angle deg vs pitch/chord ratio
expansion_loss = np.load(Path(__file__).parent / 'expansion_loss.npy')                    # blade angle deg and blade width in vs loss coeff

def rotor(t: turbine):
    """Blade row geometry and loss coefficients. Needs the velocity triangles; independent of flow rate."""
    r = t.rotor

    #######################################################################################################
    ### Blade Geometry
    #######################################################################################################
    if np.pi/2-r.blade_angle < 30*np.pi/180:
        raise ValueError(f"Rotor blade angle of {r.blade_angle*180/np.pi} deg is too aggressive (Zweifel has no correlation)\n")
    elif np.pi/2-r.blade_angle > 80*np.pi/180:
        raise ValueError(f"Rotor blade angle of {r.blade_angle*180/np.pi} deg is too shallow (Zweifel has no correlation)\nShaft Speed: {round(t.shaft_speed*30/np.pi,2)} rpm")
    else:
        pitch_chord_zweifel_ascending = pitch_chord_zweifel[::-1,:]  # reverse to be ascending order
        r.pitch_chord_ratio = CubicSpline(pitch_chord_zweifel_ascending[:,0], pitch_chord_zweifel_ascending[:,1])(-90+r.blade_angle*180/np.pi) # unitless - NASA Turbines 1974, figure 26, Zweifel impulse blades curve

    r.blade_chord = r.blade_width                           # m - in an impulse rotor these are equal
    r.blade_pitch = r.blade_chord * r.pitch_chord_ratio     # m - the ideal circumferential distance between neighboring blades
    r.blade_opening = r.blade_pitch * np.sin(r.blade_angle) # m

    r.blade_count =  4*round(2*np.pi*t.r_pitchline/r.blade_pitch/4) # unitless - number of blades along the circumference of the rotor

    r.r_tip = t.r_pitchline + r.blade_length / 2    # m
    r.r_base = t.r_pitchline - r.blade_length / 2   # m

    #######################################################################################################
    ### Loss Coefficients
    #######################################################################################################
    blade_angle_interp = np.arange(0, 141, 1)                               # deg
    blade_width_interp = np.array([0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0])      # in
    expansion_loss_interp = RegularGridInterpolator((blade_angle_interp, blade_width_interp), expansion_loss)

    r.blade_angle_in = r.blade_angle # rad - measured from the vertical
    r.blade_angle_out = r.blade_angle # rad - measured from the vertical
    r.expansion_loss_coeff = float(expansion_loss_interp(np.array([(r.blade_angle_in+r.blade_angle_out)*180.0/np.pi, r.blade_width/0.0254]))[0]) # unitless
    r.kinetic_loss_coeff = 2*r.expansion_loss_coeff**2 - 1 # unitless, equation 3

    r.incidence_loss_coeff = CubicSpline(incidence_loss[:,0], incidence_loss[:,1])(r.incidence_angle*180/np.pi) # unitless - NASA Turbines 1976 fig. 17
    r.mach_loss_coeff = CubicSpline(mach_loss[:,0], mach_loss[:,1])(l.norm(r.w_mach_norm)) # unitless - NASA Turbines 1976 fig. 18
    r.clearance_loss_coeff = -1.63*r.tip_clearance/r.blade_length+1

    r.w_out = r.w_in # this is a very bold assumption that NEEDS UPDATING
