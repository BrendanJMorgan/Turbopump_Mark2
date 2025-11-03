import numpy as np
import CoolProp as cp
import rocketprops as rp
import math
from dataclasses import dataclass, asdict

class engine: # For parameters that don't fit neatly within other categories

    g = 9.80665         # m/s2

    thrust = 15000      # N - Thrust TODO: include gg exhaust thrust into total thrust
    p_amb = 93900       # Pa - ambient pressure at 2100 feet elevation
    T_amb = 293         # K - Ambient Temperature - lowest temperature for which 75/25 ethanol is solvable by CoolProp
    
    fuel = 'JetA'       # 'JetA' 'C2H5OH(L)'
    fuelrp = 'RP1'      # rocketprop only has RP-1 not Jet-A
    oxidizer = 'LOX'    # 'O2(L)'
    proof = 1           # Proportion of fuel (remaining part is water)

    common_shaft = True # Whether the pumps share a common shaft

    mdot: float
    mdot_ox: float
    mdot_fuel: float

class tca: # Thrust Chamber Assembly

    pc = 30E5           # Pa - Chamber (Stagnation) Pressure
    stiffness = 0.25    # Pa/Pa - guess   
    OF = 2.0            # Oxidizer/Fuel Ratio (by mass)
    c_star_eff = 0.85   # Characteristic Vel Efficiency, experimental
    c_tau_eff = 0.95    # Thrust Coefficient Efficiency Factor

    # Chamber/Nozzle Geometry
    dx = 0.001                                # m - position step 
    converge_angle = 45 * math.pi / 180       # rad - half-cone convergence angle of combustion chamber end
    diverge_angle = 15 * math.pi / 180        # rad - half-cone divergence angle of nozzle
    l_star = 1                                # m - characteristic chamber length UPDATE FOR JET-A
    rc_throat = 0.025                         # m - radius of curvature around the throat
    d2_chamber = 0.125                        # m
    wall_thickness = 0.006                         # m % SLANT VS VERTICAL THICKNESS
    d1_chamber = d2_chamber - 2 * wall_thickness   # m
    r1_chamber = d1_chamber / 2               # m

    # Coolant Channels
    n_pipe1 = 16                              # number of channels along barrel
    n_pipe2 = 16                              # number of channels near throat
    n_pipe3 = 16                              # number of channels along lower nozzle section
    gap_pipe = 1/4 * 0.0254                   # m - Gap between channels (fin thickness)
    h_pipe = 1/16 * 0.0254                    # m - coolant channel height
    merge_radius = 0.45 * d1_chamber          # m - when contour is below this radius, transition to n_pipe2
    flow_direction = -1                       # 1 = forward flow (injector to nozzle), -1 = counter flow (nozzle to injector)

    # Film Cooling
    film_fraction = 0.03                      # unitless - Fraction of the fuel mass flow dedicated to film cooling orifices - typically 3%-10% (Huzel and Huang)
    v_injection = 10                          # m/s - combustion gas must have some initial velocity for injector film cooling to work mathematically

    # Properties
    k_wall = 253                              # W/m-K - thermal conductivity of aluminum (at 571 C)
    FSu_design = 2.0                          # Minimum design factor of safety, ultimate
    shaft_tensile_strength = 510E6            # Pa - stainless steel 304, ultimate
    shaft_shear_strength = shaft_tensile_strength / math.sqrt(3)  # Pa - conservative Tresca criterion

    compute_thermals = False                  # Whether to run regen thermal balance (tends to increase runtime significantly)

    mdot: float
    mdot_fuel: float
    mdot_ox: float
    gamma_avg_nozzle: float
    Tt_cc: float
    c_star_ideal: float
    p_exit: float
    isp_ideal: float
    v_exhaust_ideal: float
    R_gas: float
    A_throat: float
    Ae_At: float   

class gg: # Gas Generator
    stiffness = 0.25                    # Pa/Pa - guess
    fraction_guess = 0.05               # unitless - Fraction of total mass flow sent to the gas generator. Context: F1 = 0.030, J2 = 0.014 DO NOT CHANGE FROM 0.01 NOW
    OF = 0.4                            # unitless - OF Ratio - "[Most] operate at mixture ratios from 0.2 to 1.0, with hydrocarbons falling in the lower end, about 0.3" (NASA 1972)
    c_star_eff = 0.75                   # m/s - characteristic velocity efficiency, experimental ANY PAPERS ON THIS?
    c_tau_eff = 0.96                    # unitless - Thrust Coefficient Efficiency Factor
    orifice_number_fuel = 12            # 4 for impinging, 8 for showerhead film cooling on gg injector

    mdot_gg: float
    p_exit: float
    Tt_gg: float
    c_star_ideal: float
    isp_ideal: float

@dataclass
class pump():
    fluid: str 
    
    name: str
    p_in: float                 # Pa
    p_out: float                # Pa
    density_pump: float         # kg/m^3
    mdot_pump: float            # kg/s
    pump_shaft_speed: float     # rad/s (or 1/s; match your source units)
    vapor_pressure_pump: float  # Pa
    surface_roughness: float    # m
    slip_factor: float  # unitless - iterative initial guess

    blockage: float             # unitless - iterative initial guess - 1 is completely open, pump handbook says 0.85 is typical?
    blade_count: int            # unitless - iterative initial guess
    head_pump: float            # m
    vdot_pump: float            # m3/s - Volumetric Flow Rate
    specific_speed: float       # unitless     
    specific_speed_imperial: float # Ns pseudo unitless
    outlet_flow_coeff: float    # unitless - phi_i in pump handbook
    shutoff_coeff_op: float     # unitless - psi_op
    shutoff_coeff: float        # 1.25 for volute designs and 1.31 for diffuser designs (kept as input)
    shutoff_head_coeff: float   # unitless - psi_
    pump_efficiency: float      # unitless
    shaft_power: float          # W
    shaft_power_theory: float   # W

    def __post_init__(self):
        if self.fluid == 'ox':
            clocking = 1                        # 1 for counterclockwise and -1 for clockwise (looking down at pump inlet)
            pump_shaft_speed = 30000*math.pi/30 # rad/s - angular velocity of the pump shaft, impeller, and inducers
            gear_efficiency = 1                 # unitless - spur gears ~= 0.95; common shaft = 1
            r_hub = 0.006                       # m
            r_shaft = 0.005                     # m - portion of shaft that is stainless steel
            impeller_thickness = 0.003          # m - thickness of impeller at the exit point, not including blades
            eye_flow_coeff = 0.25               # unitless - higher means smaller impeller but larger inducer
                # phi_e in pump handbook; pg 2.29: 0.2-0.3 for impellers, ~0.1 or less for inducers
            NPSH_margin = 1.5                   # unitless - Margin of extra net positive suction head (NPSH) to be provided by inducers to prevent cavitation
            blade_number_inducer = 4            # unitless - 3 or 4 is considered good
            clearance_radial_inducer = 2E-4     # m - radial clearance between inducer blades and housing cavity
            surface_roughness = 10E-6           # m - surface roughness of the additive material of the impellers
        elif self.fluid == 'fuel':
            clocking = -1                        
            pump_shaft_speed = 30000 * math.pi / 30  
            gear_efficiency = 1                             
            r_hub = 0.006               
            r_shaft = 0.005                               
            impeller_thickness = 0.003   
            eye_flow_coeff = 0.25           
            NPSH_margin = 1.5             
            blade_number_inducer = 4          
            clearance_radial_inducer = 2E-4  
            surface_roughness = 10E-6     
        else:
            raise ValueError("fluid must be 'ox' or 'fuel'")



class turbine:
    r_pitchline_rotor = 0.070                 # m
    gear_ratio = 1                            # unitless - higher makes for a faster, smaller turbine. 1 = common shaft
    d_throat_nozzle = 0.127 * 0.0254          # m - diameter of each nozzle leading off the manifold - mdot_gg is a direct function of this and nozzle_number
    diverge_angle_gg = 15 * math.pi / 180     # rad - half-cone divergence angle of nozzle plate
    blade_width_rotor = 0.5 * 0.0254          # m - real turbines seems to be about 0.5 inches or so
    admission_fraction = 20 / 360             # unitless - fraction of the nozzle plate circle that actually has nozzles
    blade_length_rotor = 0.25 * 0.0254        # m
    tip_clearance = 0.0005                    # m
    radius_leading = 0.002                    # m - blade leading edge fillet radius. NEEDS MORE RESEARCH

    power_turbine: float
    mdot_turbine: float
    p_in: float
    p_out: float
    Tt_in: float
    Tt_out: float
    density_turbine: float
