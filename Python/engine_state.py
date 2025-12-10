import numpy as np
import CoolProp as cp
import rocketprops as rp
import math
from dataclasses import dataclass, asdict

class engine: # For parameters that don't fit neatly within other categories

    g = 9.80665         # m/s2

    thrust = 30000      # N - Thrust INCLUDE GG EXHAUST THRUST
    p_amb = 93900       # Pa - ambient pressure at 2100 feet elevation
    T_amb = 293         # K - Ambient Temperature
    
    fuel = 'JetA'       # 'JetA' 'C2H5OH(L)'
    fuelrp = 'RP1'      # rocketprop only has RP-1 not Jet-A
    oxidizer = 'LOX'    # 'O2(L)'
    proof = 1           # Proportion of fuel (remaining part is water)

    common_shaft = True # Whether the pumps share a common shaft

    mdot: float
    mdot_ox: float
    mdot_fuel: float

class tca: # Thrust Chamber Assembly

    pc = 50E5           # Pa - Chamber (Stagnation) Pressure
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
    r_chamber = 0.025                   # m - radius of combustion chamber

    pc: float
    mdot: float
    p_exit: float
    Tt: float
    c_star_ideal: float
    isp_ideal: float
    gamma: float
    R: float

class pump: # fluid is 'ox' or 'fuel'
    def __init__(self, fluid: str):
        self.fluid = fluid
        if self.fluid == 'ox':
            self.clocking = 1                        # 1 for counterclockwise and -1 for clockwise (looking down at pump inlet)
            self.shaft_speed = 20000*math.pi/30      # rad/s - angular velocity of the pump shaft, impeller, and inducers
            self.gear_efficiency = 1                 # unitless - spur gears ~= 0.95; common shaft = 1
            self.r_hub_impeller = 0.006                       # m
            self.r_shaft = 0.005                     # m - portion of shaft that is stainless steel
            self.impeller_thickness = 0.003          # m - thickness of impeller at the exit point, not including blades
            self.NPSH_margin = 1.5                   # unitless - Margin of extra net positive suction head (NPSH) to be provided by inducers to prevent cavitation
            self.blade_number_inducer = 4            # unitless - 3 or 4 is considered good
            self.clearance_radial_inducer = 2E-4     # m - radial clearance between inducer blades and housing cavity
            self.surface_roughness = 10E-6           # m - surface roughness of the additive material of the impellers
            self.hub_tip_ratio_inducer = 0.3         # unitless - ratio of hub radius to tip radius of inducer
        elif self.fluid == 'fuel':
            self.clocking = 1                        
            self.shaft_speed = 25000 * math.pi / 30  
            self.gear_efficiency = 1                             
            self.r_hub_impeller = 0.006               
            self.r_shaft = 0.005                               
            self.impeller_thickness = 0.003   
            self.eye_flow_coeff = 0.25           
            self.NPSH_margin = 1.5             
            self.blade_number_inducer = 4          
            self.clearance_radial_inducer = 2E-4  
            self.surface_roughness = 10E-6
            self.hub_tip_ratio_inducer = 0.3         # unitless - ratio of hub radius to tip radius of inducer
        else:
            raise ValueError("fluid must be 'ox' or 'fuel'")

class turbine:
    r_pitchline = 0.096                 # m
    gear_ratio = 1                            # unitless - higher makes for a faster, smaller turbine. 1 = common shaft
    d_throat_nozzle = 0.127 * 0.0254          # m - diameter of each nozzle leading off the manifold - mdot_gg is a direct function of this and nozzle_number
    diverge_angle_gg = 15 * math.pi / 180     # rad - half-cone divergence angle of nozzle plate
    blade_width = 0.5 * 0.0254          # m - real turbines seems to be about 0.5 inches or so
    admission_fraction = 20 / 360             # unitless - fraction of the nozzle plate circle that actually has nozzles
    blade_length = 0.25 * 0.0254        # m
    tip_clearance = 0.0005                    # m
    radius_leading = 0.002                    # m - blade leading edge fillet radius. NEEDS MORE RESEARCH

    power: float
    mdot: float
    p_in: float
    p_out: float
    Tt_in: float
    Tt_out: float
    density_turbine: float