import numpy as np
import CoolProp as cp
import rocketprops as rp
import math
from dataclasses import dataclass, asdict

class engine: # For parameters that don't fit neatly within other categories

    g = 9.80665         # m/s2

    thrust = 25000      # N - Thrust INCLUDE GG EXHAUST THRUST
    p_amb = 93900       # Pa - ambient pressure at 2100 feet elevation
    T_amb = 293         # K - Ambient Temperature
    
    fuel = 'JetA'       # 'JetA' 'C2H5OH(L)'
    fuelrp = 'RP1'      # rocketprop only has RP-1 not Jet-A
    oxidizer = 'LOX'    # 'O2(L)'
    proof = 1           # Proportion of fuel (remaining part is water)

    common_shaft = True # Whether the pumps share a common shaft

class tca: # Thrust Chamber Assembly

    # Chamber/Nozzle Geometry
    dx = 0.0001                                # m - position step 
    converge_angle = 45 * math.pi / 180       # rad - half-cone convergence angle of combustion chamber end
    diverge_angle = 15 * math.pi / 180        # rad - half-cone divergence angle of nozzle
    l_star = 1.5                                # m - characteristic chamber length
    rc_throat = 0.025                         # m - radius of curvature around the throat
    d2_chamber = 6.0*0.0254                        # m
    wall_thickness = 0.006                         # m % SLANT VS VERTICAL THICKNESS
    d1_chamber = d2_chamber - 2 * wall_thickness   # m
    r1_chamber = d1_chamber / 2               # m

    pc = 50E5           # Pa - Chamber (Stagnation) Pressure 
    mixture_ratio = 2.0 # Oxidizer/Fuel Ratio (by mass)
    c_star_efficiency = 0.85   # Characteristic Vel Efficiency, experimental
    c_tau_efficiency = 0.95    # Thrust Coefficient Efficiency, experimental

    p_injector_manifold = pc*1.25 # Pa - PLACEHOLDER UNTIL A PROPER INJECTOR MODULE IS WRITTEN
    compute_thermals = True                   # Whether to run regen thermal balance (tends to increase runtime significantly)

    class combustion_gas:
        p = None
        T = None    

    class film_coolant:
        fraction = np.array([0.05, 0.05])     # unitless - Fraction of the fuel mass flow dedicated to film cooling orifices - typically 3%-10% (Huzel and Huang)
        injection_x = [0, 0.25] # m - film cooling orifices around perimeter of injector and along bottom edge of chamber wall
        v_injection = 10            # m/s - combustion gas must have some initial velocity for injector film cooling to work mathematically
        injection_efficiency = 1
        class liquid:
            None
        class gas:
            None

    class regenerative_coolant:
        flow_direction = -1                      # 1 = forward flow (injector to nozzle), -1 = counter flow (nozzle to injector)

    class regenerative_jacket:
        n_pipe1 = 64                             # number of channels along barrel
        n_pipe2 = 32                             # number of channels near throat
        n_pipe3 = 64                             # number of channels along lower nozzle section
        gap_pipe = 1/8 * 0.0254                  # m - Gap between channels (fin thickness)
        h_pipe = 0.08 * 0.0254                    # m - coolant channel height
        merge_radius = 0.45*(5.5*0.0254-2*0.006) # m - when contour is below this radius, transition to n_pipe2
        cond = 401                               # W/m-K - thermal conductivity of copper



class gg: # Gas Generator
    stiffness = 0.25                    # Pa/Pa - guess
    fraction_guess = 0.05               # unitless - Fraction of total mass flow sent to the gas generator. Context: F1 = 0.030, J2 = 0.014 DO NOT CHANGE FROM 0.01 NOW
    mixture_ratio = 0.4                 # unitless - Oxidizer/Fuel Ratio - "[Most] operate at mixture ratios from 0.2 to 1.0, with hydrocarbons falling in the lower end, about 0.3" (NASA 1972)
    c_star_efficiency = 0.75            # m/s - characteristic velocity efficiency, experimental ANY PAPERS ON THIS?
    c_tau_efficiency = 0.96             # unitless - Thrust Coefficient Efficiency Factor
    orifice_number_fuel = 12            # 4 for impinging, 8 for showerhead film cooling on gg injector
    r_chamber = 0.025                   # m - radius of combustion chamber

class pump:
    def __init__(self, fluid: str):
        self.fluid = fluid
        self.surface_roughness = 10E-6
        self.shaft_speed = None

        if self.fluid == 'ox':
            self.clocking = 1                        # 1 for counterclockwise and -1 for clockwise (looking down at pump inlet)
            self.n_stages = 1
            self.inducer = inducer(fluid=self.fluid)
            self.impeller = [impeller(stage=i, fluid=self.fluid) for i in range(self.n_stages)]
            self.volute = volute(fluid=self.fluid)
            self.tank = tank(fluid=self.fluid)
        elif self.fluid == 'fuel':
            self.clocking = 1                        
            self.n_stages = 1
            self.inducer = inducer(fluid=self.fluid)
            self.impeller = [impeller(stage=i, fluid=self.fluid) for i in range(self.n_stages)]
            self.volute = volute(fluid=self.fluid)
            self.tank = tank(fluid=self.fluid)                        
        else:
            raise ValueError("fluid must be 'ox' or 'fuel'")
        
class inducer:
    def __init__(self, fluid: str):
        self.fluid = fluid
        if self.fluid == 'ox':
            self.flow_coeff = 0.1           # unitless
            self.blade_number = 3           # unitless
            self.clearance_radial = 2E-4    # m
            self.hub_tip_ratio = 0.3        # unitless - ratio of hub radius to tip radius of inducer
            self.hub_tip_ratio_out = 0.5    # unitless - ratio of hub radius to tip radius of inducer
            self.suction_margin = 1.5       # unitless - Margin of extra net positive suction head (NPSH)
            self.flow_margin = 1.2          # unitless - margin on flow rate to ensure inducer is not undersized
            self.blade_thickness = 0.001    # m - thickness of inducer blades
            self.tip_clearance = 0.001      # m
            self.length_pitch_ratio = 2.5   # 2.5 is the max empirical coefficient beyond which Gulich says there are negligible benefits
        elif self.fluid == 'fuel':
            self.flow_coeff = 0.1          
            self.blade_number = 3
            self.clearance_radial = 2E-4
            self.hub_tip_ratio = 0.3
            self.hub_tip_ratio_out = 0.5   
            self.suction_margin = 4
            self.flow_margin = 1.2      
            self.blade_thickness = 0.001  
            self.tip_clearance = 0.001
            self.length_pitch_ratio = 2.5

class impeller:
    def __init__(self, stage: int, fluid: str):
        self.fluid = fluid
        if self.fluid == 'ox':
            self.stage = stage
            self.hub_in_angle = -90 * np.pi / 180
            self.hub_out_angle = -10 * np.pi / 180
            self.shroud_in_angle = -90 * np.pi / 180
            self.shroud_out_angle = -0 * np.pi / 180
            self.thickness = 0.003
            self.NPSH_required = 20
            self.suction_margin = 1.5
            self.surface_roughness = 10E-6
        elif self.fluid == 'fuel':
            self.stage = stage
            self.hub_in_angle = -90 * np.pi / 180
            self.hub_out_angle = -10 * np.pi / 180
            self.shroud_in_angle = -90 * np.pi / 180
            self.shroud_out_angle = -0 * np.pi / 180
            self.thickness = 0.003
            self.NPSH_required = 20
            self.suction_margin = 1.5
            self.surface_roughness = 10E-6

class volute:
    def __init__(self, fluid: str):
        self.fluid = fluid
        if self.fluid == 'ox':
            self.dummy = 1.0               
        elif self.fluid == 'fuel':
            self.dummy = 1.0                     

class turbine: 
    r_pitchline = 0.08                      # m
    gear_ratio = 1                          # unitless - higher makes for a faster, smaller turbine. 1 = common shaft
    d_throat_nozzle = 0.127 * 0.0254        # m - diameter of each nozzle leading off the manifold - mdot_gg is a direct function of this and nozzle_number
    diverge_angle_gg = 15 * math.pi / 180   # rad - half-cone divergence angle of nozzle plate
    blade_width = 0.5 * 0.0254              # m - real turbines seems to be about 0.5 inches or so
    admission_fraction = 20 / 360           # unitless - fraction of the nozzle plate circle that actually has nozzles
    blade_length = 0.25 * 0.0254            # m
    tip_clearance = 0.0005                  # m
    radius_leading = 0.002                  # m - blade leading edge fillet radius. NEEDS MORE RESEARCH

class tank:
    def __init__(self, fluid: str):
        self.fluid = fluid
        if self.fluid == 'ox':
            self.p = 1.5E5                    # Pa
            self.T = 90                     # K
        elif self.fluid == 'fuel':
            self.p = 1E5                   # Pa
            self.T = 293                    # K