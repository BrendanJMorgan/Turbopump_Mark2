import os
import numpy as np
import yaml
from pint import UnitRegistry
ureg = UnitRegistry()

def baseunit(value):
    if isinstance(value, dict) and 'value' in value and 'unit' in value:
        return ureg.Quantity(value['value'], value['unit']).to_base_units().magnitude
    else:
        raise ValueError("Unsupported value type for conversion")

CONFIG_PATH = os.path.join(os.path.dirname(__file__), 'config.yaml') # Load primitive inputs from config.yaml.

with open(CONFIG_PATH, 'r') as _f:
    config = yaml.safe_load(_f)

class engine:  # For parameters that don't fit neatly within other categories
    c = config['engine']
    g               = c['g']                # m/s2
    thrust          = c['thrust']           # N - Thrust INCLUDE GG EXHAUST THRUST
    p_amb           = c['p_amb']            # Pa - ambient pressure at 2100 feet elevation
    T_amb           = c['T_amb']            # K - Ambient Temperature
    fuel            = c['fuel']             # 'JetA' 'C2H5OH(L)'
    fuelrp          = c['fuelrp']           # rocketprop only has RP-1 not Jet-A
    oxidizer        = c['oxidizer']         # 'O2(L)'
    proof           = c['proof']            # Proportion of fuel (remaining part is water)
    shaft_speed_cap = c['shaft_speed_cap']  # rad/s - max shaft speed if not suction limited
    common_shaft    = c['common_shaft']     # Whether the pumps share a common shaft
    del c

class tca:  # Thrust Chamber Assembly
    c = config['tca']

    # Chamber/Nozzle Geometry
    dx             = c['dx']              # m - position step
    converge_angle = baseunit(c['converge_angle'])  # rad - half-cone convergence angle of chamber end
    diverge_angle  = baseunit(c['diverge_angle'])   # rad - half-cone divergence angle of nozzle
    l_star         = c['l_star']          # m - characteristic chamber length
    rc_throat      = c['rc_throat']       # m - radius of curvature around the throat
    d2_chamber     = baseunit(c['d2_chamber'])  # m
    wall_thickness = c['wall_thickness']  # m - SLANT VS VERTICAL THICKNESS
    d1_chamber     = d2_chamber - 2 * wall_thickness  # m  (derived)
    r1_chamber     = d1_chamber / 2                   # m  (derived)

    pc                = c['pc']                 # Pa - Chamber (Stagnation) Pressure
    mixture_ratio     = c['mixture_ratio']      # Oxidizer/Fuel Ratio (by mass)
    c_star_efficiency = c['c_star_efficiency']  # Characteristic Vel Efficiency, experimental
    c_tau_efficiency  = c['c_tau_efficiency']   # Thrust Coefficient Efficiency, experimental

    p_injector_manifold = pc * 1.25  # Pa - PLACEHOLDER UNTIL A PROPER INJECTOR MODULE IS WRITTEN (derived)
    compute_thermals    = c['compute_thermals']  # Whether to run regen thermal balance

    class combustion_gas:
        p = None
        T = None

    class film_coolant:
        c = config['tca']['film_coolant']
        fraction             = np.array(c['fraction'])  # unitless - fuel mass fraction per film-cooling row
        injection_x          = c['injection_x']         # m - film cooling orifice positions
        v_injection          = c['v_injection']         # m/s - initial gas velocity for film model
        injection_efficiency = c['injection_efficiency']
        del c

        class liquid:
            None

        class gas:
            None

    class regenerative_coolant:
        flow_direction = config['tca']['regenerative_coolant']['flow_direction']  # 1 = forward, -1 = counter

    class regenerative_jacket:
        c = config['tca']['regenerative_jacket']
        n_pipe1      = c['n_pipe1']       # number of channels along barrel
        n_pipe2      = c['n_pipe2']       # number of channels near throat
        n_pipe3      = c['n_pipe3']       # number of channels along lower nozzle section
        gap_pipe     = baseunit(c['gap_pipe'])  # m - gap between channels (fin thickness)
        h_pipe       = baseunit(c['h_pipe'])    # m - coolant channel height
        merge_radius = c['merge_radius']  # m - contour radius to transition to n_pipe2
        cond         = c['cond']          # W/m-K - thermal conductivity of copper
        del c

    del c

class gg:  # Gas Generator
    c = config['gg']
    stiffness           = c['stiffness']            # Pa/Pa - guess
    fraction_guess      = c['fraction_guess']       # unitless - fraction of total mass flow to gg
    mixture_ratio       = c['mixture_ratio']        # unitless - Oxidizer/Fuel Ratio
    c_star_efficiency   = c['c_star_efficiency']    # characteristic velocity efficiency, experimental
    c_tau_efficiency    = c['c_tau_efficiency']     # unitless - Thrust Coefficient Efficiency Factor
    orifice_number_fuel = c['orifice_number_fuel']  # 4 impinging, 8 showerhead film cooling
    r_chamber           = c['r_chamber']            # m - radius of combustion chamber
    del c

class turbine:
    c = config['turbine']
    r_pitchline        = baseunit(c['r_pitchline'])      # m
    gear_ratio         = c['gear_ratio']         # unitless - higher = faster/smaller turbine; 1 = common shaft
    d_throat_nozzle    = baseunit(c['d_throat_nozzle'])  # m - diameter of each nozzle off the manifold
    diverge_angle_gg   = baseunit(c['diverge_angle_gg']) # rad - half-cone divergence angle of nozzle plate
    blade_width        = baseunit(c['blade_width'])      # m
    admission_fraction = c['admission_fraction']  # unitless - fraction of nozzle-plate circle with nozzles
    blade_length       = baseunit(c['blade_length'])     # m
    tip_clearance      = c['tip_clearance']      # m
    radius_leading     = c['radius_leading']     # m - blade leading edge fillet radius. NEEDS MORE RESEARCH
    del c

class pump:
    def __init__(self, fluid: str):
        self.fluid = fluid
        c = config['pump'][fluid]  # KeyError here if fluid is not configured

        self.surface_roughness = c['surface_roughness']
        self.shaft_speed = None
        self.clocking = c['clocking']    # 1 for counterclockwise and -1 for clockwise (looking down at inlet)
        self.n_stages = c['n_stages']
        self.inducer = inducer(fluid=self.fluid)
        self.impeller = [impeller(stage=i, fluid=self.fluid) for i in range(self.n_stages)]
        self.volute = volute(fluid=self.fluid)
        self.tank = tank(fluid=self.fluid)
        if 'plumbing_loss' in c:
            self.plumbing_loss = baseunit(c['plumbing_loss'])  # Pa - fuel side only

class inducer:
    def __init__(self, fluid: str):
        self.fluid = fluid
        c = config['pump'][fluid]['inducer']
        self.flow_coeff            = c['flow_coeff']             # unitless
        self.blade_number          = c['blade_number']           # unitless
        self.clearance_radial      = c['clearance_radial']       # m
        self.hub_tip_ratio         = c['hub_tip_ratio']          # unitless - inlet hub/tip radius ratio
        self.hub_tip_ratio_out     = c['hub_tip_ratio_out']      # unitless - outlet hub/tip radius ratio
        self.suction_margin_factor = c['suction_margin_factor']  # unitless - margin of extra NPSH
        self.flow_margin_factor    = c['flow_margin_factor']     # unitless - margin on flow rate
        self.blade_thickness       = c['blade_thickness']        # m - thickness of inducer blades
        self.tip_clearance         = c['tip_clearance']          # m
        self.length_pitch_ratio    = c['length_pitch_ratio']     # max empirical value per Gulich
        self.sweepback             = baseunit(c['sweepback'])    # rad - blade leading-edge sweepback

class impeller:
    def __init__(self, stage: int, fluid: str):
        self.fluid = fluid
        self.blade = blade(fluid=fluid)
        c = config['pump'][fluid]['impeller']
        self.stage             = stage
        self.NPSH_required     = c['NPSH_required']
        self.hub_in_angle      = baseunit(c['hub_in_angle'])      # rad
        self.hub_out_angle     = baseunit(c['hub_out_angle'])     # rad
        self.shroud_in_angle   = baseunit(c['shroud_in_angle'])   # rad
        self.shroud_out_angle  = baseunit(c['shroud_out_angle'])  # rad
        self.thickness         = c['thickness']
        self.suction_margin    = c['suction_margin']
        self.surface_roughness = c['surface_roughness']

class blade:
    def __init__(self, fluid: str):
        self.fluid = fluid

class volute:
    def __init__(self, fluid: str):
        self.fluid = fluid
        c = config['pump'][fluid]['volute']
        self.dummy      = c['dummy']
        self.wrap_angle = baseunit(c['wrap_angle'])  # rad - single volute is 2*pi, double volute is pi, etc.

class tank:
    def __init__(self, fluid: str):
        self.fluid = fluid
        c = config['pump'][fluid]['tank']
        self.p = c['p']  # Pa
        self.T = c['T']  # K

class link:
    def __init__(self):
      dummy = None  

class motor:
    def __init__(self, fluid: str, poles: int, ang_frequency: float, power: float, r_max: float):
        self.fluid = fluid
        self.poles = poles
        self.ang_frequency = ang_frequency
        self.power = power
        self.r_max = r_max
        c = config['machine']['motor']
        self.reactance_perunit   = c['reactance_perunit']          # unitless - per-unit synchronous reactance
        self.resistance_phase    = c['resistance_phase']      # ohm - per-phase winding resistance
        self.core_loss_fraction  = c['core_loss_fraction']    # unitless - crude core loss / throughput
        self.remanence           = c['remanence']             # T - magnet B_r (SmCo)
        self.relative_recoil_permeability = c['relative_recoil_permeability']   # unitless - mu_r
        self.magnet_length       = c['magnet_length']         # m - axial magnet thickness
        self.airgap              = c['airgap']                # m - mechanical airgap
        self.r_inner             = c['r_inner']               # m - annulus ID (centre bearing bore)
        self.shear_allowable           = c['shear_allowable']             # Pa - airgap shear ceiling

class alternator:
    def __init__(self, fluid: str, poles: int, ang_frequency: float, r_max: float):
        self.fluid = fluid
        self.poles = poles
        self.ang_frequency = ang_frequency
        self.r_max = r_max
        c = config['machine']['alternator']
        self.reactance_perunit   = c['reactance_perunit']          # unitless - per-unit synchronous reactance
        self.resistance_phase    = c['resistance_phase']      # ohm - per-phase winding resistance
        self.core_loss_fraction  = c['core_loss_fraction']    # unitless
        self.remanence           = c['remanence']             # T - magnet B_r (SmCo)
        self.relative_recoil_permeability = c['relative_recoil_permeability']   # unitless - mu_r
        self.magnet_length       = c['magnet_length']         # m
        self.airgap              = c['airgap']                # m
        self.r_inner             = c['r_inner']               # m
        self.shear_allowable           = c['shear_allowable']             # Pa

class cable:
    def __init__(self, frequency: float):
        self.frequency = frequency  # Hz - electrical frequency of the link
        c = config['machine']['cable']
        self.length           = c['length']            # m - conductor run length
        self.area             = baseunit(c['area'])    # m2 - conductor cross-section
        self.cold_fraction    = c['cold_fraction']     # unitless - length fraction at cryo temperature
        self.resistivity_warm = c['resistivity_warm']  # ohm-m - Cu near ambient
        self.resistivity_cold = c['resistivity_cold']  # ohm-m - Cu at LOX temperature (cryo credit)
        self.inductance       = c['inductance']        # H - series line inductance

