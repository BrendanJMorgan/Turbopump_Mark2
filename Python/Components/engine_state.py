import os
import numpy as np
import yaml
from pint import UnitRegistry
ureg = UnitRegistry()

# Specific speed units. Base form is the dimensionless w*sqrt(Q)/(g*H)^0.75 (rad/s, m3/s, m, standard g),
# so these are just scale factors on dimensionless. Same units serve suction specific speed (H -> NPSH).
_ss_US = (ureg('rpm').to('rad/s') * ureg('gal/min').to('m**3/s')**0.5 / (ureg('g0') * ureg('ft')).to_base_units()**0.75).to_base_units().magnitude
_ss_EU = (ureg('rpm').to('rad/s') / ureg('g0').to_base_units()**0.75).to_base_units().magnitude
ureg.define(f"specific_speed_US = {_ss_US!r} = ns_US")  # rpm, gpm, ft  (~1/2733)
ureg.define(f"specific_speed_EU = {_ss_EU!r} = ns_EU")     # rpm, m3/s, m  (Gulich n_q, n_ss)

def baseunit(value):
    """Convert a config {value, unit} entry to its base-SI magnitude. Bare numbers are rejected on purpose."""
    if isinstance(value, dict) and 'value' in value and 'unit' in value:
        magnitude = value['value']
        if isinstance(magnitude, list):
            magnitude = np.array(magnitude, dtype=float)
        return ureg.Quantity(magnitude, value['unit']).to_base_units().magnitude
    else:
        raise ValueError(f"Config value {value!r} has no unit - write it as {{value: ..., unit: ...}}")

CONFIG_PATH = os.path.join(os.path.dirname(__file__), 'config.yaml') # Load primitive inputs from config.yaml.

with open(CONFIG_PATH, 'r') as _f:
    config = yaml.safe_load(_f)

class engine:  # For parameters that don't fit neatly within other categories
    c = config['engine']
    g               = baseunit(c['g'])                # m/s2
    thrust          = baseunit(c['thrust'])           # N - Thrust INCLUDE GG EXHAUST THRUST
    p_amb           = baseunit(c['p_amb'])            # Pa - ambient pressure at 2100 feet elevation
    T_amb           = baseunit(c['T_amb'])            # K - Ambient Temperature
    fuel            = c['fuel']             # 'JetA' 'C2H5OH(L)'
    fuelrp          = c['fuelrp']           # rocketprop only has RP-1 not Jet-A
    oxidizer        = c['oxidizer']         # 'O2(L)'
    proof           = baseunit(c['proof'])            # Proportion of fuel (remaining part is water)
    shaft_speed_cap = baseunit(c['shaft_speed_cap'])  # rad/s - max shaft speed if not suction limited
    common_shaft    = c['common_shaft']     # Whether the pumps share a common shaft
    del c

class tca:  # Thrust Chamber Assembly
    c = config['tca']

    # Chamber/Nozzle Geometry
    dx             = baseunit(c['dx'])              # m - position step
    converge_angle = baseunit(c['converge_angle'])  # rad - half-cone convergence angle of chamber end
    diverge_angle  = baseunit(c['diverge_angle'])   # rad - half-cone divergence angle of nozzle
    l_star         = baseunit(c['l_star'])          # m - characteristic chamber length
    rc_throat      = baseunit(c['rc_throat'])       # m - radius of curvature around the throat
    d2_chamber     = baseunit(c['d2_chamber'])  # m
    wall_thickness = baseunit(c['wall_thickness'])  # m - SLANT VS VERTICAL THICKNESS
    d1_chamber     = d2_chamber - 2 * wall_thickness  # m  (derived)
    r1_chamber     = d1_chamber / 2                   # m  (derived)

    pc                = baseunit(c['pc'])                 # Pa - Chamber (Stagnation) Pressure
    mixture_ratio     = baseunit(c['mixture_ratio'])      # Oxidizer/Fuel Ratio (by mass)
    c_star_efficiency = baseunit(c['c_star_efficiency'])  # Characteristic Vel Efficiency, experimental
    c_tau_efficiency  = baseunit(c['c_tau_efficiency'])   # Thrust Coefficient Efficiency, experimental

    p_injector_manifold = pc * 1.25  # Pa - PLACEHOLDER UNTIL A PROPER INJECTOR MODULE IS WRITTEN (derived)
    compute_thermals    = c['compute_thermals']  # Whether to run regen thermal balance

    class combustion_gas:
        p = None
        T = None

    class film_coolant:
        c = config['tca']['film_coolant']
        fraction             = baseunit(c['fraction'])  # unitless - fuel mass fraction per film-cooling row
        injection_x          = baseunit(c['injection_x'])         # m - film cooling orifice positions
        v_injection          = baseunit(c['v_injection'])         # m/s - initial gas velocity for film model
        injection_efficiency = baseunit(c['injection_efficiency'])
        del c

        class liquid:
            None

        class gas:
            None

    class regenerative_coolant:
        c = config['tca']['regenerative_coolant']
        flow_direction   = baseunit(c['flow_direction'])    # 1 = forward, -1 = counter
        del c

    class regenerative_jacket:
        c = config['tca']['regenerative_jacket']
        n_pipe1      = baseunit(c['n_pipe1'])       # number of channels along barrel
        n_pipe2      = baseunit(c['n_pipe2'])       # number of channels near throat
        n_pipe3      = baseunit(c['n_pipe3'])       # number of channels along lower nozzle section
        gap_pipe     = baseunit(c['gap_pipe'])  # m - gap between channels (fin thickness)
        h_pipe       = baseunit(c['h_pipe'])    # m - coolant channel height
        merge_radius = baseunit(c['merge_radius'])  # m - contour radius to transition to n_pipe2
        cond         = baseunit(c['cond'])          # W/m-K - thermal conductivity of copper
        del c

    del c

class gg:  # Gas Generator
    c = config['gg']
    stiffness           = baseunit(c['stiffness'])            # Pa/Pa - guess
    fraction_guess      = baseunit(c['fraction_guess'])       # unitless - fraction of total mass flow to gg
    mixture_ratio       = baseunit(c['mixture_ratio'])        # unitless - Oxidizer/Fuel Ratio
    c_star_efficiency   = baseunit(c['c_star_efficiency'])    # characteristic velocity efficiency, experimental
    c_tau_efficiency    = baseunit(c['c_tau_efficiency'])     # unitless - Thrust Coefficient Efficiency Factor
    orifice_number_fuel = baseunit(c['orifice_number_fuel'])  # 4 impinging, 8 showerhead film cooling
    r_chamber           = baseunit(c['r_chamber'])            # m - radius of combustion chamber
    del c

class turbine:
    c = config['turbine']
    r_pitchline_max    = baseunit(c['r_pitchline_max'])  # m - manufacturing limit
    velocity_ratio     = baseunit(c['velocity_ratio'])   # unitless - target pitchline speed / spouting velocity
    gear_ratio         = baseunit(c['gear_ratio'])         # unitless - higher = faster/smaller turbine; 1 = common shaft
    d_throat_nozzle    = baseunit(c['d_throat_nozzle'])  # m - diameter of each nozzle off the manifold
    diverge_angle_gg   = baseunit(c['diverge_angle_gg']) # rad - half-cone divergence angle of nozzle plate
    blade_width        = baseunit(c['blade_width'])      # m
    admission_fraction = baseunit(c['admission_fraction'])  # unitless - fraction of nozzle-plate circle with nozzles
    blade_length       = baseunit(c['blade_length'])     # m
    tip_clearance      = baseunit(c['tip_clearance'])      # m
    radius_leading     = baseunit(c['radius_leading'])     # m - blade leading edge fillet radius. NEEDS MORE RESEARCH
    del c

class pump:
    def __init__(self, fluid: str):
        self.fluid = fluid
        c = config['pump'][fluid]  # KeyError here if fluid is not configured

        self.surface_roughness = baseunit(c['surface_roughness'])
        self.shaft_speed = None
        self.clocking = baseunit(c['clocking'])    # 1 for counterclockwise and -1 for clockwise (looking down at inlet)
        self.n_stages = baseunit(c['n_stages'])
        self.inducer = pump.inducer(fluid=self.fluid)
        self.impeller = [pump.impeller(stage=i, fluid=self.fluid) for i in range(self.n_stages)]
        if 'volute' in c:
            self.volute = pump.volute(fluid=self.fluid)
        if 'vaned_diffuser' in c:
            self.vaned_diffuser = pump.vaned_diffuser(fluid=self.fluid)
        self.tank = tank(fluid=self.fluid)
        if 'plumbing_loss' in c:
            self.plumbing_loss = baseunit(c['plumbing_loss'])  # Pa - discharge pressure allowance for line losses downstream of the pump
        if 'mainline' in c:
            self.mainline = pump.mainline(fluid=self.fluid)

    class inducer:
        def __init__(self, fluid: str):
            self.fluid = fluid
            c = config['pump'][fluid]['inducer']
            self.flow_coeff            = baseunit(c['flow_coeff'])             # unitless
            self.blade_number          = baseunit(c['blade_number'])           # unitless
            self.clearance_radial      = baseunit(c['clearance_radial'])       # m
            self.hub_tip_ratio         = baseunit(c['hub_tip_ratio'])          # unitless - inlet hub/tip radius ratio
            self.hub_tip_ratio_out     = baseunit(c['hub_tip_ratio_out'])      # unitless - outlet hub/tip radius ratio
            self.suction_margin_factor = baseunit(c['suction_margin_factor'])  # unitless - margin of extra NPSH
            self.flow_margin_factor    = baseunit(c['flow_margin_factor'])     # unitless - margin on flow rate
            self.blade_thickness       = baseunit(c['blade_thickness'])        # m - thickness of inducer blades
            self.tip_clearance         = baseunit(c['tip_clearance'])          # m
            self.length_pitch_ratio    = baseunit(c['length_pitch_ratio'])     # max empirical value per Gulich
            self.sweepback             = baseunit(c['sweepback'])    # rad - blade leading-edge sweepback
            self.incidence             = baseunit(c['incidence'])    # rad - blade angle minus flow angle at the leading edge
            del c

    class impeller:
        def __init__(self, stage: int, fluid: str):
            self.fluid = fluid
            self.blade = pump.impeller.blade(fluid=fluid)
            c = config['pump'][fluid]['impeller']
            self.stage             = stage
            self.hub_in_angle      = baseunit(c['hub_in_angle'])      # rad
            self.hub_out_angle     = baseunit(c['hub_out_angle'])     # rad
            self.shroud_in_angle   = baseunit(c['shroud_in_angle'])   # rad
            self.shroud_out_angle  = baseunit(c['shroud_out_angle'])  # rad
            self.blade_thickness_ratio = baseunit(c['blade_thickness_ratio'])  # unitless - blade thickness / outlet diameter
            self.suction_specific_speed = baseunit(c['suction_specific_speed'])  # unitless - w*sqrt(Q)/(g*NPSH3)^0.75
            self.suction_margin_factor  = baseunit(c['suction_margin_factor'])   # unitless - NPSH_available / NPSH3
            self.lambda_c          = baseunit(c['lambda_c'])          # unitless - Gulich 6.10 inlet acceleration/loss coefficient
            self.lambda_w          = baseunit(c['lambda_w'])          # unitless - Gulich 6.10 blade leading-edge coefficient (NPSH3)
            self.surface_roughness = baseunit(c['surface_roughness']) # m
            del c

        class blade:
            def __init__(self, fluid: str):
                self.fluid = fluid

    class volute:
        def __init__(self, fluid: str):
            self.fluid = fluid
            c = config['pump'][fluid]['volute']
            self.wrap_angle = baseunit(c['wrap_angle'])  # rad - single volute is 2*pi, double volute is pi, etc.
            if 'discharge_diameter' in c:
                self.discharge_diameter = baseunit(c['discharge_diameter'])  # m - discharge line diameter
            self.discharge_cone_angle = baseunit(c['discharge_cone_angle'])  # rad - total included angle of the discharge nozzle
            self.discharge_max_length = baseunit(c['discharge_max_length'])  # m - packaging limit on discharge nozzle length
            self.inlet_width_ratio  = baseunit(c['inlet_width_ratio'])   # unitless - volute inlet width / impeller outlet width
            del c

    class vaned_diffuser:
        def __init__(self, fluid: str):
            self.fluid = fluid
            c = config['pump'][fluid]['vaned_diffuser']
            self.inlet_width_ratio   = baseunit(c['inlet_width_ratio'])    # unitless - diffuser inlet width / impeller outlet width
            self.incidence           = baseunit(c['incidence'])            # rad - vane angle minus flow angle at the leading edge
            self.throat_width_factor = baseunit(c['throat_width_factor'])  # unitless - f_a3, Gulich 7.7.10
            self.outlet_factor       = baseunit(c['outlet_factor'])        # unitless - Gulich 7.9.1
            del c

    class mainline:
        def __init__(self, fluid: str):
            self.fluid = fluid
            c = config['pump'][fluid]['mainline']
            self.tube_od_options    = baseunit(c['tube_od_options'])             # m - standard tube outer diameters
            self.tube_wall_options  = baseunit(c['tube_wall_options'])           # m - matching wall thicknesses
            self.tube_material_options     = list(c['tube_material_options'])          # material names
            self.tube_density       = np.array([baseunit(m['density']) for m in c['tube_material_options'].values()])           # kg/m3
            self.allowable_stress   = np.array([baseunit(m['allowable_stress']) for m in c['tube_material_options'].values()])  # Pa
            self.dp_mass_equivalent = baseunit(c['dp_mass_equivalent'])  # kg/Pa - cost weighting of pressure drop against mass
            self.length             = baseunit(c['length'])            # m
            self.roughness          = baseunit(c['roughness'])         # m
            self.minor_losses       = {name: baseunit(k['K'])*baseunit(k['count']) for name, k in c['minor_losses'].items()}  # unitless - K per item type
            self.minor_loss_coeff   = sum(self.minor_losses.values())  # unitless - sum of K
            del c

class tank:
    def __init__(self, fluid: str):
        self.fluid = fluid
        c = config['pump'][fluid]['tank']
        self.p = baseunit(c['p'])  # Pa

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
            self.reactance_perunit   = baseunit(c['reactance_perunit'])          # unitless - per-unit synchronous reactance
            self.resistance_phase    = baseunit(c['resistance_phase'])      # ohm - per-phase winding resistance
            self.core_loss_fraction  = baseunit(c['core_loss_fraction'])    # unitless - crude core loss / throughput
            self.remanence           = baseunit(c['remanence'])             # T - magnet B_r (SmCo)
            self.relative_recoil_permeability = baseunit(c['relative_recoil_permeability'])   # unitless - mu_r
            self.magnet_length       = baseunit(c['magnet_length'])         # m - axial magnet thickness
            self.airgap              = baseunit(c['airgap'])                # m - mechanical airgap
            self.r_inner             = baseunit(c['r_inner'])               # m - annulus ID (centre bearing bore)
            self.shear_allowable     = baseunit(c['shear_allowable'])             # Pa - airgap shear ceiling

    class alternator:
        def __init__(self, fluid: str, poles: int, ang_frequency: float, r_max: float):
            self.fluid = fluid
            self.poles = poles
            self.ang_frequency = ang_frequency
            self.r_max = r_max
            c = config['machine']['alternator']
            self.reactance_perunit   = baseunit(c['reactance_perunit'])          # unitless - per-unit synchronous reactance
            self.resistance_phase    = baseunit(c['resistance_phase'])      # ohm - per-phase winding resistance
            self.core_loss_fraction  = baseunit(c['core_loss_fraction'])    # unitless
            self.remanence           = baseunit(c['remanence'])             # T - magnet B_r (SmCo)
            self.relative_recoil_permeability = baseunit(c['relative_recoil_permeability'])   # unitless - mu_r
            self.magnet_length       = baseunit(c['magnet_length'])         # m
            self.airgap              = baseunit(c['airgap'])                # m
            self.r_inner             = baseunit(c['r_inner'])               # m
            self.shear_allowable           = baseunit(c['shear_allowable'])             # Pa

    class cable:
        def __init__(self, frequency: float):
            self.frequency = frequency  # Hz - electrical frequency of the link
            c = config['machine']['cable']
            self.length           = baseunit(c['length'])            # m - conductor run length
            self.area             = baseunit(c['area'])    # m2 - conductor cross-section
            self.cold_fraction    = baseunit(c['cold_fraction'])     # unitless - length fraction at cryo temperature
            self.resistivity_warm = baseunit(c['resistivity_warm'])  # ohm-m - Cu near ambient
            self.resistivity_cold = baseunit(c['resistivity_cold'])  # ohm-m - Cu at LOX temperature (cryo credit)
            self.inductance       = baseunit(c['inductance'])        # H - series line inductance

