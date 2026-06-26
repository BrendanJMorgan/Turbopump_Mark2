import os
import numpy as np
import yaml

# ---------------------------------------------------------------------
# Load primitive inputs from config.yaml.
#
# Access is by direct key indexing throughout, so a missing or renamed
# parameter raises KeyError at import instead of silently defaulting.
# Derived quantities are computed below from the primitives so they can
# never drift out of sync with their parents.
# ---------------------------------------------------------------------
_CONFIG_PATH = os.path.join(os.path.dirname(__file__), 'config.yaml')

with open(_CONFIG_PATH, 'r') as _f:
    config = yaml.safe_load(_f)


class engine:  # For parameters that don't fit neatly within other categories
    _c = config['engine']
    g               = _c['g']                # m/s2
    thrust          = _c['thrust']           # N - Thrust INCLUDE GG EXHAUST THRUST
    p_amb           = _c['p_amb']            # Pa - ambient pressure at 2100 feet elevation
    T_amb           = _c['T_amb']            # K - Ambient Temperature
    fuel            = _c['fuel']             # 'JetA' 'C2H5OH(L)'
    fuelrp          = _c['fuelrp']           # rocketprop only has RP-1 not Jet-A
    oxidizer        = _c['oxidizer']         # 'O2(L)'
    proof           = _c['proof']            # Proportion of fuel (remaining part is water)
    shaft_speed_cap = _c['shaft_speed_cap']  # rad/s - max shaft speed if not suction limited
    common_shaft    = _c['common_shaft']     # Whether the pumps share a common shaft
    del _c


class tca:  # Thrust Chamber Assembly
    _c = config['tca']

    # Chamber/Nozzle Geometry
    dx             = _c['dx']              # m - position step
    converge_angle = _c['converge_angle']  # rad - half-cone convergence angle of chamber end
    diverge_angle  = _c['diverge_angle']   # rad - half-cone divergence angle of nozzle
    l_star         = _c['l_star']          # m - characteristic chamber length
    rc_throat      = _c['rc_throat']       # m - radius of curvature around the throat
    d2_chamber     = _c['d2_chamber']      # m
    wall_thickness = _c['wall_thickness']  # m - SLANT VS VERTICAL THICKNESS
    d1_chamber     = d2_chamber - 2 * wall_thickness  # m  (derived)
    r1_chamber     = d1_chamber / 2                   # m  (derived)

    pc                = _c['pc']                 # Pa - Chamber (Stagnation) Pressure
    mixture_ratio     = _c['mixture_ratio']      # Oxidizer/Fuel Ratio (by mass)
    c_star_efficiency = _c['c_star_efficiency']  # Characteristic Vel Efficiency, experimental
    c_tau_efficiency  = _c['c_tau_efficiency']   # Thrust Coefficient Efficiency, experimental

    p_injector_manifold = pc * 1.25  # Pa - PLACEHOLDER UNTIL A PROPER INJECTOR MODULE IS WRITTEN (derived)
    compute_thermals    = _c['compute_thermals']  # Whether to run regen thermal balance

    class combustion_gas:
        p = None
        T = None

    class film_coolant:
        _c = config['tca']['film_coolant']
        fraction             = np.array(_c['fraction'])  # unitless - fuel mass fraction per film-cooling row
        injection_x          = _c['injection_x']         # m - film cooling orifice positions
        v_injection          = _c['v_injection']         # m/s - initial gas velocity for film model
        injection_efficiency = _c['injection_efficiency']
        del _c

        class liquid:
            None

        class gas:
            None

    class regenerative_coolant:
        flow_direction = config['tca']['regenerative_coolant']['flow_direction']  # 1 = forward, -1 = counter

    class regenerative_jacket:
        _c = config['tca']['regenerative_jacket']
        n_pipe1      = _c['n_pipe1']       # number of channels along barrel
        n_pipe2      = _c['n_pipe2']       # number of channels near throat
        n_pipe3      = _c['n_pipe3']       # number of channels along lower nozzle section
        gap_pipe     = _c['gap_pipe']      # m - gap between channels (fin thickness)
        h_pipe       = _c['h_pipe']        # m - coolant channel height
        merge_radius = _c['merge_radius']  # m - contour radius to transition to n_pipe2
        cond         = _c['cond']          # W/m-K - thermal conductivity of copper
        del _c

    del _c


class gg:  # Gas Generator
    _c = config['gg']
    stiffness           = _c['stiffness']            # Pa/Pa - guess
    fraction_guess      = _c['fraction_guess']       # unitless - fraction of total mass flow to gg
    mixture_ratio       = _c['mixture_ratio']        # unitless - Oxidizer/Fuel Ratio
    c_star_efficiency   = _c['c_star_efficiency']    # characteristic velocity efficiency, experimental
    c_tau_efficiency    = _c['c_tau_efficiency']     # unitless - Thrust Coefficient Efficiency Factor
    orifice_number_fuel = _c['orifice_number_fuel']  # 4 impinging, 8 showerhead film cooling
    r_chamber           = _c['r_chamber']            # m - radius of combustion chamber
    del _c


class turbine:
    _c = config['turbine']
    r_pitchline        = _c['r_pitchline']        # m
    gear_ratio         = _c['gear_ratio']         # unitless - higher = faster/smaller turbine; 1 = common shaft
    d_throat_nozzle    = _c['d_throat_nozzle']    # m - diameter of each nozzle off the manifold
    diverge_angle_gg   = _c['diverge_angle_gg']   # rad - half-cone divergence angle of nozzle plate
    blade_width        = _c['blade_width']        # m
    admission_fraction = _c['admission_fraction']  # unitless - fraction of nozzle-plate circle with nozzles
    blade_length       = _c['blade_length']       # m
    tip_clearance      = _c['tip_clearance']      # m
    radius_leading     = _c['radius_leading']     # m - blade leading edge fillet radius. NEEDS MORE RESEARCH
    del _c


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
            self.plumbing_loss = c['plumbing_loss']  # Pa - fuel side only


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
        self.sweepback             = c['sweepback']              # rad - blade leading-edge sweepback


class impeller:
    def __init__(self, stage: int, fluid: str):
        self.fluid = fluid
        self.blade = blade(fluid=fluid)
        c = config['pump'][fluid]['impeller']
        self.stage             = stage
        self.hub_in_angle      = c['hub_in_angle']
        self.hub_out_angle     = c['hub_out_angle']
        self.shroud_in_angle   = c['shroud_in_angle']
        self.shroud_out_angle  = c['shroud_out_angle']
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
        self.wrap_angle = c['wrap_angle']  # rad - single volute is 2*pi, double volute is pi, etc.


class tank:
    def __init__(self, fluid: str):
        self.fluid = fluid
        c = config['pump'][fluid]['tank']
        self.p = c['p']  # Pa
        self.T = c['T']  # K
