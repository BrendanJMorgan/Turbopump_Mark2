import numpy as np

from engine_state import engine, pump

def inducer(p) -> None:
    # NASA SP-8052 "LIQUID ROCKET ENGINE TURBOPUMP INDUCER" 1971

    # Net Positive Suction Head (NPSH) Requirements
    k1 = 1.2
    k2 = 0.2334 + (p.pump_shaft_speed * p.r_inlet_impeller / 128.3) ** 4
    p.cavitation_limit = (p.pump_shaft_speed * p.r_inlet_impeller) ** 2 / (2 * g)  # unitless
    p.NPSH_SE = ((k1 * k2) * p.eye_flow_coeff**2 + k2) * p.cavitation_limit  # m - NPSH, shockless entry. Pump handbook pg. 2.88
    p.NPSH_required = p.NPSH_SE  # m - assuming operation at best efficiency point (BEP) and shockless

    p.NPSH_available = (p.p_amb - p.vapor_pressure_pump) / (p.density_pump * engine.g)  # m - Net positive suction head. Assuming zero plumbing water column and inlet velocity
    p.NPSH_inducer = (p.NPSH_required - p.NPSH_available) * p.NPSH_margin  # m

    # Tip Diameter
    p.vdot_corrected = p.vdot_pump / (1 - p.hub_tip_ratio_inducer**2)  # m^3/s - volumetric flow rate through open portion of inducer
    p.ss_speed = p.pump_shaft_speed * np.sqrt(p.vdot_pump) / ((g * p.NPSH_inducer) ** 0.75)  # unitless - suction specific speed
    p.flow_coeff_inducer = (
        1.3077 * np.sqrt(1 + p.hub_tip_ratio_inducer**2) / p.ss_speed
        * (1 + 0.5 * np.sqrt(1 + 6 * (1.3077 * np.sqrt(1 + p.hub_tip_ratio_inducer**2) / p.ss_speed) ** 2))
    )  # unitless - phi
    p.r_tip_inducer = 0.5 * 0.3048 * 0.37843 * (
        (p.vdot_corrected * 15580.323) / (p.pump_shaft_speed * 30 / np.pi * p.flow_coeff_inducer)
    ) ** (1 / 3)  # m
    p.r_hub_inducer = p.r_tip_inducer * p.hub_tip_ratio_inducer  # m
    p.r_pitchline_inducer = 0.5 * (1 + p.hub_tip_ratio_inducer) * p.r_tip_inducer  # m

    # Blade Construction
    v_in = np.array([0.0, 0.0])
    v_in[0] = p.pump_shaft_speed * p.r_tip_inducer  # m/s - tangential velocity at inducer pitchline tip
    v_in[1] = p.vdot_corrected / (np.pi * p.r_tip_inducer**2)  # m/s - axial velocity at inducer pitchline
    p.v_inlet_inducer = v_in

    p.flow_blade_ratio = 0.35  # unitless
    p.flow_incidence_angle = np.arctan(v_in[1] / v_in[0])  # rad - alpha
    p.blade_angle_inducer = p.flow_incidence_angle / (1 - p.flow_blade_ratio)  # rad - beta

    p.blade_lead_inducer = 2 * np.pi * p.r_tip_inducer * np.tan(p.blade_angle_inducer)  # m - blade lead - how much the blade advances per turn
    p.solidity_inducer = 2.5
    p.h_min_inducer = (
        p.blade_lead_inducer * p.solidity_inducer / p.blade_number_inducer * np.sin(p.blade_angle_inducer)
    )  # m - height of inducer

    # Axial Clearance
    # % Between inducer and impeller
    p.clearance_axial_inducer = p.blade_lead_inducer / p.blade_number_inducer  # m

    # Head
    p.v_tip_inducer = p.pump_shaft_speed * p.r_tip_inducer  # m/s
    p.head_coeff_inducer = p.NPSH_inducer * engine.g / (p.v_tip_inducer**2)  # unitless - 0.15 is the cutoff between low and high head - camber needed if >0.075

    # Hub Construction
    p.r_hub_base = p.r_hub  # m
    p.r_hub_top = 2 * p.hub_tip_ratio_inducer * p.r_tip_inducer - p.r_hub_base  # m

    # Clearance Losses - not relevant with a proper shroud
    p.r_cavity_inducer = p.r_tip_inducer + p.clearance_radial_inducer  # m
    p.ss_speed_cleared = p.ss_speed * (1 - 0.575 * np.sqrt((p.r_cavity_inducer - p.r_tip_inducer) / (p.r_tip_inducer - p.r_hub_inducer)))  # unitless
    p.head_coeff_inducer_cleared = p.head_coeff_inducer * (1 - 1.0 * np.sqrt((p.r_cavity_inducer - p.r_tip_inducer) / (p.r_tip_inducer - p.r_hub_inducer)))  # unitless

    # Other
    p.inducer_torque = p.mdot_pump * p.pump_shaft_speed * p.r_tip_inducer**2  # N*m - this is a conservative maximum
