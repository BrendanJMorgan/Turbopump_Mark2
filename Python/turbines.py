import numpy as np

def turbines(p) -> None:

    # load(expansion_loss_data.mat)  # in your project: supply this externally or convert to arrays on `p`
    blade_width_interp = np.array([0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0])  # in
    incidence_loss_interp = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # TODO: populate with your curve
    mach_loss_data = None  # TODO: your table
    load('impulse_turbine_efficiency.mat')  # TODO: replace with your accessor
    load('pitch_chord_real.mat')            # TODO: replace with your accessor

    p.turbine_shaft_speed = abs(p.pump_shaft_speed) * p.gear_ratio  # rad/s

    # XX Nozzle Plate
    mach_spouting = np.sqrt(2 / (p.gamma_BB - 1) * ((p.p_amb / p.p_BB)**((p.gamma_BB - 1) / p.gamma_BB) - 1))  # unitless
    T_spouting = p.T_gB * (p.p_amb / p.p_BB)**((p.gamma_BB - 1) / p.gamma_BB)  # K - static temperature for perfectly expanded flow
    M_BB = 0.068 * p.BB_output_oz1 * (33.91 / p.M0_k)  # molecular weight (fully combusted)   # TODO: confirm units/fit
    R_BB = 8.3145 / (M_BB / 1000.0)  # J/(kg·K) - specific gas constant
    v_spouting = np.sqrt(2 * p.gamma_BB / (p.gamma_BB - 1) * R_BB * p.T_gB * (1 - (p.p_amb / p.p_BB)**((p.gamma_BB - 1) / p.gamma_BB)))  # m/s
    v_exit_nozzle = v_spouting * p.gamma_v_ratio  # m/s - speed on nozzle exit plane

    A_throat_np = p.mdot_BB_gas / (p.rho_gB * v_exit_nozzle)  # kg/s ÷ (kg/m^3·m/s) → m^2 = cumulative throat area on the nozzle plate
    r_pitchline_nozzle = p.r_pitchline_rotor  # m  # TODO: set from geometry if different

    # XX Velocities
    # one stage impulse rotor - does not work for multistaged and/or velocity driven turbines
    v_pitchline_rotor = r_pitchline_nozzle * p.turbine_shaft_speed  # m/s - tangential rotor speed at pitchline

    v_ratio = 0.35  # unitless - “A single-row impulse stage delivers best performance at velocity ratios between 0.30 and 0.40”
    v_axial_noz = v_spouting * v_ratio            # m/s
    v_tang_noz = np.sqrt(v_spouting**2 - v_axial_noz**2)  # [tangential, axial]
    v_rel_inlet = np.hypot(v_tang_noz - v_pitchline_rotor, v_axial_noz)  # m/s
    blade_relative = np.arctan2(v_axial_noz, v_tang_noz - v_pitchline_rotor)  # rad - this is measured FROM THE VERTICAL

    pi_entrance_vector_angle = np.arctan2(v_axial_noz, v_tang_noz)  # rad
    exit_vector_angle = -blade_relative                              # rad - assume symmetric blades until proven otherwise

    blade_chord_angle = np.pi/2 - blade_relative  # rad - measure FROM THE HORIZONTAL
    pitchline_nozzle_angle = np.pi/2 - pi_entrance_vector_angle  # rad - measure FROM THE HORIZONTAL

    M_nozzle_exit = v_spouting / np.sqrt(p.gamma_BB * R_BB * T_spouting)  # unitless - mach at nozzle plate exit
    v_spouting_R = np.hypot(v_tang_noz - v_pitchline_rotor, v_axial_noz)  # unitless - spouting in corotating frame of blades

    # XX Blade Spacing
    A_min_throat = A_throat_np / p.BB  # m^2
    A_blade_gap = A_min_throat / p.admission_fraction  # m^2 = minimum throughput area needed to prevent gas backup
    # TODO: convert your exact Zweifel correlations here:
    incidence_angle_blade = (blade_chord_angle * 180 / np.pi)  # deg
    Zweifel = 0.85  # unitless (typical)
    exit_vector_angle_deg = exit_vector_angle * 180 / np.pi
    p_blade_pitch = pitch_chord(incidence_angle_blade, exit_vector_angle_deg, Zweifel)  # in; supply from your lookup
    blade_pitch_rotor = p_blade_pitch / 39.37  # m

    blade_chord_rotor = blade_pitch_rotor / np.cos(blade_chord_angle)  # ~ in an impulse rotor these are equal
    blade_opening_rotor = blade_pitch_rotor * np.sin(blade_chord_angle)  # m - circumferential distance between neighboring blades
    p.blade_count_rotor = int(np.round(2 * np.pi * r_pitchline_nozzle / blade_pitch_rotor))  # unitless - number of blades

    # XX Energies
    # Adjust turbine gas-path expansion and kinetic energy, using the respective two-dimensional
    # loss coefficients for expansion energy, kinetic energy, Mach number, and incidence angle as
    # determined from plots like those in figures 16 through 18.
    expansion_loss_coeff = expansion_loss_coeff_fn(v_ratio, blade_width_interp, incidence_angle_blade, exit_vector_angle_deg, blade_width_rotor=blade_pitch_rotor/0.0254)  # unitless - NASA turbines 1976; fig. 16
    kinetic_loss_coeff = kinetic_loss_coeff_fn(M_nozzle_exit)  # unitless; eqn 3
    incidence_loss_coeff = incidence_loss_coeff_fn(incidence_angle_blade)  # unitless - NASA turbines 1976; fig. 17
    mach_loss_coeff = mach_loss_coeff_fn(M_nozzle_exit)  # unitless - NASA turbines 1976; fig. 18

    nozzle_area = 2 * np.pi * admission_fraction(r_pitchline_nozzle) * pitchline_nozzle_angle  # TODO: your expression
    turbine_efficiency_partial = ((1 - p.wheel_ratio) * (1 - blade_pitch_rotor / (p.nozzle_cu))) / ((1 - p.wheel_ratio) * p.turbine_efficiency_full * 1.53e-9 * p.gamma_v_ratio * p.dens_gas_spouted * p.turbine_shaft_speed * (1 - p.admission_fraction) / p.admission_fraction)  # unitless - equation 14  # TODO: verify
    turbine_efficiency = turbine_efficiency_partial * clearance_loss_coeff(p.clearances_rotor)  # unitless

    specific_work_expansion = v_spouting**2 / 2 / 9.80665  # J/kg → m (placeholder for eqn 1b); replace with your eqn set
    kinetic_specific_work = v_rel_inlet**2 / 2 / 9.80665   # J/kg → m (placeholder for eqn 5b)
    specific_work_rotor = expansion_specific_work(specific_work_expansion, expansion_loss_coeff) + kinetic_specific_work * kinetic_loss_coeff * mach_loss_coeff * incidence_loss_coeff  # J/kg

    # mass flow through the GG / manifold / rotor
    p.mdot_gas_rotor = p.power_turbine / (turbine_shaft_specific_work(v_spouting, v_pitchline_rotor) * turbine_efficiency)  # kg/s

    # XX Blade Depth
    if blade_length(blade_chord_rotor) < 1.5 * p.blade_count_rotor * blade_opening_rotor:
        print("Rotor blades are too short and choking the flow")

    p.r_tip_rotor = r_pitchline_nozzle + blade_length_nozzle(p) / 2  # m
    p.r_base_rotor = r_pitchline_nozzle - blade_length_nozzle(p) / 2  # m

    # XX Flow Adjustments
    mdot_BB = p.mdot_BB  # kg/s
    mdot_fuel_BB = mdot_BB * p.fBB  # kg/s - Fuel Mass Flow Rate
    mdot_ox_BB = mdot_BB * (1 - p.fBB)  # kg/s - Oxidizer Mass Flow Rate
    p.mdot_fuel_total = p.mdot_fuel_cc + mdot_fuel_BB  # kg/s
    p.mdot_ox_total = p.mdot_ox_cc + mdot_ox_BB       # kg/s
    p.mdot_total = p.mdot_ox_total + p.mdot_fuel_total  # kg/s

    # XX Rotor Stress
    p.sigma_hoop = 2 * (p.r_tip_rotor - p.r_base_rotor) / 0.0254**2 * (p.turbine_shaft_speed * 30 / np.pi)**2  # Pa  # TODO: confirm units
    rotor_density = 7860.0  # kg/m^3 - 316 stainless  # https://bssa.org.uk/bssa_articles/elevated-temperature-physical-properties-of-stainless-steels/
    allow_stress = 0.66 * 515e6  # Pa
    rotor_strength = 0.95 * allow_stress  # Pa
    q_tip_rotor = 0.25
    p.pressure_rotor = 0.5 * rotor_density * (p.turbine_shaft_speed**2) * (q_tip_rotor**2)  # Pa
    p.sigma_rotor = p.pressure_rotor / (rotor_strength / 6.89476e3)  # ksi

    # XX Other
    p.v_spouted_gg = v_spouting / 2  # m/s
    p.spin_bottle_gg = 80.42 * (p.v_spouted_gg / 33.91)  # “this represents a cold keg of nitrogen in a K bottle. in reality flow will stop when pressure drops too low.”

