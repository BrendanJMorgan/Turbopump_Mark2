from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from engine_state import pump, turbine, link, config, baseunit

vacuum_permeability = 4.0e-7 * np.pi  # H/m

def size_machine(m, winding_factor, temperature):

    ### Requirement
    m.torque = m.power / m.ang_frequency    # N-m
    m.airgap_shear_stress = 3*m.torque / (2*np.pi * (m.r_outer**3 - m.r_inner**3)) # Pa - sigma - Hanselman 4.43 - inegrated over annulus; assumes constant flux density

    m.shear_margin = m.shear_allowable / m.airgap_shear_stress
    if m.shear_margin < 1.0:
        raise ValueError(f"machine shear stress limit (factor {m.shear_margin:.2f}): increase annulus area or reduce torque")
    
    ### Magnetic Circuit
    m.airgap_area = np.pi * (m.r_outer**2 - m.r_inner**2)   # m2 - active annulus (one side)
    m.magnet_area = m.airgap_area # PROBABLY BAD ASSUMPTION

    m.pm_reluctance = m.magnet_length / (m.magnet_permeability * m.magnet_area)  # 1/H - magnet reluctance - Hanselman 4.2
    m.airgap_reluctance = m.airgap_length / (m.airgap_permeability * m.airgap_area)  # 1/H - airgap reluctance - Hanselman 4.2

    reluctance_factor = 1 # unitless - Hanselman fig 4-3d, K_r. Empirical; adds reluctance of steel; typically 1.0 - 1.2
    leakage_factor = 1 # unitless - Hanselman, K_l. Empirical; leakage between PM's; typically 0.9 - 1.0
    m.airgap_flux_density = leakage_factor * m.pm_remanence / (1 + reluctance_factor * m.airgap_reluctance / m.pm_reluctance) # T - B_g - Hanselman 4.2

    ### Pole sizing
    m.torque_constant = m.winding_turns*m.airgap_flux_density*(m.r_outer**2-m.r_inner**2) # N*m/A - K_t - Hanselman 4.11
    m.current_rms = m.torque / m.torque_constant # A - i - Hanselman 4.11
    m.back_emf_amplitude = m.pole_count * m.torque_constant * m.ang_frequency # V - E_b - Hanselman 4.10

    ### Power efficiency
    fudge_factor = 0.8 # winding turns don't have perfectly sharp corners; wire has finite thickness; etc.
    m.mean_turn_length = fudge_factor * (2*np.pi*(m.r_outer+m.r_inner) + 2*(m.r_outer-m.r_inner)) / m.pole_count # m
    winding_length_total = m.winding_turns * m.mean_turn_length * m.pole_count # m

    copper_resistivity_curve = np.interp(temperature, [60,70,80,90,100,110,120,130,140,150,200,250,260,270,280,290,300,310,320,330,340,350,400,450,500], # ohm-m - copper resistivity vs temperature
                                   [5.426,6.968,8.335,9.506,10.493,11.322,12.021,12.616,13.129,13.577,15.209,16.199,16.328,16.431,16.509,16.560,16.583,16.577,16.542,16.476,16.379,16.252,15.165,13.420,11.238]) 
        # NIST copper OFHC; 5% error at >= 50K
    winding_resistivity = np.interp(m.temperature, copper_resistivity_curve) # ohm-m - resistivity at operating temperature

    awg_wire_gauges = [7.35E-03,6.54E-03,5.83E-03,5.19E-03,4.62E-03,4.11E-03,3.66E-03,3.26E-03,2.91E-03,2.59E-03,2.30E-03,2.05E-03,1.83E-03,1.63E-03,1.45E-03,1.29E-03,1.15E-03,1.02E-03,9.12E-04,8.12E-04,7.23E-04,6.44E-04,5.73E-04,5.11E-04,4.55E-04,4.05E-04,3.61E-04,3.21E-04,2.86E-04,2.55E-04,2.27E-04,2.02E-04,1.80E-04,1.60E-04,1.43E-04,1.27E-04,1.13E-04,1.01E-04,8.97E-05,7.99E-05]

    slot_width = 0.003
    slot_depth = 0.010

    wire_packing_factor = np.pi/4 # unitless - theoretical maximum is 0.907

    idx = (np.abs(awg_wire_gauges - slot_width/2*slot_depth/m.winding_turns)).argmin()
    wire_diameter = awg_wire_gauges[idx]

    wire_dc_resistance = winding_resistivity * winding_length_total / (np.pi/4*wire_diameter**2) # ohm

    copper_relative_permeability = 0.999994 # 
    skin_depth = np.sqrt(2*winding_resistivity/m.ang_frequency/copper_relative_permeability/vacuum_permeability) # m - Hanselman 9.49
    skin_ratio = wire_diameter/skin_depth # unitless - skin effect ratio - Hanselman 9.48
    f = (np.sinh(2*skin_ratio) + np.sin(2*skin_ratio)) / (np.cosh(2*skin_ratio) - np.cos(2*skin_ratio)) # unitless - Hanselman 9.46
    g = (np.sinh(2*skin_ratio) - np.sin(2*skin_ratio)) / (np.cosh(2*skin_ratio) - np.cos(2*skin_ratio)) # unitless - Hanselman 9.47

    wire_row_count = np.floor(slot_depth / wire_diameter) # unitless
    ac_resistance_multiplier = f + 2/3*(wire_row_count**2-1)*g # unitless - Hanselman 9.45
        
    winding_resistance = wire_dc_resistance * ac_resistance_multiplier # ohm
    winding_loss = m.current_rms**2 * winding_resistance # W

    m.motor_constant = m.torque_constant / np.sqrt(winding_loss) # N*m/sqrt(W) - K_m

    coil_self_inductance =
    coil_mutual_inductance = 
    airgap_inductance = 
    slot_leakage_inductance =
    end_turn_inductance = 

    m.inductance_coil = coil_self_inductance + coil_mutual_inductance + airgap_inductance + slot_leakage_inductance + end_turn_inductance # H - L_phase - Hanselman 4.31
    m.inductance_phase = m.coil_count*m.inductance_coil # H - assumes 

    ### Rotorthrust production
    m.tangential_flux_density = vacuum_permeability * m.airgap_shear_stress / m.airgap_flux_density   # T - armature (load) field; from sigma_xz = B_n*B_t/mu0
    m.axial_force = m.airgap_area * (m.airgap_flux_density**2 - m.tangential_flux_density**2) / (2.0 * vacuum_permeability)   # N - full sigma_zz traction - Maxwell stress tensor with E = 0

def synchronous_links(link: link):
    # bus_voltage        = config['machine']['bus_voltage']         # V - line-neutral rms design bus
    winding_factor     = baseunit(config['machine']['winding_factor'])      # unitless - fundamental k_w
    pullout_margin_min = baseunit(config['machine']['pullout_margin_min'])  # unitless

    motor = link.motor
    alternator = link.alternator
    cable = link.cable

    alternator.ang_frequency_electrical = 
    # Shared electrical frequency (identical on both machines by construction)
    assert np.isclose(alternator.ang_frequency, motor.ang_frequency), \
        "electrical frequency mismatch: pole lock is inconsistent"
    alternator.frequency_electrical = alternator.ang_frequency_electrical / (2.0 * np.pi) # Hz
    motor.frequency_electrical = motor.ang_frequency_electrical / (2.0 * np.pi) # Hz

    # Machine geometry / sizing
    size_machine(motor, winding_factor, 90)
    size_machine(alternator, winding_factor, 300)

    # Power flow through the link
    cable.power = motor.power / motor.efficiency       # W - real power crossing the cable
    alternator.power = cable.power / alternator.efficiency       # W - alternator mechanical input

    cable.current_rms = motor.current_rms # A
   
    cable.resistance = cable.power / cable.current_rms**2 # ohm
    cable.inductance # property of geometry
    cable.reactance = cable.ang_frequency * cable.inductance # unitless

    cable.impedance = cable.resistance + 1j*cable.reactance # ohm
    link.impedance = alternator.impedance + cable.impedance + motor.impedance # ohm

    load_angle = np.arcsin(link.reactance_perunit) # rad - pullout at delta of 90 deg  ->  margin = 1 / sin(delta).
    # pullout_margin = 1.0 / link.reactance_perunit # unitless
    # if pullout_margin < pullout_margin_min:
    #     raise ValueError(f"Pull-out margin {pullout_margin:.2f} < required {pullout_margin_min:.2f}. "
    #         f"The link will slip poles under any load perturbation.")

    # # Line current: reactance-only series link, E_alt ~ E_mot ~ bus_voltage.
    # #   |I| = |E_alt - E_mot*exp(-j*delta)| / X = 2 V sin(delta/2) / X
    # link.reactance= link.reactance_perunit * cable.base_impedance # ohm
    # cable.current = 2.0 * bus_voltage * np.sin(load_angle / 2.0) / link.reactance # A - line current through the link

    # # --- 5. Losses, heat split, efficiency ---------------------------
    # resistance_warm = cable.resistivity_warm * (1 - cable.cold_fraction) * cable.length / cable.area
    # resistance_cold = cable.resistivity_cold * cable.cold_fraction * cable.length / cable.area
    # cable.resistance = resistance_warm + resistance_cold

    # motor.copper_loss = 3.0 * cable.current**2 * motor.resistance_phase
    # alternator.copper_loss = 3.0 * cable.current**2 * alternator.resistance_phase
    # cable.loss        = 3.0 * cable.current**2 * cable.resistance

    # motor.core_loss = motor.core_loss_fraction * motor.power
    # alternator.core_loss = alternator.core_loss_fraction * cable.power

    # heat_into_lox  = motor.copper_loss + motor.core_loss + 3.0 * cable.current**2 * resistance_cold
    # heat_into_fuel = alternator.copper_loss + alternator.core_loss + 3.0 * cable.current**2 * resistance_warm
    # link.efficiency = motor.power / alternator.power

    print(f"Link: {link.efficiency*100:.1f}% eff, {pullout_margin:.2f} pullout margin, "
          f"{motor.shear_margin:.2f} motor shear margin, {alternator.shear_margin:.2f} alternator shear margin, "
          f"{motor.axial_stiffness:.1f} N/m motor axial stiffness, {alternator.axial_stiffness:.1f} N/m alternator axial stiffness, "
          f"{cable.current:.1f} A line current, {cable.loss/1e3:.1f} kW cable loss, {heat_into_lox/1e3:.1f} kW heat into LOX, {heat_into_fuel/1e3:.1f} kW heat into fuel")
