from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from engine_state import pump, turbine, link, motor, alternator, cable, config

vacuum_permeability = 4.0e-7 * np.pi  # H/m

def size_machine(m, bus_voltage, winding_factor):

    ### Requirement
    m.torque = m.power / m.ang_frequency    # N-m
    m.airgap_shear_stress = 2*np.pi/(3*m.torque) * (m.r_outer**3 - m.r_inner**3) # Pa - sigma - Hanselman 4.43 - inegrated over annulus; assumes constant flux density

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
    m.torque_constant = m.pole_count*m.winding_turns*m.airgap_flux_density*(m.r_outer**2 - m.r_inner**2) # N*m - K_t - Hanselman 4.11
    m.current_rms = m.torque / m.torque_constant # A - i - Hanselman 4.11
    m.back_emf_rms = m.torque_constant * m.ang_frequency # V - E_b - Hanselman 4.10

    ### Power efficiency

    ### Rotorthrust production
    m.tangential_flux_density = vacuum_permeability * m.airgap_shear_stress / m.airgap_flux_density   # T - armature (load) field; from sigma_xz = B_n*B_t/mu0
    m.axial_force = m.airgap_area * (m.airgap_flux_density**2 - m.tangential_flux_density**2) / (2.0 * vacuum_permeability)   # N - full sigma_zz traction - Maxwell stress tensor with E = 0

def synchronous_links(link: link):
    # bus_voltage        = config['machine']['bus_voltage']         # V - line-neutral rms design bus
    winding_factor     = config['machine']['winding_factor']      # unitless - fundamental k_w
    pullout_margin_min = config['machine']['pullout_margin_min']  # unitless

    motor = link.motor
    alternator = link.alternator
    cable = link.cable

    # Shared electrical frequency (identical on both machines by construction)
    assert np.isclose(alternator.ang_frequency, motor.ang_frequency), \
        "electrical frequency mismatch: pole lock is inconsistent"
    alternator.frequency = alternator.ang_frequency / (2.0 * np.pi) # Hz
    motor.frequency = motor.ang_frequency / (2.0 * np.pi) # Hz

    # Machine geometry / sizing
    size_machine(motor, winding_factor)
    size_machine(alternator, winding_factor)

    # Power flow through the link
    cable.power = motor.power / motor.efficiency       # W - real power crossing the cable
    alternator.power = cable.power / alternator.efficiency       # W - alternator mechanical input
   
    # cable.base_impedance    = bus_voltage**2 / cable.power # ohm
    # cable.reactance_perunit = alternator.ang_frequency * cable.inductance / cable.base_impedance # unitless
    # link.reactance_perunit = alternator.reactance_perunit + cable.reactance_perunit + motor.reactance_perunit # ohm - P = 3 V^2 sin(delta) / X   ->  sin(delta) = x_total_pu

    # if link.reactance_perunit >= 1.0:
    #     raise ValueError(f"Series reactance {link.reactance_perunit:.2f} pu >= 1: the link cannot transmit "
    #         f"rated power at any angle. Raise bus voltage (fewer turns) or split the power.")
    
    # load_angle = np.arcsin(link.reactance_perunit) # rad - pullout at delta of 90 deg  ->  margin = 1 / sin(delta).
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
