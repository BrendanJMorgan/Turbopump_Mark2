from dataclasses import dataclass
import numpy as np
import math
from engine_state import engine, tca, gg, pump, turbine, link, config, baseunit

from pumps import pumps
from mainline import mainline
from gas_generator import gas_generator
from turbines import turbines
# from synchronous_links import synchronous_links

from CoolProp.CoolProp import PropsSI
from rocketprops.rocket_prop import get_prop as rprop

def powerhead():

    rprop_fuel = rprop(engine.fuelrp)
    rprop_oxidizer = rprop(engine.oxidizer)

    # Total Mass Flow Rate INITIAL GUESSES
    gg.mdot = tca.mdot*gg.fraction_guess/(1-gg.fraction_guess)                  # kg/s - mass flow rate through gas generator
    gg.mdot_fuel = gg.mdot*(1/(1+gg.mixture_ratio))                             # kg/s - fuel mass flow rate through gas generator
    gg.mdot_ox   = gg.mdot*(gg.mixture_ratio/(1+gg.mixture_ratio))              # kg/s - oxidizer mass flow rate through gas generator
    engine.mdot_fuel_total = tca.mdot_fuel + gg.mdot_fuel                       # kg/s - total fuel mass flow rate
    engine.mdot_ox_total   = tca.mdot_ox   + gg.mdot_ox                         # kg/s - total oxidizer mass flow rate
    # TODO: add a verification and/or iteration on the gg fraction guess at end of powerhead module

    # Liquid Oxygen Pump
    ox_pump = pump(fluid='ox')
    ox_pump.p_out = tca.p_injector_manifold + ox_pump.plumbing_loss             # Pa - add margin for plumbing losses
    ox_pump.p_in  = ox_pump.tank.p                                              # Pa   
    ox_pump.T_in  = rprop_oxidizer.TdegRAtPsat(engine.p_amb/6894.76)/1.8        # K - tank temperature - assumes saturated propellant during fill
    ox_pump.density = 1000*rprop_oxidizer.SGLiqAtTdegR(ox_pump.T_in*1.8)        # kg/m3 - density of LOX at inlet
    ox_pump.viscosity = 0.1*rprop_oxidizer.ViscAtTdegR(ox_pump.T_in*1.8)        # Pa-s - viscosity of LOX at inlet
    ox_pump.mdot = engine.mdot_ox_total                                         # kg/s     
    ox_pump.pvap_inlet = 6894.76*rprop_oxidizer.PvapAtTdegR(ox_pump.T_in*1.8)   # Pa - vapor pressure of LOX at tank temperature
  
    # Fuel Pump
    fuel_pump = pump(fluid='fuel')
    fuel_pump.p_in = fuel_pump.tank.p                                           # Pa
    fuel_pump.T_in = engine.T_amb                                               # K - tank temperature    
    fuel_pump.density   = 1000*rprop_fuel.SGLiqAtTdegR(fuel_pump.T_in*1.8)      # kg/m3 - density of fuel at inlet
    fuel_pump.viscosity = 0.1*rprop_fuel.ViscAtTdegR(fuel_pump.T_in*1.8)        # Pa-s - viscosity of fuel at inlet
    fuel_pump.mdot  = engine.mdot_fuel_total                                    # kg/s
    fuel_pump.pvap_inlet = 6894.76*rprop_fuel.PvapAtTdegR(fuel_pump.T_in*1.8)   # Pa - vapor pressure of RP1 at tank temperature
    fuel_pump.mainline.p_out = tca.regenerative_coolant.p[0 if tca.regenerative_coolant.flow_direction == 1 else -1] # Pa - regen jacket inlet
    fuel_pump.p_out = mainline(fuel_pump) # Pa
    
    # Find Shaft Speed suction limited caps
    ox_pump.shaft_speed = fuel_pump.shaft_speed = None
    pumps(ox_pump)
    pumps(fuel_pump)
    # Each pump's own max allowable speed, clipped to the mechanical cap
    ox_cap   = min(ox_pump.shaft_speed,   engine.shaft_speed_cap)
    fuel_cap = min(fuel_pump.shaft_speed, engine.shaft_speed_cap)

    max_pole_pairs = baseunit(config['machine']['max_pole_pairs'])  # max pole pairs to allow. Set to 1 for a common shaft
    best_utilization = (-1.0, 1, 1, 0.0)
    for a in range(1, max_pole_pairs + 1):
        for b in range(1, max_pole_pairs + 1):
            if math.gcd(a, b) != 1:
                continue
            ang_frequency_electric = min(ox_cap*a, fuel_cap*b)          # base cable frequency
            utilization  = min(ang_frequency_electric/a / ox_cap, ang_frequency_electric/b / fuel_cap)
            if utilization > best_utilization[0]:
                best_utilization = (utilization, a, b, ang_frequency_electric)

    _, pole_pairs_ox, pole_pairs_fuel, ang_frequency_electric = best_utilization
    ox_pump.shaft_speed   = ang_frequency_electric/pole_pairs_ox     # rad/s - omega_mech = omega_elec / pole pairs
    fuel_pump.shaft_speed = ang_frequency_electric/pole_pairs_fuel   # rad/s

    # Actual Runs
    pumps(ox_pump)
    pumps(fuel_pump)

    # Re-pick the fuel mainline now that the volute throat is known (nozzle mass and loss), then re-run the pump
    for iteration in range(5):
        old_tube = fuel_pump.mainline.tube_id
        fuel_pump.p_out = mainline(fuel_pump) # Pa
        pumps(fuel_pump)
        if fuel_pump.mainline.tube_id == old_tube:
            break
    else:
        print(f"Fuel mainline tube choice did not settle (last {fuel_pump.mainline.tube_id*1000:.1f} mm ID)")

    powerhead_link = link()
    powerhead_link.motor = link.motor(fluid=ox_pump.fluid, poles=2*pole_pairs_ox, ang_frequency=ox_pump.shaft_speed, power=ox_pump.shaft_power, r_max=0.95*ox_pump.impeller[0].r_outlet)
    powerhead_link.alternator = link.alternator(fluid=fuel_pump.fluid, poles=2*pole_pairs_fuel, ang_frequency=fuel_pump.shaft_speed, r_max=0.95*fuel_pump.impeller[0].r_outlet)
    powerhead_link.cable = link.cable(frequency=ang_frequency_electric/(2.0*np.pi))
    powerhead_link.cable.ang_frequency_electric = ang_frequency_electric
    powerhead_link.motor.pole_count = 2*pole_pairs_ox
    powerhead_link.alternator.pole_count = 2*pole_pairs_fuel

    powerhead_link.motor.r_outer=0.95*ox_pump.impeller[0].r_outlet
    powerhead_link.alternator.r_outer=0.95*fuel_pump.impeller[0].r_outlet

    # synchronous_links(powerhead_link)

    gg.pc = np.min([fuel_pump.p_out, ox_pump.p_out]) * gg.stiffness  # Pa - combustion chamber pressure guess
    gas_generator()

    # Turbine
    if engine.common_shaft == True:
        common_turbine = turbine()
        # TODO: power loss into bearings and seals
        common_turbine.shaft_efficiency = 1 # unitless - efficiency of power transmission through the shaft(s)
        common_turbine.shaft_power = (ox_pump.shaft_power+fuel_pump.shaft_power)/common_turbine.shaft_efficiency  
            # W - power that must be produced by the turbine
        common_turbine.shaft_speed = fuel_pump.shaft_speed  # rad/s - both pumps on same shaft
        turbines(common_turbine)

    else:
        ox_turbine = turbine()
        fuel_turbine = turbine()
        ox_turbine.shaft_efficiency = turbine.gear_efficiency # unitless - efficiency of power transmission through the shaft(s)
        fuel_turbine.shaft_efficiency = turbine.gear_efficiency # unitless - efficiency of power transmission through the shaft(s)
        ox_turbine.shaft_power = ox_pump.shaft_power/ox_turbine.shaft_efficiency  # W - power that must be produced by the turbine
        fuel_turbine.shaft_power = fuel_pump.shaft_power/fuel_turbine.shaft_efficiency  # W - power that must be produced by the turbine
        turbines(turbine=ox_turbine)
        turbines(turbine=fuel_turbine)


    powerhead.ox_pump = ox_pump
    powerhead.fuel_pump = fuel_pump
    powerhead.link = powerhead_link
    if engine.common_shaft == True:
        powerhead.turbine = common_turbine
    else:
        powerhead.ox_turbine = ox_turbine
        powerhead.fuel_turbine = fuel_turbine
