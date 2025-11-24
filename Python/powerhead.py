from dataclasses import dataclass
import numpy as np
from engine_state import engine, tca, gg, pump, turbine

from pumps import pumps
from gas_generator import gas_generator
from turbines import turbines

from CoolProp.CoolProp import PropsSI
from rocketprops.rocket_prop import get_prop as rprop

def powerhead():
    # Total Mass Flow Rate INITIAL GUESSES
    gg.mdot = tca.mdot*gg.fraction_guess/(1-gg.fraction_guess)   # kg/s - mass flow rate through gas generator
    gg.mdot_fuel = gg.mdot*(1/(1+gg.OF))                                # kg/s - fuel mass flow rate through gas generator
    gg.mdot_ox   = gg.mdot*(gg.OF/(1+gg.OF))                         # kg/s - oxidizer mass flow rate through gas generator
    engine.mdot_fuel_total = tca.mdot_fuel + gg.mdot_fuel            # kg/s - total fuel mass flow rate
    engine.mdot_ox_total   = tca.mdot_ox   + gg.mdot_ox              # kg/s - total oxidizer mass flow rate
    # TODO: add a verification and/or iteration on the gg fraction guess at end of powerhead module

    # Liquid Oxygen Pump
    ox_pump = pump(fluid='ox')
    ox_pump.p_out = tca.pc*(1+tca.stiffness) + 10*6894.76                                           # Pa - add 10 psi of margin for plumbing losses
    ox_pump.p_in  = engine.p_amb                                                                    # Pa   
    ox_pump.T_in  = PropsSI("T","P",ox_pump.p_in,"Q",0,"Oxygen")                                    # K - tank temperature
    ox_pump.density   = PropsSI("D","Q",0,"P",ox_pump.p_in,"Oxygen")                                # kg/m3 - density of LOX at inlet
    ox_pump.mdot  = engine.mdot_ox_total                                                            # kg/s     
    ox_pump.pvap_inlet = PropsSI("P","T",ox_pump.T_in,"Q",0,"Oxygen")                           # Pa - vapor pressure of LOX at tank temperature
  
    print(ox_pump.clocking)
    pumps(ox_pump) 

    # Fuel Pump
    fuel_pump = pump(fluid='fuel')
    fuel_pump.p_out = fuel_pump.p_cool_1 + 10*6894.76                                               # Pa - add 10 psi of margin for plumbing losses
    fuel_pump.p_in  = fuel_pump.p_amb                                                               # Pa
    fuel_pump.T_in  = fuel_pump.T_amb                                                               # K - tank temperature    
    fuel_pump.density   = 1000*rprop(engine.fuelrp).SGLiqAtTdegR(engine.T_amb*1.8)                  # kg/m3 - density of fuel at inlet
    fuel_pump.mdot  = engine.mdot_fuel_total                                                        # kg/s
    fuel_pump.pvap_inlet = 0.001*rprop(engine.fuelrp).rp1_pvap_at_TdegR(fuel_pump.T_amb*1.8)    # Pa - vapor pressure of RP1 at tank temperature
    
    pumps(fuel_pump)

    gas_generator()

    # Turbine
    if engine.common_shaft == True:
        common_turbine = turbine()
        # TODO: power loss into bearings and seals
        common_turbine.shaft_efficiency = turbine.gear_efficiency # unitless - efficiency of power transmission through the shaft(s)
        common_turbine.shaft_power = (ox_pump.shaft_power+fuel_pump.shaft_power)/common_turbine.shaft_efficiency  
            # W - power that must be produced by the turbine
        turbines(turbine=common_turbine)
    else:
        ox_turbine = turbine()
        fuel_turbine = turbine()
        ox_turbine.shaft_efficiency = turbine.gear_efficiency # unitless - efficiency of power transmission through the shaft(s)
        fuel_turbine.shaft_efficiency = turbine.gear_efficiency # unitless - efficiency of power transmission through the shaft(s)
        ox_turbine.shaft_power = ox_pump.shaft_power/ox_turbine.shaft_efficiency  # W - power that must be produced by the turbine
        fuel_turbine.shaft_power = fuel_pump.shaft_power/fuel_turbine.shaft_efficiency  # W - power that must be produced by the turbine
        turbines(turbine=ox_turbine)
        turbines(turbine=fuel_turbine)
