import numpy as np

from fullflow import *
from fullplot import Axis, generate_map
from thermoprop import *

from engine_state import engine, tca, gg, pump, turbine, motor, alternator, cable, link

import pint
ureg = pint.UnitRegistry()


mark2network = Network("Mark 2 Network")

ambient = 93900 * ureg.pascal

gg.pc = 50.0 * ureg.bar

# Liquid Oxygen Pump
ox_pump = pump(fluid='ox')
ox_pump.p_out = 100.0 * ureg.bar                   # Pa - add 10 psi of margin for plumbing losses
ox_pump.p_in  = 5.0 * ureg.bar                                           # Pa - inlet pressure
ox_pump.T_in  = rprop_oxidizer.TdegRAtPsat(engine.p_amb/6894.76)/1.8     # K - tank temperature - assumes saturated propellant during fill
ox_pump.density = 1000*rprop_oxidizer.SGLiqAtTdegR(ox_pump.T_in*1.8)      # kg/m3 - density of LOX at inlet
ox_pump.viscosity = 0.1*rprop_oxidizer.ViscAtTdegR(ox_pump.T_in*1.8)     # Pa-s - viscosity of LOX at inlet
ox_pump.mdot = engine.mdot_ox_total                                      # kg/s     
ox_pump.pvap_inlet = 6894.76*rprop_oxidizer.PvapAtTdegR(ox_pump.T_in*1.8)        # Pa - vapor pressure of LOX at tank temperature

# Fuel Pump
fuel_pump = pump(fluid='fuel')
fuel_pump.p_out = 100.0 * ureg.bar                   # Pa - add 10 psi of margin for plumbing losses
fuel_pump.p_in = 5.0 * ureg.bar                                                      # Pa
fuel_pump.T_in = engine.T_amb                                                       # K - tank temperature    
fuel_pump.density   = 1000*rprop_fuel.SGLiqAtTdegR(fuel_pump.T_in*1.8)      # kg/m3 - density of fuel at inlet
fuel_pump.viscosity = 0.1*rprop_fuel.ViscAtTdegR(fuel_pump.T_in*1.8)     # Pa-s - viscosity of fuel at inlet
fuel_pump.mdot  = engine.mdot_fuel_total                                            # kg/s
fuel_pump.pvap_inlet = 6894.76*rprop_fuel.PvapAtTdegR(fuel_pump.T_in*1.8) # Pa - vapor pressure of RP1 at tank temperature


turbine.p_downstream = ambient
fuel_pump.shaft_speed = 30000 * ureg.rpm
fuel_pump.mdot = 5.0 * ureg.kilogram / ureg.second
fuel_pump.p_inlet = 5.0 * ureg.bar
fuel_pump.p_outlet = 100.0 * ureg.bar
fuel_pump.efficiency = 0.72

ox_pump.mdot = 5.0 * ureg.kilogram / ureg.second
ox_pump.p_inlet = 5.0 * ureg.bar
ox_pump.p_outlet = 100.0 * ureg.bar
ox_pump.shaft_speed = 20000 * ureg.rpm

ox_pump.efficiency = 0.72
turbine.efficiency = 0.55

