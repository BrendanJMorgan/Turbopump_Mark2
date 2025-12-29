import numpy as np

from combustion_chamber import combustion_chamber
from cc_geometry import cc_geometry
from coolant_flow import coolant_flow
# from thermal_balance import thermal_balance  # assumes it returns updated T_wall_hot
import rocketcea.py_cea as py_cea

def thrust_chamber_assembly():
    combustion_chamber()
    cc_geometry()
    coolant_flow()

    # if inputs.compute_thermals == True:
    #     T_wall_hot = 800 * np.ones(len(inputs.x)) # K - Initial guess for the hot wall temperatures 
    #     thermal_balance(inputs, cc_outputs, geometry_outputs, cc_gas_outputs, coolant_outputs) # Perform thermal balance iteration
    #     thermal_balance(inputs, cc_outputs, geometry_outputs, cc_gas_outputs, coolant_outputs) # Run thermal balance a second time for greater convergence 

    # Clear CEA output arrays to prevent data carryover between runs
    for attr in dir(py_cea.prtout):
        obj = getattr(py_cea.prtout, attr)
        if isinstance(obj, np.ndarray) and obj.dtype.kind == 'f':
            obj[:] = np.nan


