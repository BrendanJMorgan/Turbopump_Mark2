import numpy as np

from combustion_chamber import combustion_chamber
from combustion_chamber_geometry import combustion_chamber_geometry
from combustion_chamber_cooling2 import combustion_chamber_cooling2
from coolant_flow import coolant_flow
from thermal_balance import thermal_balance
# from thermal_balance import thermal_balance  # assumes it returns updated T_wall_hot
import rocketcea.py_cea as py_cea

from engine_state import tca

def thrust_chamber_assembly():
    combustion_chamber()
    combustion_chamber_geometry()
    tca.regenerative_jacket.T_hot = 0.5 * tca.Tt * np.ones(len(tca.x)) # K - Initial guess for the hot wall temperatures 
    combustion_chamber_cooling2() # Perform thermal balance iteration

    # Clear CEA output arrays to prevent data carryover between runs
    for attr in dir(py_cea.prtout):
        obj = getattr(py_cea.prtout, attr)
        if isinstance(obj, np.ndarray) and obj.dtype.kind == 'f':
            obj[:] = np.nan


