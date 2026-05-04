import numpy as np

from combustion_chamber import combustion_chamber
from combustion_chamber_geometry import combustion_chamber_geometry
from combustion_chamber_cooling import combustion_chamber_cooling
from coolant_flow import coolant_flow
from thermal_balance import thermal_balance
# from thermal_balance import thermal_balance  # assumes it returns updated T_wall_hot
import rocketcea.py_cea as py_cea

from engine_state import tca

def thrust_chamber_assembly():
    combustion_chamber()
    combustion_chamber_geometry()
     
    combustion_chamber_cooling() # Perform thermal balance iteration

    # Clear CEA output arrays to prevent data carryover between runs
    for attr in dir(py_cea.prtout):
        obj = getattr(py_cea.prtout, attr)
        if isinstance(obj, np.ndarray) and obj.dtype.kind == 'f':
            obj[:] = np.nan


