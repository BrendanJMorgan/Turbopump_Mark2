import numpy as np
from rocketcea.cea_obj import CEA_Obj as rcea
import rocketcea.py_cea as py_cea
from rocketcea.units import add_user_units 
from engine_state import engine, gg
from types import SimpleNamespace as _ns

def gas_generator():

    gg.p_exit = engine.p_amb

    gg.mdot_fuel = gg.mdot*(1/(1+gg.mixture_ratio))                 # kg/s - Fuel Mass Flow Rate
    gg.mdot_ox = gg.mdot*(gg.mixture_ratio/(1+gg.mixture_ratio))    # kg/s - Oxidizer Mass Flow Rate  

    # NASA CEA run ------------------------------------------------
    cea_gg = rcea(oxName=engine.oxidizer, fuelName=engine.fuel)  # NASA CEA Run
    cea_gg.get_full_cea_output(Pc=gg.pc/1E5, MR=gg.mixture_ratio, eps=0, subar=None, PcOvPe=gg.pc/gg.p_exit,
                                frozen=0, frozenAtThroat=0, short_output=0, show_transport=1,
                                pc_units='bar', output='calories', show_mass_frac=False)
    gg.Tt = py_cea.prtout.ttt[0]                                    # K - stagnation temperature at chamber
    if gg.Tt > 1100:
        raise ValueError(f"Warning: Chamber temperature of {gg.Tt} K may be too high for typical materials\n")
    gg.c_star_ideal = float(py_cea.rockt.cstr) / 3.28084 # m/s (CEA reports c* in ft/s) - Characteristic Velocity from infinite area combustor
    gg.Ae_At = py_cea.rockt.aeat[3]
    gg.gamma_avg = (py_cea.prtout.gammas[2] + py_cea.prtout.gammas[3]) / 2 # Ratio of Specific Heats, estimate
    gg.R_gas = 8314.5 / (py_cea.prtout.wm[2])  # J/kg-K (wm is g/mol, so R_u in J/kmol-K) - Specific Gas Constant (throat)
    
    # Performance Parameters ------------------------------------------------
    gg.c_star   = gg.c_star_ideal * gg.c_star_efficiency    # m/s - Characteristic Velocity
    gg.A_throat = gg.mdot * np.sqrt(gg.R_gas*gg.Tt/gg.gamma_avg) / gg.pc * (0.5*(gg.gamma_avg+1))**(0.5*(gg.gamma_avg+1.0)/(gg.gamma_avg-1.0))  # m2 - throat area (CEA form)


