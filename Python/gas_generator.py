import numpy as np
from rocketcea.cea_obj import CEA_Obj as rcea
from engine_state import engine, gg
from types import SimpleNamespace as _ns

def gas_generator():
    
    # NASA CEA run ------------------------------------------------
    cea_gg = rcea(oxName=engine.oxidizer, fuelName=engine.fuel)  # NASA CEA Run
    gg.p_exit = engine.p_amb  # Pa - exit pressure - assumes perfectly expanded flow

    gg.Ae_At = cea_gg.get_eps_at_PcOvPe(Pc=gg.pc/6894.76, MR=gg.OF, PcOvPe=gg.pc/gg.p_exit, frozen=0)    # unitless - Nozzle Expansion / Area Ratio
    gg.mw = 0.001*cea_gg.get_exit_MolWt_gamma(Pc=gg.pc/6894.76, MR=gg.OF, eps=gg.Ae_At)[0]             # kg/mol - molecular weight (fully combusted)
    gg.R = 8.3145 / (gg.mw / 1000.0)  # J/(kg*K) - specific gas constant

    gg.gamma = 1*cea_gg.get_exit_MolWt_gamma(Pc=gg.pc/6894.76, MR=gg.OF, eps=gg.Ae_At)[1]       # Ratio of Specific Heats, exit

    # Performance Parameters ------------------------------------------------
    [gg.isp_ideal,b,c] = cea_gg.get_IvacCstrTc(Pc=gg.pc/6894.76, MR=gg.OF, eps=gg.Ae_At) # s - specific impulse
    gg.c_star_ideal = b/3.28084 # m/s - Characteristic Velocity
    gg.Tc = c/1.8 # K - chamber temperature
    gg.p_exit = engine.p_amb                                                                # Pa - exit pressure - assumes perfectly expanded flow
    gg.v_exhaust_ideal = gg.isp_ideal * engine.g                                            # m/s - ideal exhaust velocity
    gg.c_star   = gg.c_star_ideal * gg.c_star_eff    # m/s - Characteristic Velocity

    gg.A_throat = gg.mdot * np.sqrt(gg.R*gg.Tc/gg.gamma) / gg.pc * (0.5*(gg.gamma+1))**(0.5*(gg.gamma+1.0)/(gg.gamma-1.0))  # m2 - throat area (CEA form)

    gg.p_injector = gg.pc * cea_gg.get_Pinj_over_Pcomb(Pc=gg.pc/6894.76, MR=gg.OF, fac_CR=np.pi*gg.r_chamber**2 / gg.A_throat) # Pa - injector dp criterion

    gg.density = cea_gg.get_Chamber_Density(Pc=gg.pc/6894.76, MR=gg.OF, eps=gg.Ae_At)  # kg/m3 - approximate density at chamber conditions


