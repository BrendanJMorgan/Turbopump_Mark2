import numpy as np
import math
from rocketcea.cea_obj import CEA_Obj as rcea
from engine_state import engine, tca

def combustion_chamber():   

    cea_cc = rcea(oxName=engine.oxidizer, fuelName=engine.fuel) # NASA CEA Run

    tca.Ae_At = cea_cc.get_eps_at_PcOvPe(Pc=tca.pc/6894.76, MR=tca.OF, PcOvPe=tca.pc/engine.p_amb, frozen=0) # unitless - Nozzle Expansion Ratio

    gamma_avg_nozzle = np.mean([cea_cc.get_Throat_MolWt_gamma(Pc=tca.pc/6894.76, MR=tca.OF, eps=tca.Ae_At)[1], 
      cea_cc.get_exit_MolWt_gamma(Pc=tca.pc/6894.76, MR=tca.OF, eps=tca.Ae_At)[1]]); # Ratio of Specific Heats, estimate
    tca.gamma_avg_nozzle = gamma_avg_nozzle
    tca.R_gas = 8.3145 / (0.001*cea_cc.get_Throat_MolWt_gamma(Pc=tca.pc/6894.76, MR=tca.OF, eps=tca.Ae_At)[0]);  # J/kg-K - Specific Gas Constant (throat)

    # Performance Parameters -----------------------------------------------------------------------------------
    [tca.isp_ideal,b,c] = cea_cc.get_IvacCstrTc(Pc=tca.pc/6894.76, MR=tca.OF, eps=tca.Ae_At) # s - specific impulse
    tca.c_star_ideal = b/3.28084 # m/s - Characteristic Velocity
    tca.Tt_cc = c/1.8 # K - chamber temperature
    tca.p_exit = engine.p_amb; # Pa - exit pressure - assumes perfectly expanded flow
    tca.v_exhaust_ideal = tca.isp_ideal*engine.g; # m/s - ideal exhaust velocity

    lambda_cone = 0.5*(1+np.cos(tca.diverge_angle))  # Divergence correction factor for conical nozzle

    tca.c_star = tca.c_star_ideal * tca.c_star_eff     # m/s - Characteristic Velocity
    c_tau_ideal = np.sqrt( (2*gamma_avg_nozzle**2/(gamma_avg_nozzle-1) * (2/(gamma_avg_nozzle+1))**((gamma_avg_nozzle+1)/(gamma_avg_nozzle-1)) * \
                    (1-(tca.p_exit/tca.pc)**((gamma_avg_nozzle-1)/gamma_avg_nozzle) ) ) ) + tca.Ae_At*(tca.p_exit-engine.p_amb)/tca.pc; # Ideal Thrust Coefficient
    tca.c_tau = c_tau_ideal*lambda_cone*tca.c_tau_eff                # Thrust Coefficient
    tca.A_throat = engine.thrust / (tca.pc * tca.c_tau * tca.c_star_eff)   # m2 - Throat Area

    p_injector = tca.pc*cea_cc.get_Pinj_over_Pcomb(Pc=tca.pc/6894.76, MR=tca.OF, fac_CR=math.pi*tca.r1_chamber**2 / tca.A_throat) # Pa
    v_exhaust = tca.c_star * tca.c_tau  # m/s - Exhaust Velocity

    # Mass Flow Rates ----------------------------------------------------------------------------------------------
    tca.mdot = tca.A_throat / (np.sqrt(tca.Tt_cc)/tca.pc * np.sqrt(tca.R_gas/gamma_avg_nozzle) * 
                          ((gamma_avg_nozzle+1)/2)**((gamma_avg_nozzle+1)/(2*(gamma_avg_nozzle-1)))); # m2 - throat area
    tca.mdot_fuel = tca.mdot*(1/(1+tca.OF));            # kg/s - Fuel Mass Flow Rate
    tca.mdot_ox = tca.mdot*(tca.OF/(1+tca.OF));   # kg/s - Oxidizer Mass Flow Rate