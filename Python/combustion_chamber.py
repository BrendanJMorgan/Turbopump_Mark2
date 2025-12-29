import numpy as np
from rocketcea.cea_obj import CEA_Obj as rcea
import rocketcea.py_cea as py_cea
from rocketcea.units import add_user_units 
from engine_state import engine, tca

def combustion_chamber():   
    tca.p_exit = engine.p_amb  # Pa - exit pressure - assumes perfectly expanded flow
    gamma_guess = 1.11 # Initial guess for gamma
    c_tau_guess = tca.c_tau_eff*np.sqrt( (2*gamma_guess**2/(gamma_guess-1) * (2/(gamma_guess+1))**((gamma_guess+1)/(gamma_guess-1)) * (1-(tca.p_exit/tca.pc)**((gamma_guess-1)/gamma_guess) ) ) )  # Initial guess for thrust coefficient 
    tca.A_throat = engine.thrust / (tca.pc*c_tau_guess)  # m2 - Throat Area
    tca.r_throat = np.sqrt(tca.A_throat/np.pi)  # m - Throat Radius

    cea_cc = rcea(oxName=engine.oxidizer, fuelName=engine.fuel) # NASA CEA Run
    cea_cc.get_full_cea_output(Pc=tca.pc/1E5, MR=tca.OF, eps=0, subar=None, PcOvPe=tca.pc/tca.p_exit,
                            frozen=0, frozenAtThroat=0, short_output=0, show_transport=1,
                            pc_units='bar', output='calories', show_mass_frac=False,
                            fac_CR= (tca.r1_chamber/tca.r_throat)**2)  # Finite area combustor
    tca.c_star_ideal = float(py_cea.rockt.cstr) # m/s - Characteristic Velocity from infinite area combustor
    tca.Ae_At = py_cea.rockt.aeat[3]
    tca.gamma_avg_nozzle = (py_cea.prtout.gammas[2] + py_cea.prtout.gammas[3]) / 2 # Ratio of Specific Heats, estimate
    tca.R_gas = 8.3145 / (py_cea.prtout.wm[2])  # J/kg-K - Specific Gas Constant (throat)

    # Performance Parameters -----------------------------------------------------------------------------------
    # [tca.isp_ideal, tca.c_star_ideal, tca.Tt] = tca.cea_cc.get_IvacCstrTc(Pc=tca.pc, MR=tca.OF, eps=tca.Ae_At) # s, m/s, K - specific impulse, characteristic velocity, chamber temperature

    lambda_cone = 0.5*(1+np.cos(tca.diverge_angle))  # Divergence correction factor for conical nozzle

    tca.c_star = tca.c_star_ideal * tca.c_star_eff     # m/s - Characteristic Velocity
    c_tau_ideal = np.sqrt( (2*tca.gamma_avg_nozzle**2/(tca.gamma_avg_nozzle-1) * (2/(tca.gamma_avg_nozzle+1))**((tca.gamma_avg_nozzle+1)/(tca.gamma_avg_nozzle-1)) * \
                    (1-(tca.p_exit/tca.pc)**((tca.gamma_avg_nozzle-1)/tca.gamma_avg_nozzle) ) ) ) + tca.Ae_At*(tca.p_exit-engine.p_amb)/tca.pc # Ideal Thrust Coefficient
    tca.c_tau = c_tau_ideal*lambda_cone*tca.c_tau_eff                # Thrust Coefficient
    tca.A_throat = engine.thrust / (tca.pc * tca.c_tau)  # m2 - Throat Area

    tca.v_exhaust = tca.c_star * tca.c_tau  # m/s - Exhaust Velocity

    # Mass Flow Rates ----------------------------------------------------------------------------------------------
    tca.Tt = py_cea.prtout.ttt[0]  # K - stagnation temperature at chamber
    tca.mdot = tca.A_throat / (np.sqrt(tca.Tt)/tca.pc * np.sqrt(tca.R_gas/tca.gamma_avg_nozzle) * 
                          ((tca.gamma_avg_nozzle+1)/2)**((tca.gamma_avg_nozzle+1)/(2*(tca.gamma_avg_nozzle-1)))) # m2 - throat area
    tca.mdot_fuel = tca.mdot*(1/(1+tca.OF))            # kg/s - Fuel Mass Flow Rate
    tca.mdot_ox = tca.mdot*(tca.OF/(1+tca.OF))   # kg/s - Oxidizer Mass Flow Rate

