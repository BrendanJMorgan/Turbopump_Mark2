import numpy as np
from rocketcea.cea_obj import CEA_Obj as rcea
from engine_state import engine, tca, gg
from types import SimpleNamespace as _ns

def gas_generator():

    # NASA CEA run ------------------------------------------------
    cea_gg = rcea(oxName=engine.oxidizer, fuelName=engine.fuel)  # NASA CEA Run

    gg.Ae_At = cea_gg.get_eps_at_PcOvPe(Pc=gg.pc/6894.76, MR=gg.OF, PcOvPe=gg.pc/engine.p_amb, frozen=0)            # unitless - Nozzle Expansion / Area Ratio
    gg.mw = 0.001*cea_gg.get_Chamber_MolWt_gamma(Pc=gg.pc/6894.76, MR=gg.OF, PcOvPe=gg.pc/engine.p_amb, frozen=0)   # kg/mol - molecular weight (fully combusted)
    gg.R = 8.3145 / (gg.mw / 1000.0)  # J/(kg*K) - specific gas constant

    gg.gamma_nozzle = np.mean([
        cea_gg.get_Throat_MolWt_gamma(Pc=gg.pc/6894.76, MR=gg.OF, eps=gg.Ae_At)[1],     # Ratio of Specific Heats, throat
        cea_gg.get_exit_MolWt_gamma(Pc=gg.pc/6894.76, MR=gg.OF, eps=gg.Ae_At)[1],       # Ratio of Specific Heats, exit
    ])

    # Performance Parameters ------------------------------------------------
    [gg.isp_ideal, _] = cea_gg.get_IvacCstrTc(Pc=gg.pc/6894.76, MR=gg.OF, eps=gg.Ae_At)     # s - specific impulse (vacuum curve-fit)
    gg.c_star_ideal = cea_gg.get_Cstar(Pc=gg.pc/6894.76, MR=gg.OF, eps=gg.Ae_At)            # m/s - Characteristic Velocity (match your style)
    gg.Tt_gg = cea_gg.get_Tc(Pc=gg.pc/6894.76, MR=gg.OF)                                    # K - chamber temperature
    gg.p_exit = engine.p_amb                                                                # Pa - exit pressure - assumes perfectly expanded flow
    gg.v_exhaust_ideal = gg.isp_ideal * engine.g                                            # m/s - ideal exhaust velocity

    lambda_cone = 0.5 * (1 + np.cos(gg.diverge_angle))  # Divergence correction factor for conical nozzle
    gg.c_star   = gg.c_star_ideal * gg.c_star_eff_gg    # m/s - Characteristic Velocity

    
    c_tau_ideal = (np.sqrt((2 * (gg.gamma_nozzle**2) / (gg.gamma_nozzle - 1)) * ((2 / (gg.gamma_nozzle + 1))**((gg.gamma_nozzle + 1)/(gg.gamma_nozzle - 1)) *
        (1 - (gg.p_exit/gg.pc)**((gg.gamma_nozzle - 1)/gg.gamma_nozzle)))) + gg.Ae_At * (gg.p_exit - engine.p_amb) / gg.pc)
        # unitless - Ideal Thrust Coefficient

    gg.c_tau = c_tau_ideal * lambda_cone * gg.c_tau                               # unitless - Thrust Coefficient
    gg.A_throat = gg.mdot * gg.R_gas * np.sqrt(gg.Tt_gg) / (gg.pc * np.sqrt(gg.gamma_nozzle) *
                    (( (gg.gamma_nozzle+1)/2 )**(-(gg.gamma_nozzle+1)/(2*(gg.gamma_nozzle-1)))))  # m2 - throat area (CEA form)

    p_cc_injector = gg.pc * cea_gg.get_Pinj_over_Pcomb(Pc=gg.pc/6894.76, MR=gg.OF,
                        fac_CR=np.pi*gg.r1_chamber**2 / gg.A_throat)                  # Pa - injector dp criterion
    gg.v_exhaust = gg.c_star * gg.c_tau                                                        # m/s - Exhaust Velocity



