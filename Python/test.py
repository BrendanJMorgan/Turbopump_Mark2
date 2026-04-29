import numpy as np
from rocketcea.cea_obj import CEA_Obj as rcea
import rocketcea.py_cea as py_cea
from rocketcea.units import add_user_units 
from engine_state import engine, tca

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
                                        
print(py_cea.comp)
print("~~~~~~~~~~~~~~")
print(py_cea.prtout.totn)