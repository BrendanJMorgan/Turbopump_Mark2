## Objective: finds the temperature distributions of the combustion chamber and nozzle walls

from os import error
import numpy as np
import rocketcea.py_cea as py_cea

from engine_state import engine, tca
from mixture import mixture
from rocketcea.cea_obj import CEA_Obj as rcea
from rocketprops.rocket_prop import get_prop as rprop

def thermal_balance():

    gas = tca.combustion_gas
    jacket = tca.regenerative_jacket
    cool = tca.regenerative_coolant
    film = tca.film_coolant

    rprop_fuel = rprop(engine.fuelrp)

    # Combustion Thermal Environment
    gas.T_free = tca.Tt / (1 + (gas.gamma - 1) / 2 * gas.M**2)  # Free-Stream Temperature

    # Wall Conduction
    fin_thickness = jacket.gap_pipe*np.ones_like(tca.x)  # m
    fin_height = jacket.h_pipe     # m
    fin_length = tca.x[-1]      # m

    gas.A = 2 * np.pi * tca.r1 * tca.dx / jacket.n_pipe  # m2 - Area of tca.dx per coolant channel
    t_wall = tca.r2 - tca.r1                          # m - Wall thickness (projected vertically)

    # Thermal Balance (preallocations)
    # products = engine.products
    # fractions = engine.fractions_at(tca.x)  # shape: (n_species, len(x))

    jacket.T_cold=np.zeros_like(tca.x)            # K - wall temperature on coolant side
    cool.T =    engine.T_amb*np.ones_like(tca.x)   # K - regent coolant temperature
    gas.T_adiabatic_wall =engine.T_amb*np.ones_like(tca.x)   # K - adiabatic wall temperature
    gas.cp =    np.zeros_like(tca.x)            # J/kg-K - specific heat capacity (at constant pressure) of combustion gases
    gas.visc =  np.zeros_like(tca.x)            # Pa-s - dynamic viscosity of combustion gases
    gas.cond =  np.zeros_like(tca.x)            # W/m-K - thermal conductivity of combustion gases
    gas.h =     np.zeros_like(tca.x)            # W / m2-K - convective heat transfer coefficient
    gas.q =     np.zeros_like(tca.x)            # W - heat flux from combustion gases to wall
    gas.prandtl=np.zeros_like(tca.x)            # Prandtl number of combustion gases
    gas.dens =  np.zeros_like(tca.x)            # kg/m3 - density of combustion gases
    cool.nusselt =   np.zeros_like(tca.x)            # Nusselt number of coolant flow
    cool.h =    np.zeros_like(tca.x)            # W / m2-K - convective heat transfer coefficient of coolant flow

    warned = False # warning flag for laminar coolant flow

    
    gas.T_boundary = gas.T_free * (1 + 0.032 * gas.M**2 + 0.58 * (jacket.T_hot / gas.T_free - 1)) # Reference temperature for gas transport properties
    
    #######################################################################################################
    ### Gas Mixture Properties
    #######################################################################################################

    gas.products = []
    gas.mass_fractions = []
    for k in range(py_cea.indx.ngc):
        raw_name = py_cea.cdata.prod[k].decode("utf-8").strip() # finds species name from CEA output

        gas.mol_weight = py_cea.therm.mw[k - 1]  # kg/mol - Fortran 1-indexed
        mass_fractions = np.zeros_like(tca.x)
        for i in range(4):
            mole_fraction = py_cea.comp.en[k - 1, i]
            mass_fractions[tca.cea_station_indices[i]] = float(mole_fraction * gas.mol_weight)  # mass fraction of species k at station i

        if sum(mass_fractions) > 0: # add to lists if species is present in the mixture
            gas.products.append(raw_name)
            gas.mass_fractions.append(interp_zeros(mass_fractions))

    gas.R = 8.3145 / gas.mol_weight  # J/kg-K - specific gas constant of combustion gases
    gas.velocity = gas.M * np.sqrt(gas.gamma * gas.R * gas.T_free)  # m/s - velocity of combustion gases at each station
    
    gas.mass_fractions = np.array(gas.mass_fractions)  # shape: (n_species, 4)
    
    for i in range(len(tca.x)):
        gas.cp[i], gas.visc[i], gas.cond[i], gas.dens[i] = mixture(
            gas.products, 
            gas.mass_fractions[:, i], 
            gas.T_boundary[i], 
            gas.p[i]
        )
        gas.prandtl = np.divide(gas.cp * gas.visc, gas.cond, out=np.zeros_like(gas.cp), where=gas.cond > 0)  # Prandtl Number

    # Bartz relation
    sigma = 1 / (
        (0.5 * (jacket.T_hot / tca.Tt) * (1 + (gas.gamma - 1) / 2 * gas.M**2) + 0.5) ** 0.68
        * (1 + (gas.gamma - 1) / 2 * gas.M**2) ** 0.12
    )

    gas.h = (
        (0.026 / tca.d_throat ** 0.2)
        * (gas.visc ** 0.2 * gas.cp / gas.prandtl ** 0.6)
        * (gas.p / tca.c_star) ** 0.8
        * (tca.d_throat / tca.rc_throat) ** 0.1
        * (tca.r_throat / tca.r1) ** (2*0.9)
        * sigma
    )
    # K - Adiabatic Wall Temperature. 0.5 for laminar flow, 0.33 for turbulemt
    r = gas.prandtl ** 0.33  # recovery factor
    gas.T_adiabatic_wall = gas.T_free * (1 + r * (gas.gamma - 1) / 2 * gas.M ** 2)


    #######################################################################################################
    ### Film Cooling
    #######################################################################################################

    ## Film Parameters
    film.mdot = film.fraction * tca.mdot
    gas.mdot = (1 - np.sum(film.fraction)) * tca.mdot
    T_injection = engine.T_amb  # K - PLACEHOLDER UNTIL REPLACED WITH REGEN TEMPERATURE AT THAT LOCATION
    film.T = T_injection * np.ones_like(tca.x)  # K - film coolant temperature distribution
    # film.phase = np.zeros_like(tca.x) # Film phase map: 0 means film coolant is a gas, and 1 means it is a liquid at that location
    film.T_sat = np.zeros_like(tca.x)  # K - saturation temperature of film coolant at each injection location
    
    # Precompute film coolant properties at every station
    film.cp = np.array([4180 * rprop_fuel.CpAtTdegR(T * 1.8) for T in film.T]) # J/kg-K
    film.cond = np.array([0.144 * rprop_fuel.CondAtTdegR(T * 1.8) for T in film.T]) # W/m-K
    film.heat_vap = np.array([2326*rprop_fuel.HvapAtTdegR(T * 1.8) for T in film.T]) # J/kg
    film.T_sat = np.array([rprop_fuel.TdegRAtPsat(p / 6894.76) / 1.8 for p in gas.p]) # K

    # film.velocity = gas.velocity # m/s - BOLD ASSUMPTION
    drdx = np.gradient(tca.r1, tca.x) # unitless - wall slope
    dA = 2 * np.pi * tca.r1 * np.sqrt(1 + drdx**2) * np.gradient(tca.x) # m2 - station area

#############################

    film.liquid_mdot = np.zeros_like(tca.x)
    film.liquid_mdot[np.round(np.divide(film.injection_x, tca.dx)).astype(int)] = film.mdot
    film.gas_mdot = np.zeros_like(tca.x)
    film.phase = np.zeros(len(tca.x), dtype=int)  # 0 = no film, 1 = liquid, 2 = gas

    for i in range(1, len(tca.x)):

        # Carry forward from previous station, mixing temperatures with any new injection at this station
        film.liquid_mdot[i] = film.liquid_mdot[i-1]
        film.gas_mdot[i] = film.gas_mdot[i-1]
        # if total_liquid > 0:
        #     film.T[i] = (film.liquid_mdot[i-1] * film.T[i-1] + new_inject_liquid * film.T[i]) / total_liquid
        # elif film.gas_mdot[i-1] > 0:
        #     film.T[i] = film.T[i-1]

        # Heat absorbed by film from combustion gas
        q = gas.h[i] * dA[i] * (gas.T_adiabatic_wall[i] - film.T[i]) # W

        # Raise liquid film temperature; vaporize if saturation temperature is reached
        if film.liquid_mdot[i] > 0:
            T_new = film.T[i] + q / (film.liquid_mdot[i] * film.cp[i])
            if T_new >= film.T_sat[i]:
                q_sensible = film.liquid_mdot[i] * film.cp[i] * (film.T_sat[i] - film.T[i])
                dm_vap = min(max(q - q_sensible, 0) / film.heat_vap[i], film.liquid_mdot[i])
                film.liquid_mdot[i] -= dm_vap
                film.gas_mdot[i]    += dm_vap
                film.T[i] = film.T_sat[i]
            else:
                film.T[i] = T_new

        # Raise gas film temperature (only when no liquid film remains at this station)
        elif film.gas_mdot[i] > 0:
            film.T[i] += q / (film.gas_mdot[i] * film.cp[i])

        # Set phase flag for downstream wall heat balance
        if film.liquid_mdot[i] > 0:
            film.phase[i] = 1
        elif film.gas_mdot[i] > 0:
            film.phase[i] = 2


#############################

    # for i in range(len(film.injection_x)):
    #     liquid_start = int(np.floor(film.injection_x[i] / tca.dx))
    #     # March the liquid film temperature downstream
    #     film.T[liquid_start] = T_injection

    #     vaporization_start = len(tca.x) - 1  # default if never saturates
        
    #     for j in range(liquid_start + 1, len(tca.x)):
    #         # Heat transfer to liquid film (convection from gas to film surface)
    #         q_to_film = gas.h[j] * dA[j] * (gas.T_adiabatic_wall[j] - film.T[j-1])
            
    #         # Temperature rise of liquid film
    #         dT_film = q_to_film / (film.mdot[i] * film.cp[j-1])
    #         film.T[j] = film.T[j-1] + dT_film
            
    #         # Check if film reaches saturation
    #         if film.T[j] >= film.T_sat[j]:
    #             film.T[j] = film.T_sat[j]
    #             vaporization_start = j # start boiling from this station onward
    #             break
        
    #     dt = 1
    #     boiling_energy = film.injection_efficiency * film.mdot[i] * film.heat_vap[vaporization_start] * dt # J

    #     accumulated_energy = np.cumsum(gas.h[vaporization_start:] * (gas.T_adiabatic_wall[vaporization_start:]-film.T_sat[vaporization_start:]) * dA[vaporization_start:]) * dt # J

    #     liquid_end_indices = np.where(accumulated_energy >= boiling_energy)[0]
    #     if len(liquid_end_indices) == 0:
    #         liquid_end = len(tca.x) - 1
    #     else:
    #         liquid_end = liquid_end_indices[0] + vaporization_start

    #     # Apply liquid film region
    #     film.phase[liquid_start:liquid_end + 1] = 1
    #     film.T[vaporization_start:liquid_end + 1] = film.T_sat[vaporization_start]

    #     # # March gaseous film downstream
    #     # for j in range(liquid_end+1, len(tca.x)):
    #     #     dT = (tca.dx * 1.628 * (2 * np.pi * tca.r1[j]) * gas.h[j] / film.cp[j]
    #     #         * (gas.velocity[j] / film.v_injection) * (1 / (2 * np.pi * tca.r1[j]))
    #     #         * (film.cp[j] / film.cond[j]) ** 0.125 * film.mdot[i] ** -0.875
    #     #         * (T_adiabatic_wall[j] - film.T[j - 1]))
    #     #     film.T[j] = film.T[j - 1] + dT

    #######################################################################################################
    ### Regenerative Cooling
    #######################################################################################################

    # March coolant (regen) along flow direction
    if cool.flow_direction == 1:
        march_indices = range(0, len(tca.x))
    else:
        march_indices = range(len(tca.x) - 1, -1, -1)

    # Initialize coolant inlet temperature
    T_input = engine.T_amb # K - PLACEHOLDER UNTIL REPLACED WITH FUEL MAINLINE TEMPERATURE

    cool.cp = np.zeros_like(tca.x)
    cool.visc = np.zeros_like(tca.x)
    cool.cond = np.zeros_like(tca.x)
    cool.prandtl = np.zeros_like(tca.x)
    cool.density = np.zeros_like(tca.x)
    cool.nusselt = np.zeros_like(tca.x)
    cool.A = np.zeros_like(tca.x)
    jacket.q = np.zeros_like(tca.x)

    for i in march_indices:
        # Regenerative Coolant Properties (evaluate at current temperature)
        cool.cp[i] = 4186.8*rprop_fuel.CpAtTdegR(cool.T[i] * 1.8)  # J/kg-K
        cool.visc[i] = 0.1*rprop_fuel.ViscAtTdegR((cool.T[i] * 1.8))     # Pa-s
        cool.cond[i] = 1.7307*rprop_fuel.CondAtTdegR(cool.T[i] * 1.8)  # W/m-K
        cool.prandtl[i] = cool.cp[i] * cool.visc[i] / cool.cond[i]  # unitless - Prandtl number
        cool.density[i] = 1000 * rprop_fuel.SGLiqAtTdegR(cool.T[i] * 1.8) # kg/m3

        cool.velocity[i] = cool.mdot / (cool.density[i] * jacket.n_pipe[i] * jacket.w_pipe[i] * jacket.h_pipe[i])
        cool.Re[i] = cool.density[i] * cool.velocity[i] * jacket.d_hydraulic[i] / cool.visc[i]
        cool.f_darcy[i] = (0.79 * np.log(cool.Re[i]) - 1.64) ** (-2.0)

        # Nusselt number correlations (Gnielinski / laminar interpolation)
        if cool.Re[i] < 3000:  # Really should be < 2300, but using laminar flow in transition region for conservatism
            lam_x = np.array([1, 1.43, 2, 3, 4, 8, 1e10])
            lam_y = np.array([3.61, 3.73, 4.12, 4.79, 5.33, 6.49, 8.23])
            # approximate laminar Nu in rectangular duct by interpolation on h_pipe/w_pipe
            ratio = np.maximum(jacket.h_pipe[i] / np.maximum(jacket.w_pipe[i], 1e-9), 1e-6)
            cool.nusselt[i] = np.interp(ratio, lam_x, lam_y)
            if not warned:
                print(f"Warning: regen coolant goes laminar at {tca.x[i]} m from injector ({tca.x_exit-tca.x[i]} m from exit plane)\n")
                warned = True
        elif (cool.prandtl[i] >= 0.7 and cool.prandtl[i] <= 2000 and cool.Re[i] >= 3000 and cool.Re[i] <= 5e6):
            cool.nusselt[i] = (cool.f_darcy[i]/8.0) * (cool.Re[i]-1000.0) * cool.prandtl[i] / (1+12.7*(cool.f_darcy[i]/8.0)**0.5 * (cool.prandtl[i]**(2/3)-1))
        else:
            raise ValueError("No correlation for coolant Reynold's and/or Prandtl number")

        # Convective heat transfer coefficient
        cool.h[i] = cool.nusselt[i] * cool.cond[i] / jacket.d_hydraulic[i]  # W/m2-K

        # Fin efficiency
        m_coeff = np.sqrt(2 * cool.h[i] / jacket.cond * fin_thickness[i])  # coefficient for fin efficiency equation
        fin_eff = np.tanh(m_coeff * (jacket.h_pipe[i] + fin_thickness[i]) / 2) / (m_coeff * (jacket.h_pipe[i] + fin_thickness[i]) / 2)
        cool.A[i] = (2 * fin_height[i] * fin_eff + jacket.w_pipe[i]) * tca.dx  # m2 - adjusted contact area of coolant on channel walls

        # Heat Balancing
        if film.phase[i] == 1:    # Liquid State
            jacket.q[i] = (film.T[i] - cool.T[i]) / ( 1/(cool.h[i]*cool.A[i]) + 1/(jacket.cond*gas.A[i]/t_wall[i]) )  # W
        elif film.phase[i] == 2:  # Gaseous State
            jacket.q[i] = (film.T[i] - cool.T[i]) / ( 1/(cool.h[i]*cool.A[i]) + 1/(jacket.cond*gas.A[i]/t_wall[i]) + 1/(gas.h[i]*gas.A[i]) )  # W - NEED TO CHANGE H TO FILM CONVECTION COEFFICIENT ON GAS SIDE
        else:                     # No Film
            jacket.q[i] = (gas.T_adiabatic_wall[i] - cool.T[i]) / ( 1/(cool.h[i]*cool.A[i]) + 1/(jacket.cond*gas.A[i]/t_wall[i]) + 1/(gas.h[i]*gas.A[i]) )  # W

        # Split wall temperatures (hot/cold) across the wall thickness with 1-D conduction
        jacket.T_cold[i] = jacket.q[i] / (cool.h[i] * cool.A[i]) + cool.T[i]  # K
        jacket.T_hot[i] = jacket.q[i] / (jacket.cond * gas.A[i] / t_wall[i]) + jacket.T_cold[i]  # K

        # March coolant temperature
        dT = jacket.q[i] / cool.cp[i] / (tca.mdot_fuel/jacket.n_pipe[i]) # Coolant Temperature change at one channel BLEED MDOT AT FILM INJECTION
        nxt = i + cool.flow_direction
        if 0 <= nxt < len(tca.x):
            cool.T[nxt] = cool.T[i] + dT

        # SHOULD ADD VAPOR MARGIN
        # Check if the regen coolant boils. Gaseous coolant will be extremely ineffective, as well as harder to analyze
        boiling_cool = rprop_fuel.TdegRAtPsat(cool.p[i]/6894.76)/1.8 # K - fuel saturation temperature
        if cool.T[i] > boiling_cool:
            raise ValueError(f"Coolant starts boiling at {tca.x[i]} m from injector ({tca.x_exit-tca.x[i]} m from exit plane)")
        
def interp_zeros(arr): # replaces zeros in sparse array with linearly interpolated values
    return np.interp(np.arange(len(arr)), (nz := np.nonzero(arr)[0]), arr[nz])