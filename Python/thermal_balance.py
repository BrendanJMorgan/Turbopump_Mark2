## Objective: finds the temperature distributions of the combustion chamber and nozzle walls

from os import error
import numpy as np

from engine_state import inputs, tca
from rocketcea.cea_obj import CEA_Obj as rcea
from rocketprops.rocket_prop import get_prop as rprop

def thermal_balance():

    # Combustion Thermal Environment
    T_free = tca.Tc / (1 + (tca.gamma - 1) / 2 * tca.M**2)  # Free-Stream Temperature

    # Wall Conduction
    fin_thickness = tca.w_pipe    # m
    fin_height = tca.h_pipe       # m
    fin_length = tca.x[-1]      # m

    A_gas = 2 * np.pi * tca.r1 * tca.dx / tca.n_pipe   # m2 - Area of tca.dx per coolant channel
    t_wall = tca.r2 - tca.r1                          # m - Wall thickness (projected vertically)

    # Thermal Balance (preallocations)
    products = inputs.products
    fractions = inputs.fractions_at(tca.x)  # shape: (n_species, len(x))

    T_wall_cold=np.zeros_like(tca.x)            # K - wall temperature on coolant side
    Tab =       np.zeros_like(tca.x)            # K
    T_cool =    tca.T_amb*np.ones_like(tca.x) # K
    T_ab =      tca.T_amb*np.ones_like(tca.x) # K
    cp_gas =    np.zeros_like(tca.x)            # J/kg-K - specific heat capacity (at constant pressure) of combustion gases
    visc_gas =  np.zeros_like(tca.x)            # Pa-s - dynamic viscosity of combustion gases
    cond_gas =  np.zeros_like(tca.x)            # W/m-K - thermal conductivity of combustion gases
    h_gas =     np.zeros_like(tca.x)            # W / m2-K - convective heat transfer coefficient
    q_gas =     np.zeros_like(tca.x)
    Pr_gas =    np.zeros_like(tca.x)
    Nu_cool =   np.zeros_like(tca.x)
    h_cool =    np.zeros_like(tca.x)

    warned = False # warning flag for laminar coolant flow

    # Reference temperature for gas transport properties
    T_ref = T_free * (1 + 0.032 * tca.M**2 + 0.58 * (T_wall_cold / np.maximum(T_free, 1e-9) - 1))  # T_wall_hot not known yet start with cold side
    # Get mixture properties along x
    for i in range(len(tca.x)):
        cp_gas[i], visc_gas[i], cond_gas[i] = mixture(products, fractions[:, i], float(T_ref[i]), float(tca.p[i]))
    Pr_gas = np.divide(cp_gas * visc_gas, cond_gas, out=np.zeros_like(cp_gas), where=cond_gas > 0)  # Prandtl Number

    # Bartz relation
    sigma = 1 / ((0.5*(T_wall_cold/tca.Tc)*((1+(tca.gamma-1)/2*tca.M**2)+0.5))**0.68 * ((1+(tca.gamma-1)/2*tca.M**2)**0.12))
    h_gas = (
        (0.026 / (tca.d_throat ** 0.2))
        * (visc_gas ** 0.2)
        * (cp_gas ** 0.6)
        * ((tca.p[i] / tca.c_star_eff) ** 0.8)
        * ((tca.d_throat / tca.rc_throat) ** 0.1)
        * sigma
    )

    # K - Adiabatic Wall Temperature. 0.33 ==> 0.5 IF GAS
    T_ab = tca.Tc * (1 + Pr_gas ** 0.33 * ((tca.gamma - 1) / 2) * tca.M ** 2) / (1 + ((tca.gamma - 1) / 2) * tca.M ** 2)
    # K - recovery temperature of combustion gas (fluid state if it slowed to stagnation)
    T_recovery = tca.c_star_eff**2 * (T_free + Pr_gas ** (1/3) * (tca.Tc - T_free))





    #######################################################################################################
    ### Film Cooling
    #######################################################################################################

    ## Film Parameters
    mdot_film = inputs.film_fraction * tca.mdot_cc * np.ones(2)  # kept your [0.5,0.5]*... idea if split tune as needed
    mdot_gas = (1 - inputs.film_fraction) * tca.mdot_cc
    film_injection_x = [0, tca.l_chamber] # m - film cooling orifices around perimeter of injector and along bottom edge of chamber wall
    T_injection = tca.T_amb  # K - PLACEHOLDER UNTIL REPLACED WITH REGEN TEMPERATURE AT THAT LOCATION
    T_sat_film = inputs.rprop_Tsat_film()  # K - saturation temp of film liquid (e.tca., ethanol/water mix)
    T_film = tca.T_amb * np.ones_like(tca.x)  # K - film coolant temperature distribution
    film_phase = np.zeros_like(tca.x) # Film phase map: 0 means film coolant is a gas, and 1 means it is a liquid at that location

    # Liquid Film Coolant Distribution
    for i in range(len(film_injection_x)):
        # Liquid Film Coolant Distribution
        boiling_energy = tca.injection_efficiency * mdot_film[i] * (cp_film * (T_sat_film - T_injection) + heat_vap_film)  # J - energy required to boil all liquid film coolant into gas
        liquid_start = int(np.floor(film_injection_x[i] / tca.dx)) + 1  # location at which liquid film coolant is injected
        
        # J - total heat energy outputted by combustion gas up to each point
        accumulated_energy = np.cumsum(
            h_gas[liquid_start:] * (T_recovery[liquid_start:] - T_sat_film * np.ones(len(tca.x) - liquid_start)) * dA[liquid_start:]
        )
        
        # location at which film coolant boils into a gas
        liquid_end_indices = np.where(accumulated_energy >= boiling_energy)[0]
        if len(liquid_end_indices) == 0:
            liquid_end = len(accumulated_energy) + liquid_start - 1
        else:
            liquid_end = liquid_end_indices[0] + liquid_start - 1
        
        film_phase[liquid_start:liquid_end+1] = 1  # 1 for liquid and 0 for gas
        T_film[liquid_start:liquid_end+1] = T_sat_film  # temperature of the liquid coolant is always the saturation temperature
        
        # Gaseous Film Coolant Distribution
        for j in range(liquid_end + 1, len(tca.x)):
            # cp_film = 4180 * rocketprop('RP1').CpAtTdegR(T_film[j-1] * 1.8)  # J/kg-K
            # k_film = 0.144 * rocketprop('RP1').CondAtTdegR(T_film[j-1] * 1.8)  # W/m-K
            # cp_film = PropsSI('C', 'P', p_gas[j], 'T', T_film[j-1], 'Ethanol')  # J/kg-K
            # k_film = PropsSI('L', 'P', p_gas[j], 'T', T_film[j-1], 'Ethanol')  # W/m-K
            
            dT = (tca.dx * 1.628 * (2 * np.pi * r1[j]) * h_gas[j] / cp_film 
                * ((v_gas[j]) / v_injection) * (1 / (2 * np.pi * tca.r1[j])) 
                * (cp_film / k_film) ** 0.125 * (mdot_film[i]) ** -0.875 
                * (T_recovery[j] - T_film[j-1]))  # K
            
            T_film[j] = T_film[j-1] + dT




    #######################################################################################################
    ### Regenerative Cooling
    #######################################################################################################

    # March coolant (regen) along flow direction
    flow_dir = inputs.flow_direction  # +1 or -1
    if flow_dir == 1:
        march_indices = range(0, len(tca.x))
    else:
        march_indices = range(len(tca.x) - 1, -1, -1)

    # Initialize coolant inlet temperature
    T_input = tca.T_amb

    for i in march_indices:
        ## Regenerative Coolant Properties (evaluate at current T_cool)
        cp_cool = inputs.rprop_cp_cool(T_cool[i])        # J/kg-K
        visc_cool = inputs.rprop_mu_cool(T_cool[i])      # Pa-s
        k_cool = inputs.rprop_k_cool(T_cool[i])          # W/m-K
        Pr_cool = cp_cool * visc_cool / k_cool           # Prandtl Number

        # Nusselt number correlations (Gnielinski / laminar interpolation)
        if tca.Re_cool[i] < 3000:  # Really should be < 2300, but using laminar flow in transition region for conservatism
            lam_x = np.array([1, 1.43, 2, 3, 4, 8, 1e10, 3.61, 3.73, 4.12, 4.79, 5.33, 6.49, 8.23])
            lam_y = lam_x.copy()   # placeholder—MATLAB used two rows: provide your table here
            # approximate laminar Nu in rectangular duct by interpolation on h_pipe/w_pipe
            ratio = np.maximum(tca.h_pipe[i] / np.maximum(tca.w_pipe[i], 1e-9), 1e-6)
            Nu_cool[i] = np.interp(ratio, lam_x, lam_y)
            if not warned:
                print(f"Warning: regen coolant goes laminar at {tca.x[i]} m from injector ({tca.x_exit-tca.x[i]} m from exit plane)\n")
                warned = True
        elif (Pr_cool >= 0.7 and Pr_cool <= 2000 and tca.Re_cool[i] >= 3000 and tca.Re_cool[i] <= 5e6):
            Nu_cool[i] = (tca.f_cool/8.0) * (tca.Re_cool[i]-1000.0) * Pr_cool / (1+12.7*(tca.f_cool/8.0)**0.5 * (Pr_cool**(2/3)-1))
        else:
            raise ValueError("No correlation for coolant Reynold's and/or Prandtl number")

        # Convective heat transfer coefficient
        h_cool[i] = Nu_cool[i] * k_cool / tca.d_hydraulic[i]  # W/m2-K

        # Fin efficiency
        m_coeff = np.sqrt(2 * h_cool[i] / tca.k_wall * fin_thickness[i])  # coefficient for fin efficiency equation
        fin_eff = np.tanh(m_coeff * (tca.h_pipe[i] + fin_thickness[i]) / 2) / (m_coeff * (tca.h_pipe[i] + fin_thickness[i]) / 2)
        A_cool = (2 * fin_height[i] * fin_eff * tca.w_pipe[i]) * tca.dx[i]  # m2 - adjusted contact area of coolant on channel walls

        # Heat Balancing
        if film_phase[i] == 1:    # Liquid State
            # q_wall: conduction from film to coolant through wall and gas-side film convection
            q_wall = (T_film[i] - T_cool[i]) / (1 / (h_cool[i] * A_cool) + 1 / (tca.k_wall * A_gas[i] / max(t_wall[i], 1e-12)))  # W
        else:                      # Gaseous State
            q_wall = (T_film[i] - T_cool[i]) / (1 / (h_gas[i]*A_gas[i]) + 1 / (tca.k_wall * A_gas[i] / t_wall[i]) + 1 / (h_cool[i] * A_cool)
            )  # W

        # Split wall temperatures (hot/cold) across the wall thickness with 1-D conduction
        T_wall_cold[i] = q_wall / (h_cool[i] * A_cool) + T_cool[i]  # K
        T_wall_hot_i = q_wall / (tca.k_wall * A_gas[i] / t_wall[i]) + T_wall_cold[i]  # K

        # March coolant temperature
        dT = q_wall * tca.n_pipe[i] / cp_cool * tca.mdot_fuel_cc  # Coolant Temperature change at one channel
        nxt = i + flow_dir
        if 0 <= nxt < len(tca.x):
            T_cool[nxt] = T_cool[i] + dT

        # SHOULD ADD VAPOR MARGIN
        # Check if the regen coolant boils. Gaseous coolant will be extremely ineffective, as well as harder to analyze
        boiling_cool = rprop.TdegRAtPsat(tca.p_cool(i))/1.8 # K - fuel saturation temperature
        if T_cool(i) > boiling_cool:
            raise ValueError(f"Coolant starts boiling at {tca.x(i)} m from injector ({tca.x_exit-tca.x(i)} m from exit plane)")