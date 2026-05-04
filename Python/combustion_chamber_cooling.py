from os import error
import numpy as np
import rocketcea.py_cea as py_cea

from engine_state import engine, tca
from mixture import mixture
from rocketcea.cea_obj import CEA_Obj as rcea
from rocketprops.rocket_prop import get_prop as rprop

from fuel import n_dodecane

def combustion_chamber_cooling():

    gas = tca.combustion_gas
    jacket = tca.regenerative_jacket
    cool = tca.regenerative_coolant
    film = tca.film_coolant

    # rprop_fuel = rprop(engine.fuelrp)
    fuel = n_dodecane()

    #######################################################################################################
    ### Film Coolant Precomputing
    #######################################################################################################
    film.liquid.injection_mdot = np.zeros_like(tca.x) # kg/s - mass flow rate of film coolant injected at each station
    film.liquid.injection_mdot[np.round(np.divide(film.injection_x, tca.dx)).astype(int)] = film.fraction*tca.mdot_fuel # [kg/s]
    
    
    film.T_critical = fuel.Tc  # K - critical temperature of film
    film.p_critical = fuel.Pc * 6894.76  # Pa - critical pressure of film

    film.mw = fuel.M  # g/mol - molecular weight of RP-1
    film.liquid.velocity = np.zeros_like(tca.x) # m/s - velocity of liquid film coolant at each station

    film.T_sat = np.zeros_like(tca.x) # K - saturation temperature of film coolant at each station
    for attr in ['velocity', 'mdot', 'cp', 'cond', 'heat_vap', 'density', 'resistance', 'thickness']:
        setattr(film.liquid, attr, np.zeros_like(tca.x)) 
    film.liquid.T = engine.T_amb*np.ones_like(tca.x) # K - initialize liquid film temperature to ambient  
    for attr in ['velocity', 'mdot', 'cp', 'cond', 'p', 'p_sat', 'density', 'resistance', 'vaporized_mass', 'thickness', 'T']:
        setattr(film.gas, attr, np.zeros_like(tca.x))
    film.gas.T = engine.T_amb*np.ones_like(tca.x) # K - initialize gas film temperature to ambient

    #######################################################################################################
    ### Regenerative Coolant Precomputing
    #######################################################################################################
    for attr in ['dp1', 'dp2', 'dp3', 'dp', 'cp', 'cond', 'prandtl', 'nusselt', 'h', 'viscosity', 'density', 'velocity', 'Re', 'f_darcy', 'A', 'dyn_viscosity', 'kin_viscosity', 'mdot', 'resistance']:
        setattr(cool, attr, np.zeros_like(tca.x))
    for attr in ['T_cold', 'T_hot', 'resistance']:
        setattr(jacket, attr, np.zeros_like(tca.x))
    
    jacket.T_hot = 0.5 * tca.Tt * np.ones(len(tca.x)) # K - Initial guess for the hot wall temperatures
    cool.p = np.full(len(tca.x), float(engine.p_amb))     # Pa - coolant pressure (initialized at p_amb)
    jacket.d_hydraulic = 4.0 * jacket.w_pipe * jacket.h_pipe / (2.0 * jacket.w_pipe + 2.0 * jacket.h_pipe)  # m - hydraulic diameter
    jacket.pipe_area = jacket.w_pipe * jacket.h_pipe # m2 - cross-sectional area of one channel
    jacket.h_pipe = jacket.h_pipe*np.ones_like(tca.x)
    
    cool.mdot = np.full_like(tca.x, tca.mdot_fuel) # kg/s - decreases as film is injected
    film_injected_cumulative = np.cumsum(film.liquid.injection_mdot)

    # March coolant (regen) along flow direction
    if cool.flow_direction == 1:
        flow_entrance_index = 0
        flow_exit_index = len(tca.x) - 1
        cool.mdot = tca.mdot_fuel - film_injected_cumulative
    else:
        flow_entrance_index = len(tca.x) - 1
        flow_exit_index = 0
        cool.mdot = tca.mdot_fuel - (film_injected_cumulative[-1] - film_injected_cumulative)
    next = cool.flow_direction  # index offset for next station
    cool.T = 400*np.ones_like(tca.x) # K - presets coolant temperatures to get reasonable properties on a first pass

    #######################################################################################################
    ### Combustion Gas Precomputing
    #######################################################################################################
    
    for attr in ['q1', 'q2', 'q3']:
        setattr(tca, attr, np.zeros_like(tca.x))
    for attr in ['resistance']:
        setattr(gas, attr, np.zeros_like(tca.x))

    # Combustion Thermal Environment
    gas.T_free = tca.Tt / (1 + (gas.gamma - 1) / 2 * gas.M**2)  # Free-Stream Temperature

    # Wall Conduction
    fin_thickness = jacket.gap_pipe*np.ones_like(tca.x)  # m
    fin_height = jacket.h_pipe     # m
    fin_length = tca.x[-1]      # m

    gas.A = 2 * np.pi * tca.r1 * tca.dx  # m2 - Area of tca.dx
    t_wall = tca.r2 - tca.r1                          # m - Wall thickness (projected vertically)

    for attr in ['T_cold', 'T', 'T_adiabatic_wall', 'cp', 'viscosity', 'cond', 'h', 'q', 'prandtl', 'density', 'nusselt']:
        setattr(gas, attr, np.zeros_like(tca.x))
    warned = False # warning flag for laminar coolant flow

    gas.T_boundary = gas.T_free * (1 + 0.032 * gas.M**2 + 0.58 * (jacket.T_hot / gas.T_free - 1)) # Reference temperature for gas transport properties

    gas.mdot = (1 - np.sum(film.fraction)) * tca.mdot_fuel + tca.mdot_ox
    drdx = np.gradient(tca.r1, tca.x) # unitless - wall slope
    dA = 2 * np.pi * tca.r1 * np.sqrt(1 + drdx**2) * np.gradient(tca.x) # m2 - station area

    #######################################################################################################
    ### Gas Mixture Precomputing
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
        gas.cp[i], gas.viscosity[i], gas.cond[i], gas.density[i] = mixture(
            gas.products, 
            gas.mass_fractions[:, i], 
            gas.T_boundary[i], 
            gas.p[i]
        )
        gas.prandtl = np.divide(gas.cp * gas.viscosity, gas.cond, out=np.zeros_like(gas.cp), where=gas.cond > 0)  # Prandtl Number

    # Bartz relation
    sigma = 1 / (
        (0.5 * (jacket.T_hot / tca.Tt) * (1 + (gas.gamma - 1) / 2 * gas.M**2) + 0.5) ** 0.68
        * (1 + (gas.gamma - 1) / 2 * gas.M**2) ** 0.12
    )

    gas.h = (
        (0.026 / tca.d_throat ** 0.2)
        * (gas.viscosity ** 0.2 * gas.cp / gas.prandtl ** 0.6)
        * (gas.p / tca.c_star) ** 0.8
        * (tca.d_throat / tca.rc_throat) ** 0.1
        * (tca.r_throat / tca.r1) ** (2*0.9)
        * sigma
    )
    
    r = gas.prandtl ** 0.33  # recovery factor - 0.5 for laminar flow, 0.33 for turbulent
    gas.T_adiabatic_wall = gas.T_free * (1 + r * (gas.gamma - 1) / 2 * gas.M ** 2) # K - Adiabatic Wall Temperature. 

    #######################################################################################################
    ### Main Loop
    #######################################################################################################
    for j in range(2):
        for i in range(len(tca.x)):

            #######################################################################################################
            ### Regenerative Coolant
            #######################################################################################################

            cool.cp[i]            = fuel.cp(cool.T[i], cool.p[i])
            cool.cond[i]          = fuel.k(cool.T[i], cool.p[i])
            cool.density[i]       = fuel.rho(cool.T[i], cool.p[i])
            cool.dyn_viscosity[i] = fuel.mu(cool.T[i], cool.p[i])
            cool.kin_viscosity[i] = cool.dyn_viscosity[i] / cool.density[i] # m2/s - kinematic viscosity of RP1
            cool.prandtl[i] = cool.cp[i] * cool.dyn_viscosity[i] / cool.cond[i]  # unitless - Prandtl number
            cool.velocity[i] = cool.mdot[i] / (cool.density[i] * jacket.n_pipe[i] * jacket.w_pipe[i] * jacket.h_pipe[i]) # m/s
            cool.Re[i] = cool.density[i] * cool.velocity[i] * jacket.d_hydraulic[i] / cool.dyn_viscosity[i] # unitless
            cool.f_darcy[i] = (0.79 * np.log(cool.Re[i]) - 1.64) ** (-2.0) # unitless

            # Nusselt number correlations (Gnielinski / laminar interpolation)
            if cool.Re[i] < 3000:  # Really should be < 2300, but using laminar flow in transition region for conservatism
                lam_x = np.array([1, 1.43, 2, 3, 4, 8, 1e10]) 
                lam_y = np.array([3.61, 3.73, 4.12, 4.79, 5.33, 6.49, 8.23])
                # approximate laminar Nu in rectangular duct by interpolation on h_pipe/w_pipe
                cool.nusselt[i] = np.interp(jacket.h_pipe[i] / jacket.w_pipe[i], lam_x, lam_y) # PROBABLY SHOULD NOT BE A LINEAR INTERPOLATION
                if not warned:
                    print(f"Warning: regen coolant goes laminar at {tca.x[i]} m from injector ({tca.x_exit-tca.x[i]} m from exit plane)\n")
                    warned = True
            elif (cool.prandtl[i] >= 0.7 and cool.prandtl[i] <= 2000 and cool.Re[i] >= 3000 and cool.Re[i] <= 5e6): # turbulent case
                cool.nusselt[i] = (cool.f_darcy[i]/8.0) * (cool.Re[i]-1000.0) * cool.prandtl[i] / (1+12.7*(cool.f_darcy[i]/8.0)**0.5 * (cool.prandtl[i]**(2/3)-1))
            else:
                raise ValueError(f"No correlation for coolant Re = {cool.Re[i]} and/or Pr = {cool.prandtl[i]} number\ncool.T = {cool.T[i]} K, cool.p = {cool.p[i]} Pa, cool.density = {cool.density[i]} kg/m3, cool.viscosity = {cool.dyn_viscosity[i]} Pa-s\n")                               

            # Convective heat transfer coefficient
            cool.h[i] = cool.nusselt[i] * cool.cond[i] / jacket.d_hydraulic[i]  # W/m2-K

            # Fin efficiency
            m_coeff = np.sqrt(2 * cool.h[i] / jacket.cond * fin_thickness[i])  # coefficient for fin efficiency equation
            fin_eff = np.tanh(m_coeff * (jacket.h_pipe[i] + fin_thickness[i]) / 2) / (m_coeff * (jacket.h_pipe[i] + fin_thickness[i]) / 2)
            cool.A[i] = (2 * fin_height[i] * fin_eff + jacket.w_pipe[i]) * tca.dx * jacket.n_pipe[i]  # m2 - adjusted contact area of coolant on channel walls

            # Major Losses (channel wall roughness)
            cool.dp1[i] = - 0.5 * cool.f_darcy[i] * cool.density[i] * cool.velocity[i]**2 * tca.dx / jacket.d_hydraulic[i] # Pa - major loss

            # Minor Losses (diameter / channel-count changes)
            if jacket.n_pipe[i+next] == jacket.n_pipe[i]:  # most of the time: same number of channels
                d_hydraulic_mod = jacket.d_hydraulic[i+next]
                A1 = jacket.w_pipe[i] * jacket.h_pipe[i] # channel area before
                A2 = jacket.w_pipe[i+next] * jacket.h_pipe[i+next] # channel area after
            elif jacket.n_pipe[i+next] < jacket.n_pipe[i]:  # channels merge
                d_hydraulic_mod = 4.0 * 0.5 * jacket.w_pipe[i+next] * jacket.h_pipe[i+next] / (2.0 * 0.5 * jacket.w_pipe[i+next] + 2.0 * jacket.h_pipe[i+next])  # half-channel
                A1 = jacket.w_pipe[i] * jacket.h_pipe[i] # channel area before merge
                A2 = 0.5 * jacket.w_pipe[i+next] * jacket.h_pipe[i+next] # channel area after merge
            elif jacket.n_pipe[i+next] > jacket.n_pipe[i]:  # channels divide
                d_hydraulic_mod = 4.0 * 0.5 * jacket.w_pipe[i] * jacket.h_pipe[i] / (2.0 * 0.5 * jacket.w_pipe[i] + 2.0 * jacket.h_pipe[i])  # half-channel
                A1 = 0.5 * jacket.w_pipe[i] * jacket.h_pipe[i] # channel area before divide
                A2 = jacket.w_pipe[i+next] * jacket.h_pipe[i+next] # channel area after divide

            # Contraction / expansion coefficient (empirical)
            d_ratio = d_hydraulic_mod / jacket.d_hydraulic[i]
            if d_ratio < 1.0: # Contraction
                K = 0.5 - 0.167 * d_ratio - 0.125 * d_ratio ** 2 - 0.208 * d_ratio ** 3 
            else: # Expansion
                K = (d_ratio ** -2.0 - 1.0) ** 2.0
            cool.dp2[i] = - 0.5 * K * cool.density[i] * cool.velocity[i] ** 2  # Pa - Minor Losses (diameter changes)
            cool.dp3[i] = - 0.5*cool.density[i]*(cool.velocity[i+next]**2 - cool.velocity[i]**2) # Pa - Acceleration Losses

            cool.dp[i] = cool.dp1[i] + cool.dp2[i] + cool.dp3[i]  # Pa - Total Losses

            # March pressure to the next station
            if i == 0:
                cool.p[i] = tca.p_injector_manifold
            else:
                cool.p[i] = cool.p[i-1] - 1*cool.dp[i-1]

            cool.velocity[i] = cool.mdot[i] / (cool.density[i] * jacket.n_pipe[i] * jacket.pipe_area[i]) # m/s - fluid bulk speed at each station
            cool.Re[i] = cool.velocity[i] * jacket.d_hydraulic[i] / cool.kin_viscosity[i] # unitless - Reynolds number
            cool.f_darcy[i] = (0.79 * np.log(cool.Re[i]) - 1.64) ** (-2.0) # unitless - friction factor, smooth pipe approximation REVIEW IF THIS IS ACCURATE

            #######################################################################################################
            ### Film Coolant
            #######################################################################################################

            film.liquid.velocity[i] = gas.velocity[i] # MODIFY WITH ENTRAINMENT MODEL
            dt = tca.dx / film.liquid.velocity[i] # s - time for gas to travel one station

            film.liquid.resistance[i] = 0
            film.gas.resistance[i] = 0

            if i == 0:
                i_prev = 0
            else:
                i_prev = i-1
            
            film.liquid.mdot[i] = film.liquid.mdot[i_prev] + film.liquid.injection_mdot[i] - film.gas.vaporized_mass[i_prev]  # kg/s
            film.gas.mdot[i] = film.gas.mdot[i_prev] + film.gas.vaporized_mass[i_prev] # kg/s

            if film.liquid.injection_mdot[i] > 0:  # Injection point
                if film.liquid.mdot[i_prev] > 0:
                    film.liquid.T[i] = (film.liquid.mdot[i_prev] * film.liquid.T[i] + film.liquid.injection_mdot[i] * cool.T[i]) / (film.liquid.mdot[i_prev] + film.liquid.injection_mdot[i])
                else:
                    film.liquid.T[i] = cool.T[i]  # First injection, no upstream

            # film T_sat at injection
            film.T_sat[i] = fuel.Tsat(gas.p[i]) if not fuel.is_supercritical(gas.p[i]) else fuel.Tc

            if film.liquid.mdot[i] > 0:
                # film liquid properties
                film.liquid.cp[i]       = fuel.cp(film.liquid.T[i], gas.p[i])
                film.liquid.cond[i]     = fuel.k(film.liquid.T[i], gas.p[i])
                film.liquid.density[i]  = fuel.rho(film.liquid.T[i], gas.p[i])
                film.liquid.heat_vap[i] = fuel.hvap(film.liquid.T[i])  # 0 above Tc, handled

            if film.gas.mdot[i] > 0:
                # film gas properties
                film.gas.cp[i]      = fuel.cp(film.gas.T[i], gas.p[i])
                film.gas.cond[i]    = fuel.k(film.gas.T[i], gas.p[i])
                film.gas.density[i] = fuel.rho(film.gas.T[i], gas.p[i])

            film.liquid.thickness[i] = film.liquid.mdot[i] / (film.liquid.density[i]  * gas.velocity[i] * 2*np.pi*tca.r1[i] )  # m VELOCITY PROBABLY WRONG
            film.gas.thickness[i] = film.gas.mdot[i]  / (film.gas.density[i] * gas.velocity[i] * 2*np.pi*tca.r1[i] )  # m
            film.gas.resistance[i] = film.gas.thickness[i] / (film.gas.cond[i] * gas.A[i])              # K/W - gaseous RP-1 film conduction
            film.liquid.resistance[i] = film.liquid.thickness[i] / (film.liquid.cond[i] * gas.A[i]) # K/W

            # Thermal resistances
            gas.resistance[i] = 1 / (gas.h[i] * gas.A[i])                      # K/W - combustion gas boundary layer convection
            jacket.resistance[i] = t_wall[i] / (jacket.cond * gas.A[i])            # K/W - wall conduction
            cool.resistance[i] = 1 / (cool.h[i] * cool.A[i])                    # K/W - regen coolant convection

            # Heat balance - q's between the resistive layers. T's are bulk temps of the layers
            if film.gas.mdot[i] > 0 and film.liquid.mdot[i] > 0: # two-phase film coolant
                tca.q1[i] = (gas.T_adiabatic_wall[i] - film.gas.T[i]) / (gas.resistance[i]+film.gas.resistance[i]/2)  # W
                tca.q2[i] = (film.gas.T[i] - film.liquid.T[i]) / (film.gas.resistance[i]/2+film.liquid.resistance[i]/2)  # W
                tca.q3[i] = (film.liquid.T[i] - cool.T[i]) / (film.liquid.resistance[i]/2 + jacket.resistance[i]/2 + cool.resistance[i]/2)  # W
            elif film.gas.mdot[i] > 0 and film.liquid.mdot[i] == 0: # fully vaporized film coolant
                tca.q1[i] = (gas.T_adiabatic_wall[i] - film.gas.T[i]) / (gas.resistance[i]+film.gas.resistance[i]/2)  # W
                tca.q2[i] = (film.gas.T[i] - cool.T[i]) / (film.gas.resistance[i]/2 + jacket.resistance[i]/2 + cool.resistance[i]/2)  # W
                tca.q3[i] = tca.q2[i]
            elif film.gas.mdot[i] == 0 and film.liquid.mdot[i] > 0: # fully liquid film coolant
                tca.q1[i] = (gas.T_adiabatic_wall[i] - film.liquid.T[i]) / (gas.resistance[i]+film.liquid.resistance[i]/2)  # W
                tca.q2[i] = tca.q1[i]
                tca.q3[i] = (film.liquid.T[i] - cool.T[i]) / (film.liquid.resistance[i]/2 + jacket.resistance[i]/2 + cool.resistance[i]/2)  # W
            else:
                tca.q1[i] = (gas.T_adiabatic_wall[i] - cool.T[i]) / (gas.resistance[i] + jacket.resistance[i] + cool.resistance[i]/2)  # W
                tca.q2[i] = tca.q1[i]
                tca.q3[i] = tca.q1[i]

            jacket.T_cold[i] = cool.T[i] + tca.q3[i] * cool.resistance[i]/2  # K - wall temperature on coolant side
            jacket.T_hot[i] = jacket.T_cold[i] + tca.q3[i] * jacket.resistance[i]  # K - wall temperature on hot side

            if i == flow_entrance_index:
                cool.T[i] = 400 # K - CHANGE TO PUMP DISCHARGE TEMPERATURE
            else:
                cool.T[i] = cool.T[i-1] + next * tca.q3[i] / (cool.mdot[i] * cool.cp[i]) # K

            # March film temperatures only when the phase exists
            if film.gas.mdot[i] > 0 and film.gas.cp[i] > 0:
                film.gas.T[i] = film.gas.T[i-1] + (tca.q1[i]-tca.q2[i])/(film.gas.mdot[i]*film.gas.cp[i])
            else:
                film.gas.T[i] = film.liquid.T[i] if film.liquid.mdot[i] > 0 else cool.T[i]
                # Set to something physically sensible so next station's property call doesn't NaN

            if film.liquid.mdot[i] > 0 and film.liquid.cp[i] > 0:
                film.liquid.T[i+1] = film.liquid.T[i] + (tca.q2[i]-tca.q3[i])/(film.liquid.mdot[i]*film.liquid.cp[i])
            # else: leave as initialized (engine.T_amb from precompute)  

            
            if film.liquid.T[i] > film.T_sat[i]:
                if not fuel.is_supercritical(gas.p[i]):
                    # Subcritical boiling: latent heat sink, T held at Tsat
                    film.gas.vaporized_mass[i] = (
                        (film.liquid.T[i] - film.T_sat[i]) * film.liquid.cp[i]
                        * film.liquid.mdot[i] * dt / film.liquid.heat_vap[i]
                    )
                    film.liquid.T[i] = film.T_sat[i]
                # Supercritical: no phase change, T marches freely.
            
            # Check if the regen coolant boils. Gaseous coolant will be extremely ineffective.
            vapor_margin = 10 # K - vapor margin (only meaningful subcritically)
            if not fuel.is_supercritical(cool.p[i]):
                Tsat = fuel.Tsat(cool.p[i])
                if cool.T[i] > Tsat - vapor_margin:
                    raise ValueError(...)
                
def interp_zeros(arr):
    """Replace zeros in a sparse array with linearly interpolated values."""
    indices = np.arange(len(arr))
    nonzero = np.nonzero(arr)[0]
    return np.interp(indices, nonzero, arr[nonzero])