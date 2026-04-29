import numpy as np

from engine_state import engine, tca
from rocketcea.cea_obj import CEA_Obj as rcea
from rocketprops.rocket_prop import get_prop as rprop

def coolant_flow(): # https://www.cryo-rocket.com/regenerative-cooling-circuit/6.1-ltmcc-cooling-channels/

    gas = tca.combustion_gas
    jacket = tca.regenerative_jacket
    cool = tca.regenerative_coolant
    film = tca.film_coolant

    rprop_fuel = rprop(engine.fuelrp)

    cool.density = 1000*rprop_fuel.SGLiqAtTdegR(engine.T_amb*1.8) # kg/m3 - density of RP1
    cool.dyn_visc = 0.1*rprop_fuel.ViscAtTdegR(engine.T_amb*1.8) # Pa-s - dynamic viscosity of RP1
    cool.kin_visc = cool.dyn_visc / cool.density # m2/s - kinematic viscosity of RP1
    cool.mdot = tca.mdot_fuel # kg/s

    cool.p = np.full(len(tca.x), float(engine.p_amb))     # Pa - coolant pressure (initialized at p_amb)
    jacket.d_hydraulic = 4.0 * jacket.w_pipe * jacket.h_pipe / (2.0 * jacket.w_pipe + 2.0 * jacket.h_pipe)  # m - hydraulic diameter

    jacket.pipe_area = jacket.w_pipe * jacket.h_pipe # m2 - cross-sectional area of one channel
    cool.velocity = cool.mdot / (cool.density * jacket.n_pipe * jacket.pipe_area) # m/s - fluid bulk speed at each station
    cool.Re = cool.velocity * jacket.d_hydraulic / cool.kin_visc # unitless - Reynolds number
    cool.f_darcy = (0.79 * np.log(cool.Re) - 1.64) ** (-2.0) # unitless - Darcy friction factor

    # marching order depending on flow direction
    if cool.flow_direction == 1:
        idx_range = range(0, len(tca.x) - 1)
        flow_exit_index = len(tca.x) - 1
    else:
        idx_range = range(len(tca.x) - 1, 0, -1)
        flow_exit_index = 0

    next = cool.flow_direction  # index offset for next station
    jacket.h_pipe = jacket.h_pipe*np.ones_like(tca.x)

    cool.dp1 = np.zeros_like(tca.x) # Pa - major losses
    cool.dp2 = np.zeros_like(tca.x) # Pa - minor losses (diam
    cool.dp3 = np.zeros_like(tca.x) # Pa - acceleration losses
    cool.dp = np.zeros_like(tca.x)  # Pa - total losses

    for k, i in enumerate(idx_range): # loop along the flow

        # Major Losses (channel wall roughness)
        cool.dp1[i] = - 0.5 * cool.f_darcy[i] * cool.density * cool.velocity[i] ** 2 * tca.dx / jacket.d_hydraulic[i] # Pa - major loss

        # Minor Losses (diameter / channel-count changes)
        if jacket.n_pipe[i+next] == jacket.n_pipe[i]:  # most of the time: same number of channels
            d_hydraulic_mod = jacket.d_hydraulic[i+next]
            d_ratio = d_hydraulic_mod / jacket.d_hydraulic[i]
            A1 = jacket.w_pipe[i] * jacket.h_pipe[i] # channel area before
            A2 = jacket.w_pipe[i+next] * jacket.h_pipe[i+next] # channel area after
        elif jacket.n_pipe[i+next] < jacket.n_pipe[i]:  # channels merge
            d_hydraulic_mod = 4.0 * 0.5 * jacket.w_pipe[i+next] * jacket.h_pipe[i+next] / (2.0 * 0.5 * jacket.w_pipe[i+next] + 2.0 * jacket.h_pipe[i+next])  # half-channel
            d_ratio = d_hydraulic_mod / jacket.d_hydraulic[i]
            A1 = jacket.w_pipe[i] * jacket.h_pipe[i] # channel area before merge
            A2 = 0.5 * jacket.w_pipe[i+next] * jacket.h_pipe[i+next] # channel area after merge
        elif jacket.n_pipe[i+next] > jacket.n_pipe[i]:  # channels divide
            d_hydraulic_mod = 4.0 * 0.5 * jacket.w_pipe[i] * jacket.h_pipe[i] / (2.0 * 0.5 * jacket.w_pipe[i] + 2.0 * jacket.h_pipe[i])  # half-channel
            d_ratio = jacket.d_hydraulic[i+next]/d_hydraulic_mod
            A1 = 0.5 * jacket.w_pipe[i] * jacket.h_pipe[i] # channel area before divide
            A2 = jacket.w_pipe[i+next] * jacket.h_pipe[i+next] # channel area after divide

        # Contraction / expansion coefficient (empirical)
        if d_ratio < 1.0:
            # Contraction
            r = float(d_ratio)
            K = 0.5 - 0.167 * r - 0.125 * r ** 2 - 0.208 * r ** 3 
        else:
            # Expansion
            r = float(d_ratio)
            K = (r ** -2.0 - 1.0) ** 2.0

        cool.dp2[i] = - 0.5 * K * cool.density * cool.velocity[i] ** 2  # Pa - Minor Losses (diameter changes)
        cool.dp3[i] = - 0.5*cool.density*(cool.velocity[i+next]**2 - cool.velocity[i]**2) # Pa - Acceleration Losses

        cool.dp[i] = cool.dp1[i] + cool.dp2[i] + cool.dp3[i]  # Pa - Total Losses

        # March pressure to the next station 
        cool.p[i+next] = cool.p[i] + cool.dp[i]

        cool.velocity[i+next] = cool.mdot / (cool.density * jacket.n_pipe[i+next] * jacket.pipe_area[i+next]) # m/s - fluid bulk speed at each station
        cool.Re[i+next] = cool.velocity[i+next] * jacket.d_hydraulic[i+next] / cool.kin_visc # unitless - Reynolds number
        cool.f_darcy[i+next] = (0.79 * np.log(cool.Re[i+next]) - 1.64) ** (-2.0) # unitless - friction factor, smooth pipe approximation REVIEW IF THIS IS ACCURATE

    cool.p = cool.p + gas.p[0] * (1.0 + tca.stiffness) - cool.p[flow_exit_index]

