import numpy as np

from engine_state import engine, tca
from rocketcea.cea_obj import CEA_Obj as rcea
from rocketprops.rocket_prop import get_prop as rprop

def coolant_flow(): # https://www.cryo-rocket.com/regenerative-cooling-circuit/6.1-ltmcc-cooling-channels/
    
    density_cool = 1000*rprop('RP1').SGLiqAtTdegR(engine.T_amb*1.8) # kg/m3 - density of RP1
    dyn_visc_cool = 0.1*rprop('RP1').ViscAtTdegR(engine.T_amb*1.8) # Pa-s - dynamic viscosity of RP1
    kin_visc_cool = dyn_visc_cool / density_cool # m2/s - kinematic viscosity of RP1
    mdot_cool = tca.mdot_fuel # kg/s

    tca.p_cool = np.full(len(tca.x), float(engine.p_amb))     # Pa - coolant pressure (initialized at p_amb)
    d_hydraulic = 4.0 * tca.w_pipe * tca.h_pipe / (2.0 * tca.w_pipe + 2.0 * tca.h_pipe)  # m - hydraulic diameter

    area = tca.w_pipe * tca.h_pipe # m2 - cross-sectional area of one channel
    tca.v_cool = mdot_cool / (density_cool * tca.n_pipe * area) # m/s - fluid bulk speed at each station
    Re_cool = tca.v_cool * d_hydraulic / kin_visc_cool # unitless - Reynolds number
    tca.f_cool = (0.79 * np.log(Re_cool) - 1.64) ** (-2.0) # unitless - Darcy friction factor

    # marching order depending on flow direction
    if tca.flow_direction == 1:
        idx_range = range(0, len(tca.x) - 1)
        flow_exit_index = len(tca.x) - 1
    else:
        idx_range = range(len(tca.x) - 1, 0, -1)
        flow_exit_index = 0

    next = tca.flow_direction  # index offset for next station
    tca.h_pipe = tca.h_pipe*np.ones_like(tca.x)

    for k, i in enumerate(idx_range): # loop along the flow

        # Major Losses (channel wall roughness)
        dp1 = 0.5 * tca.f_cool[i] * density_cool * tca.v_cool[i] ** 2 * tca.dx / d_hydraulic[i] # Pa - major loss

        # Minor Losses (diameter / channel-count changes)
        if tca.n_pipe[i+next] == tca.n_pipe[i]:  # most of the time: same number of channels
            d_hydraulic_mod = d_hydraulic[i+next]
            d_ratio = d_hydraulic_mod / d_hydraulic[i]
            A1 = tca.w_pipe[i] * tca.h_pipe[i] # channel area before
            A2 = tca.w_pipe[i+next] * tca.h_pipe[i+next] # channel area after
        elif tca.n_pipe[i+next] < tca.n_pipe[i]:  # channels merge
            d_hydraulic_mod = 4.0 * 0.5 * tca.w_pipe[i+next] * tca.h_pipe[i+next] / (2.0 * 0.5 * tca.w_pipe[i+next] + 2.0 * tca.h_pipe[i+next])  # half-channel
            d_ratio = d_hydraulic_mod / d_hydraulic[i]
            A1 = tca.w_pipe[i] * tca.h_pipe[i] # channel area before merge
            A2 = 0.5 * tca.w_pipe[i+next] * tca.h_pipe[i+next] # channel area after merge
        elif tca.n_pipe[i+next] > tca.n_pipe[i]:  # channels divide
            d_hydraulic_mod = 4.0 * 0.5 * tca.w_pipe[i] * tca.h_pipe[i] / (2.0 * 0.5 * tca.w_pipe[i] + 2.0 * tca.h_pipe[i])  # half-channel
            d_ratio = d_hydraulic[i+next]/d_hydraulic_mod
            A1 = 0.5 * tca.w_pipe[i] * tca.h_pipe[i] # channel area before divide
            A2 = tca.w_pipe[i+next] * tca.h_pipe[i+next] # channel area after divide


        # Contraction / expansion coefficient (empirical)
        if d_ratio < 1.0:
            # Contraction
            r = float(d_ratio)
            K = 0.5 - 0.167 * r - 0.125 * r ** 2 - 0.208 * r ** 3 
        else:
            # Expansion
            r = float(d_ratio)
            K = (r ** -2.0 - 1.0) ** 2.0

        dp2 = 0.5 * K * density_cool * tca.v_cool[i] ** 2  # Pa - Minor Losses (diameter changes)
        dp3 = 2/(tca.n_pipe[i]**2*density_cool*(A1+A2))*(1/A1-1/A2) * mdot_cool**2 # Pa - Acceleration Losses

        dp = dp1 + dp2 + dp3  # Pa - Total Losses

        # March pressure to the next station
        tca.p_cool[i+next] = tca.p_cool[i] - dp

        tca.v_cool[i+next] = mdot_cool / (density_cool * tca.n_pipe[i+next] * area[i+next]) # m/s - fluid bulk speed at each station
        Re_cool[i+next] = tca.v_cool[i+next] * d_hydraulic[i+next] / kin_visc_cool # unitless - Reynolds number
        tca.f_cool[i+next] = (0.79 * np.log(Re_cool[i+next]) - 1.64) ** (-2.0) # unitless - friction factor, smooth pipe approximation REVIEW IF THIS IS ACCURATE

    tca.p_cool = tca.p_cool + tca.p[0] * (1.0 + tca.stiffness) - tca.p_cool[flow_exit_index]

