
from dataclasses import dataclass
import numpy as np
from scipy.integrate import cumulative_trapezoid
import rocketcea.py_cea as py_cea
from engine_state import engine, tca

def combustion_chamber_geometry():
    gas = tca.combustion_gas
    jacket = tca.regenerative_jacket
    cool = tca.regenerative_coolant
    film = tca.film_coolant
    
    # Cross Section Areas ------------------------------------------------------------------------------------------------
    tca.r_throat = np.sqrt(tca.A_throat / np.pi)    # m - Throat Radius
    tca.d_throat = 2 * tca.r_throat                 # m - Throat Diameter
    tca.A_exit   = tca.A_throat * tca.Ae_At         # m2 - Nozzle Exit Area
    tca.d_exit   = np.sqrt(4 * tca.A_exit / np.pi)  # m - Nozzle Exit Diameter

    d1_chamber = tca.d1_chamber                          # m - chamber diameter
    Ac_At = np.pi * (d1_chamber ** 2) / (4 * tca.A_throat)      # unitless - contraction Ratio
    if Ac_At < 2.0:
        raise ValueError(f"Contraction ratio is {round(((d1_chamber / tca.d_throat) ** 2), 2)}, below recommended minimum of 2\n")
    elif Ac_At > 5.0:
        raise ValueError(f"Contraction ratio is {round(((d1_chamber / tca.d_throat) ** 2), 2)}, above recommended maximum of 5\n")

    # Chamber Length -----------------------------------------------------------------------------------------------
    V_frustrum = np.pi/3 * (tca.r1_chamber**2 + tca.r1_chamber*tca.r_throat + tca.r_throat**2) * (tca.r1_chamber-tca.r_throat)/np.tan(tca.converge_angle) # m3
    tca.l_chamber = (tca.l_star*tca.A_throat - V_frustrum) / (np.pi*tca.r1_chamber**2) # m
    V_chamber = tca.l_chamber*np.pi*tca.r1_chamber**2 # m3
    x_combustor = tca.l_chamber # m

    # Contours ----------------------------------------------------------------------------------------------------
    x1_throat = 0.5*(tca.d1_chamber-tca.d_throat)/np.tan(tca.converge_angle) + tca.l_chamber
        # m - Projected point, onto central axis, of the converging straight line contour IF throat had no curvature
    x2_throat = x1_throat + tca.r_throat*(np.sin(tca.converge_angle)+1/np.tan(tca.converge_angle)*(np.cos(tca.converge_angle)-1))
        # m - Central point of the throat
    x3_throat = x2_throat + tca.r_throat*(np.sin(tca.diverge_angle)+1/np.tan(tca.diverge_angle)*(np.cos(tca.diverge_angle)-1))
        # m - Projected point, onto central axis, of the diverging straight line contour IF throat had no curvature

    tca.x_exit = x3_throat + 0.5*(tca.d_exit-tca.d_throat)/np.tan(tca.diverge_angle) # Position of nozzle exit, aka length of the entire chamber + nozzle SUSSSSSSSSS
    tca.x = np.arange(0.0, tca.x_exit, tca.dx) # m - positional domain

    # Define inner contour
    tca.r1 = np.ones_like(tca.x) * 0.5 * d1_chamber # Chamber

    # Convergence
    index_l = int(np.floor(tca.l_chamber/tca.dx))
    tca.r1[index_l:] = 0.5 * (d1_chamber - tca.d_throat) * ( x1_throat - tca.x[index_l:]) / (x1_throat - tca.x[index_l]) + 0.5*tca.d_throat # m

    # Divergence
    index_x3 = int(np.floor(x3_throat/tca.dx))
    tca.r1[index_x3:] = 0.5 * (tca.d_exit - tca.d_throat) * ( (tca.x[index_x3:] - tca.x[index_x3]) / (tca.x[-1] - tca.x[index_x3]) ) + 0.5 * tca.d_throat # m

    # Throat Arc
    x_arc = np.arange(-tca.r_throat*np.sin(tca.converge_angle), tca.r_throat*np.sin(tca.diverge_angle), tca.dx)  # m

    throat_arc = 0.5 * tca.d_throat - np.sqrt(tca.r_throat**2 - x_arc**2) + tca.r_throat # m
    index_throat = np.argmin(throat_arc) # index of throat minimum
    index_x2 = int(np.ceil(x2_throat/tca.dx))
    tca.r1[index_x2-index_throat : index_x2 + len(x_arc) - index_throat] = throat_arc

    # Define outer contour
    thickness = tca.wall_thickness # m
    tca.r2 = tca.r1 + thickness # m

    # Pipe bounds (find first/last index where tca.r1 <= merge_radius_1)
    mask_le = tca.r1 <= jacket.merge_radius
    pipe_bound1 = int(np.argmax(mask_le)) if np.any(mask_le) else 0
    pipe_bound2 = int(len(tca.x) - 1 - np.argmax(mask_le[::-1])) if np.any(mask_le) else len(tca.x) - 1

    jacket.n_pipe = np.zeros_like(tca.x, dtype=int)
    jacket.n_pipe[:pipe_bound1] = jacket.n_pipe1
    jacket.n_pipe[pipe_bound1:pipe_bound2] = jacket.n_pipe2
    jacket.n_pipe[pipe_bound2:] = jacket.n_pipe3
    if np.min(jacket.n_pipe) <= 0:
        raise ValueError("Geometric conditions for coolant channels are invalid")
    if pipe_bound2 < len(tca.x):
        jacket.n_pipe[pipe_bound2:] = jacket.n_pipe3

    # Coolant channel width -----------------------------------------------------------------------------
    jacket.w_pipe = 2 * np.pi * tca.r2 / jacket.n_pipe - jacket.gap_pipe              # m - coolant channel width
    if np.any(jacket.w_pipe <= 0):
        raise ValueError("Sections of coolant channels have zero thickness")

    #######################################################################################################
    ### Flow Properties Along Combustion Chamber and Nozzle
    #######################################################################################################

    n_chamber_end = int(np.round(tca.l_chamber/tca.dx))
    n_throat = int(np.round(x2_throat/tca.dx))

    tca.cea_station_indices = [0, n_chamber_end, n_throat, len(tca.x)-1] # indices corresponding to chamber inlet, chamber exit, throat, and nozzle exit for storing CEA outputs

    gas.gamma = np.zeros_like(tca.x)  # Ratio of Specific Heats
    gas.gamma[0] = py_cea.prtout.gammas[0]  # Ratio of Specific Heats at chamber
    gas.gamma[n_chamber_end] = py_cea.prtout.gammas[1]  # Ratio of Specific Heats at chamber exit
    gas.gamma[n_throat] = py_cea.prtout.gammas[2]  # Ratio of Specific Heats at throat
    gas.gamma[-1] = py_cea.prtout.gammas[3]  # Ratio of Specific Heats at nozzle exit
    gas.gamma = interp_zeros(gas.gamma)

    # gas.cp = np.zeros_like(tca.x)  # J/kg-K - Specific Heat Capacity at constant pressure
    # gas.cp[0] = 1000*py_cea.trpts.cpeql[0]  # J/kg-K - Specific Heat Capacity at constant pressure at injector
    # gas.cp[n_chamber_end] = 1000*py_cea.trpts.cpeql[1]  # J/kg-K - Specific Heat Capacity at constant pressure at chamber exit
    # gas.cp[n_throat] = 1000*py_cea.trpts.cpeql[2]  # J/kg-K - Specific Heat Capacity at constant pressure at throat
    # gas.cp[-1] = py_cea.trpts.cpeql[3]  # J/kg-K - Specific Heat Capacity at constant pressure at nozzle exit
    # gas.cp = interp_zeros(gas.cp)

    gas.M = np.zeros_like(tca.x)  # Mach Number                         
    for i in range(0, len(tca.x)):
        if i < n_throat: supersonic = False
        else: supersonic = True
        gas.M[i] = solve_mach((tca.r1[i]/tca.r_throat)**2, gamma=tca.gamma_avg_nozzle, supersonic=supersonic)  # Mach number at nozzle exit
    
    gas.M = interp_zeros(gas.M)

    gas.p = np.zeros_like(tca.x)  # Pa - Static Pressure
    gas.p[0] = 1E5*py_cea.prtout.ppp[0]  # Pa - chamber pressure
    gas.p[n_chamber_end] = 1E5*py_cea.prtout.ppp[1]  # Pa - chamber exit pressure
    gas.p[n_throat] = 1E5*py_cea.prtout.ppp[2]  # Pa - throat pressure
    gas.p[-1] = 1E5*py_cea.prtout.ppp[3]  # Pa - nozzle exit pressure
    gas.p = interp_zeros(gas.p)

def solve_mach(A_ratio, gamma, supersonic, tol=1e-10):
    gp1 = gamma + 1
    gm1 = gamma - 1
    exp = gp1 / (2 * gm1)  # (gamma+1)/(2*(gamma-1))

    M = (1 + np.sqrt(A_ratio - 1)) if supersonic else (0.3 / A_ratio)
    for _ in range(50):
        M2 = M**2
        t = 1 + 0.5 * gm1 * M2
        f = (1/M) * (2/gp1 * t)**exp - A_ratio
        # df/dM via product rule: d/dM[ M^{-1} * (2/gp1 * t)^exp ]
        df = -1/M2 * (2/gp1 * t)**exp + (1/M) * exp * (2/gp1 * t)**(exp-1) * (2/gp1) * gm1 * M
        # simplifies to: (2/gp1*t)^exp * (-1/M^2 + exp*gm1/(t))
        M -= f / df
        if abs(f) < tol:
            break
    return M

def interp_zeros(arr):
    """Replace zeros in a sparse array with linearly interpolated values."""
    indices = np.arange(len(arr))
    nonzero = np.nonzero(arr)[0]
    return np.interp(indices, nonzero, arr[nonzero])
