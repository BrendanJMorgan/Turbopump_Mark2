
from dataclasses import dataclass
import numpy as np
from scipy.integrate import cumulative_trapezoid
import rocketcea.py_cea as py_cea
from engine_state import engine, tca

def cc_geometry():

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
    l_chamber = (tca.l_star*tca.A_throat - V_frustrum) / (np.pi*tca.r1_chamber**2) # m
    V_chamber = l_chamber*np.pi*tca.r1_chamber**2 # m3
    x_combustor = l_chamber # m

    # Contours ----------------------------------------------------------------------------------------------------
    x1_throat = 0.5*(tca.d1_chamber-tca.d_throat)/np.tan(tca.converge_angle) + l_chamber
        # m - Projected point, onto central axis, of the converging straight line contour IF throat had no curvature
    x2_throat = x1_throat + tca.r_throat*(np.sin(tca.converge_angle)+1/np.tan(tca.converge_angle)*(np.cos(tca.converge_angle)-1))
        # m - Central point of the throat
    x3_throat = x2_throat + tca.r_throat*(np.sin(tca.diverge_angle)+1/np.tan(tca.diverge_angle)*(np.cos(tca.diverge_angle)-1))
        # m - Projected point, onto central axis, of the diverging straight line contour IF throat had no curvature

    x_exit = x3_throat + 0.5*(tca.d_exit-tca.d_throat)/np.tan(tca.diverge_angle) # Position of nozzle exit, aka length of the entire chamber + nozzle SUSSSSSSSSS
    tca.x = np.arange(0.0, x_exit, tca.dx) # m - positional domain

    # Define inner contour
    tca.r1 = np.ones_like(tca.x) * 0.5 * d1_chamber # Chamber

    # Convergence
    index_l = int(np.floor(l_chamber/tca.dx))
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
    mask_le = tca.r1 <= tca.merge_radius
    pipe_bound1 = int(np.argmax(mask_le)) if np.any(mask_le) else 0
    pipe_bound2 = int(len(tca.x) - 1 - np.argmax(mask_le[::-1])) if np.any(mask_le) else len(tca.x) - 1

    tca.n_pipe = np.zeros_like(tca.x, dtype=int)
    tca.n_pipe[:pipe_bound1] = tca.n_pipe1
    tca.n_pipe[pipe_bound1:pipe_bound2] = tca.n_pipe2
    tca.n_pipe[pipe_bound2:] = tca.n_pipe3
    if np.min(tca.n_pipe) <= 0:
        raise ValueError("Geometric conditions for coolant channels are invalid")
    if pipe_bound2 < len(tca.x):
        tca.n_pipe[pipe_bound2:] = tca.n_pipe3

    # Coolant channel width -----------------------------------------------------------------------------
    gap_pipe = tca.gap_pipe
    tca.w_pipe = 2 * np.pi * tca.r2 / tca.n_pipe - gap_pipe              # m - coolant channel width
    if np.any(tca.w_pipe <= 0):
        raise ValueError("Sections of coolant channels have zero thickness")

    # Cumulative Area ----------------------------------------------------------------------------------------
    # drdx = np.gradient(tca.r1, tca.dx) # m/m                         
    # dA = 2 * np.pi * tca.r1 * np.sqrt(1 + drdx**2) * tca.dx 
    # A_cum = cumulative_trapezoid(dA, tca.x, initial=0.0)     




    #######################################################################################################
    ### Flow Properties Along Combustion Chamber + Nozzle
    #######################################################################################################

    n_chamber_end = int(np.round(l_chamber/tca.dx))
    n_throat = int(np.round(x2_throat/tca.dx))

    tca.gamma = np.zeros_like(tca.x)  # Ratio of Specific Heats
    tca.gamma[0] = py_cea.prtout.gammas[0]  # Ratio of Specific Heats at chamber
    tca.gamma[n_chamber_end] = py_cea.prtout.gammas[1]  # Ratio of Specific Heats at chamber exit
    tca.gamma[n_throat] = py_cea.prtout.gammas[2]  # Ratio of Specific Heats at throat
    tca.gamma[-1] = py_cea.prtout.gammas[3]  # Ratio of Specific Heats at nozzle exit
    tca.gamma = interp_zeros(tca.gamma)

    tca.cp = np.zeros_like(tca.x)  # J/kg-K - Specific Heat Capacity at constant pressure
    tca.cp[0] = 1000*py_cea.trpts.cpeql[0]  # J/kg-K - Specific Heat Capacity at constant pressure at injector
    tca.cp[n_chamber_end] = 1000*py_cea.trpts.cpeql[1]  # J/kg-K - Specific Heat Capacity at constant pressure at chamber exit
    tca.cp[n_throat] = 1000*py_cea.trpts.cpeql[2]  # J/kg-K - Specific Heat Capacity at constant pressure at throat
    tca.cp[-1] = py_cea.trpts.cpeql[3]  # J/kg-K - Specific Heat Capacity at constant pressure at nozzle exit
    tca.cp = interp_zeros(tca.cp)

    tca.M = np.zeros_like(tca.x)  # Mach Number                         
    tca.M[n_chamber_end] = solve_mach(tca.Ae_At, gamma=tca.gamma_avg_nozzle, supersonic=False)  # Mach number at chamber exit
    tca.M[0] = tca.M[n_chamber_end]  # Mach number at injector face (assumed constant in chamber) COULD BE IMPROVED WITH A RAYLEIGH FLOW CORRECTION
    for i in range(n_chamber_end+1, 1, n_throat-1):
        tca.M[i] = solve_mach((tca.r1[i]/tca.r_throat)**2, gamma=tca.gamma_avg_nozzle, supersonic=False)  # Mach number at nozzle exit
    tca.M[n_throat] = 1 # Mach number at throat
    for i in range(n_throat+1, 1, -1):
        tca.M[i] = solve_mach((tca.r1[i]/tca.r_throat)**2, gamma=tca.gamma_avg_nozzle, supersonic=True)  # Mach number at nozzle exit
    tca.M = interp_zeros(tca.M)

    tca.p = np.zeros_like(tca.x)  # Pa - Static Pressure
    tca.p[0] = 1E5*py_cea.prtout.ppp[0]  # Pa - chamber pressure
    tca.p[n_chamber_end] = 1E5*py_cea.prtout.ppp[1]  # Pa - chamber exit pressure
    tca.p[n_throat] = 1E5*py_cea.prtout.ppp[2]  # Pa - throat pressure
    tca.p[-1] = 1E5*py_cea.prtout.ppp[3]  # Pa - nozzle exit pressure
    tca.p = interp_zeros(tca.p)

def solve_mach(A_ratio, gamma=1.4, supersonic=True, tol=1e-10):
    M = (1 + np.sqrt(A_ratio - 1)) if supersonic else (0.3 / A_ratio)
    for _ in range(50): # Newton-Raphson
        t = 1 + 0.5*(gamma-1)*M**2
        f = ((2/(gamma+1))*t)**((gamma+1)/(2*(gamma-1))) / M - A_ratio
        df = f/M * ((gamma+1)*M**2/(2*t) - 1)
        M -= f/df
        if abs(f) < tol: break
    return M

def interp_zeros(arr): # replaces zeros in sparse array with linearly interpolated values
    return np.interp(np.arange(len(arr)), (nz := np.nonzero(arr)[0]), arr[nz])
