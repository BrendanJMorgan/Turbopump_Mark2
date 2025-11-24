
from dataclasses import dataclass
import numpy as np
from scipy.integrate import cumulative_trapezoid
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
        raise ValueError("Contraction ratio is #g, below recommended minimum of 2\\n", (d1_chamber / tca.d_throat) ** 2)
    elif Ac_At > 5.0:
        raise ValueError("Contraction ratio is #g, above recommended maximum of 5\\n", (d1_chamber / tca.d_throat) ** 2)

    # Chamber Length -----------------------------------------------------------------------------------------------
    V_frustrum = np.pi/3 * (tca.r1_chamber**2 + tca.r1_chamber*tca.r_throat + tca.r_throat**2) * \
                 (tca.r1_chamber-tca.r_throat)/np.tan(tca.converge_angle); # m3
    l_chamber = (tca.l_star*tca.A_throat - V_frustrum) / (np.pi*tca.r1_chamber**2); # m
    V_chamber = l_chamber*np.pi*tca.r1_chamber**2; # m3
    x_combustor = l_chamber; # m

    # Contours ----------------------------------------------------------------------------------------------------
    x1_throat = 0.5*(tca.d1_chamber-tca.d_throat)/np.tan(tca.converge_angle) + l_chamber; 
        # m - Projected point, onto central axis, of the converging straight line contour IF throat had no curvature
    x2_throat = x1_throat + tca.r_throat*(np.sin(tca.converge_angle)+1/np.tan(tca.converge_angle)*(np.cos(tca.converge_angle)-1));
        # m - Central point of the throat
    x3_throat = x2_throat + tca.r_throat*(np.sin(tca.diverge_angle)+1/np.tan(tca.diverge_angle)*(np.cos(tca.diverge_angle)-1));
        # m - Projected point, onto central axis, of the diverging straight line contour IF throat had no curvature

    x_exit = x3_throat + 0.5*(tca.d_exit-tca.d_throat)/np.tan(tca.diverge_angle); # Position of nozzle exit, aka length of the entire chamber + nozzle SUSSSSSSSSS
    tca.x = np.arange(0.0, x_exit + 0.5 * tca.dx, tca.dx) # m - positional domain

    # Define inner contour
    r1 = np.ones_like(tca.x) * 0.5 * d1_chamber # Chamber

    # Convergence
    index_l = int(np.floor(l_chamber/tca.dx))
    r1[index_l:] = 0.5 * (d1_chamber - tca.d_throat) * ( (tca.x[index_l:] - tca.x[index_l]) / (x1_throat - tca.x[index_l]) ) + 0.5 * tca.d_throat # m

    # Divergence
    index_x3 = int(np.floor(x3_throat/tca.dx))
    r1[index_x3:] = 0.5 * (tca.d_exit - tca.d_throat) * ( (tca.x[index_x3:] - tca.x[index_x3]) / (tca.x[-1] - tca.x[index_x3]) ) + 0.5 * tca.d_throat # m

    # Throat Arc
    x_arc = -tca.r_throat * np.sin(tca.converge_angle) + np.arange(0, int(tca.r_throat * np.sin(tca.diverge_angle) / tca.dx)) * tca.dx
    throat_arc = 0.5 * tca.d_throat - np.sqrt(tca.r_throat**2 - x_arc**2) + tca.r_throat
    index_throat = np.argmin(throat_arc)  # or your exact criterion
    index_x2 = int(np.ceil(x2_throat/tca.dx))
    r1[index_x2-index_throat : index_x2 + len(x_arc) - index_throat] = throat_arc

    # Define outer contour
    thickness = tca.wall_thickness # m
    r2 = r1 + thickness # m

    # Pipe bounds (find first/last index where r1 <= merge_radius_1)
    mask_le = r1 <= tca.merge_radius
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

    # Coolant channel width check -----------------------------------------------------------------------------
    gap_pipe = tca.gap_pipe
    tca.w_pipe = 2 * np.pi * r2 / tca.n_pipe - gap_pipe              # m - coolant channel width
    if np.any(tca.w_pipe <= 0):
        raise ValueError("Sections of coolant channels have zero thickness")

    # Cumulative Area ----------------------------------------------------------------------------------------
    drdx = np.gradient(r1, tca.dx)                            # MATLAB gradient(r1,dx)
    dA = 2 * np.pi * r1 * np.sqrt(1 + drdx**2) * tca.dx
    A_cum = cumulative_trapezoid(dA, tca.x, initial=0.0)      # MATLAB cumtrapz(x, dA)

