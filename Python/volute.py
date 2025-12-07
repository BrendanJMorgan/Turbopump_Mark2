import numpy as np

from engine_state import engine, pump

def volute(p) -> None:
    # Pump Handbook page 2.68

    p.r_tongue = p.r_outlet_impeller * (1.03 + 0.0025 * p.specific_speed + 7e-8 * p.density * p.head_rise)  # m - Gulich table 10.2, pg 729
    p.t_tongue = 0.04 * p.r_outlet_impeller  # m - Gulich pg 493
    p.A_throat = (p.vdot / (3.8 * p.r_outlet_impeller * p.v_tangential_outlet) * (1 + np.sqrt(1 + 15.2 * p.r_outlet_impeller * p.v_tangential_outlet * (p.r_tongue + p.t_tongue) / p.vdot))) ** 2  # m^2
    p.r_throat = p.r_outlet_impeller + 0.5 * np.sqrt(p.A_throat) + p.t_tongue  # m

    p.theta_volute = np.linspace(0.0, 2.0 * np.pi, 1000)  # rad
    p.A_volute = (p.vdot * p.theta_volute / (7.6 * np.pi * p.r_outlet_impeller * p.v_tangential_outlet) #* (1 + np.sqrt(1 + 30.4 * np.pi * p.r_outlet_impeller * p.r_tongue * p.v_tangential_outlet / (p.vdot * p.theta_volute))) ) ** 2  # m^2
    p.A_volute[np.isnan(p.A_volute)] = 0.0

    p.h_volute = np.sqrt(0.25 * p.t_tongue ** 2 + p.A_volute[-1]) - 0.5 * p.t_tongue  # m - height of volute gap - set so the outlet is a square
    p.r_volute = p.r_tongue + p.A_volute / p.h_volute  # m - radius of volute wall
    p.d_volute_outlet = 2.0 * np.sqrt(p.A_volute[-1] / np.pi)  # m - diameter of volute outlet hole
    p.r_volute_outlet = p.r_volute[-1] + p.h_volute / 2.0  # m

    p.v_throat = p.vdot / p.A_throat  # m/s - fluid outlet velocity into plumbing

    p.volute_curve = np.column_stack([p.r_volute * np.cos(p.clocking * p.theta_volute),
                                      p.r_volute * np.sin(p.clocking * p.theta_volute)])  # [m, m]

    # "...considerable freedom to configure the cross section without risking major losses ... flat cross sections result in less intense secondary flow than circular cross
    # sections, thereby generating fewer losses ... a ratio of width B to height H in the range of B/H = 2 to 3 may be considered the optimum." - Gulich pg 492
