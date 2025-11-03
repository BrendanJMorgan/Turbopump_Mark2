from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from engine_state import engine, pump
from typing import List
import math

def impellers(p: pump):
    # Impeller Contours - pg 2.29-2.32 in pump handbook
    r_eye_margin = 1.15 + 0.2116*(p.specific_speed-0.2836); # unitless - margin correction factor, Gulich 7.1.4; for a normal impeller (not a suction impeller)
    swirl_number = 1; # NEEDS TO BE A FUNCTION OF INDUCER OUTLET CONDITIONS (ITERATIVE)
    r_inlet_impeller = r_eye_margin * np.sqrt(p.r_hub^2 + 0.0026418 * p.outlet_flow_coeff* p.specific_speed^1.33 / swirl_number^0.67) 
        # m - impeller inlet radius, Gulich eqn 7.1.4

    r_outlet_impeller = 1 / abs(p.pump_shaft_speed) * np.sqrt (engine.g*p.head_pump / p.head_coeff); # m - impeller exit radius, Gulich eqn 7.1.3

    # The exit width b2 ------------------------------------------------------------------------------
    w_exit = 2*r_outlet_impeller * (0.017 + 0.1386*p.specific_speed - 0.022387*p.specific_speed^2 + 0.0013767*p.specific_speed^3); 
        # m - outlet width; empirical; Gulich eqn 7.1

    # this needs boundary layer and leakage effect - pg 2.60?
    # could parameterize a width at each station with a blockage vector, linear ramp of area?

    # Hub and shroud profiles ------------------------------------------------------------------
    r_min = 0.5*r_inlet_impeller; # m - minimum allowable radius of curvature

    min_radius = 1/max(p.min_shroud_curvature) # m - minimum radius of curvature on shroud contour
    if min_radius < 0.99*r_min:
        print("\nPump shroud curvature reaches %.2f mm radius. Need above ~%.2f mm (half of eye radius) to avoid cavitation.", min_radius*1000, r_min*1000)
