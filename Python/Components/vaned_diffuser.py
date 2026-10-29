import numpy as np

from engine_state import engine, pump

def vaned_diffuser(p: pump) -> None:

    dif = p.vaned_diffuser
    imp = p.impeller[0]  # Currently only single stage pumps are supported

    ### Calculation of the flow rate
    vdot_seals = 0 # m3/s - PLACEHOLDER until leakage is modelled
    dif.vdot = p.vdot - vdot_seals # m3/s

    ### Inlet diameter
    nq = 52.9186*imp.specific_speed # unitless - Gulich nq, only for his correlations
    if nq <= 40.0:  # Gulich table 10.2
        if imp.head_rise < 100:
            dif.r_inlet = imp.r_outlet * 1.015
        else:
            dif.r_inlet = imp.r_outlet * (1.015 + 0.08*(p.density*imp.head_rise/1e6-0.1)**0.8) # m
    else:
        dif.r_inlet = imp.r_outlet * (1.04 + 0.001*(nq-40)) # m

    ### Inlet width
    dif.width_inlet = imp.width_outlet * dif.inlet_width_ratio # m

    ### Inlet velocity
    dif.v_tangential_inlet = imp.v_tangential_outlet * imp.r_outlet / dif.r_inlet  # m/s - tangential velocity at diffuser inlet
    dif.v_merid_inlet = dif.vdot / (np.pi * 2 * dif.r_inlet * dif.width_inlet) # m/s - meridional velocity at diffuser inlet (vaneless gap, no blockage)

    ###  Inlet angle and leading edge profile
    dif.flow_angle_inlet = np.arctan2(dif.v_merid_inlet, dif.v_tangential_inlet) # rad - angle between leading edge and tangential azimuth
    dif.blade_angle_inlet = dif.flow_angle_inlet + dif.incidence # rad - angle between leading edge and tangential azimuth

    ### Vane number
    allowed_vane_counts = [0,0,0,0,0,0,0,5,5,5,6,7,7,7,7] # unitless - Gulich table D7.4 - the index is the impeller blade count
    dif.vane_count = allowed_vane_counts[imp.blade.count] if imp.blade.count < len(allowed_vane_counts) else 0 # unitless
    if not dif.vane_count > 0:
        raise ValueError(f"Impeller blade count {imp.blade.count} is not currently supported for a vaned diffuser")

    ### Throat width
    X = np.exp(dif.vdot/(dif.width_inlet*imp.v_tangential_outlet*imp.r_outlet*dif.vane_count))
    dif.throat_width_inlet = dif.throat_width_factor*dif.r_inlet*(X-1) # m - Gulich eqn 7.7.10

    ### Outer diameter
    dif.r_outlet = dif.r_inlet*(dif.outlet_factor + 0.01*nq) # m - Gulich 7.9.1

    ### Diffusing channel


    ### Criteria checks


    ### Return vanes


    dif.total_loss = 0.0 # unitless - PLACEHOLDER until diffuser losses are modelled
