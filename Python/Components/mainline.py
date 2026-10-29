import numpy as np
from engine_state import pump
from volute import discharge_cp

def mainline(p: pump):
    ml = p.mainline

    #######################################################################################################
    ### Option Matrices
    #######################################################################################################
    tube_id_options = ml.tube_od_options - 2*ml.tube_wall_options # m
    inner_area = np.pi/4*tube_id_options**2 # m2
    bulk_velocity = p.mdot / (p.density*inner_area) # m/s

    reynolds = p.density*bulk_velocity*tube_id_options/p.viscosity # unitless
    f_darcy = (-1.8*np.log10((ml.roughness/tube_id_options/3.7)**1.11 + 6.9/reynolds))**-2 # unitless - Haaland

    q = 0.5*p.density*bulk_velocity**2 # Pa - dynamic pressure
    dp_major = f_darcy*ml.length/tube_id_options*q # Pa - Darcy-Weisbach
    dp_minor = ml.minor_loss_coeff*q # Pa - bends, valve, manifold entry
    dp = dp_major + dp_minor # Pa

    mass = ml.tube_density[:, None]*np.pi/4*(ml.tube_od_options**2 - tube_id_options**2)*ml.length # kg
    hoop_stress = (ml.p_out + dp)*ml.tube_od_options / (2*ml.tube_wall_options) # Pa - thin wall, at pump discharge pressure

    #######################################################################################################
    ### Volute Discharge Nozzle
    #######################################################################################################
    # Only once the volute throat is known (after the first pump run). Nozzle loss is already in pump efficiency,
    # so it only steers the tube choice here and isn't added to p_out.
    nozzle_mass = np.zeros_like(mass) # kg
    nozzle_dp = np.zeros_like(dp) # Pa
    vol = vars(p).get('volute')
    if vol is not None and hasattr(vol, 'A_throat'):
        area_ratio = inner_area/vol.A_throat # unitless
        r_throat = vol.d_throat/2 # m
        r_exit = tube_id_options/2 # m
        nozzle_length = np.maximum(r_exit - r_throat, 0)/np.tan(vol.discharge_cone_angle/2) # m
        slant = np.sqrt(nozzle_length**2 + (r_exit - r_throat)**2) # m
        nozzle_mass = ml.tube_density[:, None]*np.pi*(r_throat + r_exit)*slant*ml.tube_wall_options # kg - thin conical shell, same wall as the tube
        cp = np.array([nozzle_cp(L/r_throat, A) for L, A in zip(nozzle_length, area_ratio)]) # unitless
        nozzle_dp = 0.5*p.density*vol.c_throat**2*(1 - cp - 1/area_ratio**2) # Pa - Gulich 3.8.22
        nozzle_dp[area_ratio < 1] = np.inf # contracting nozzle not allowed
        nozzle_dp[nozzle_length > vol.discharge_max_length] = np.inf # packaging limit

    #######################################################################################################
    ### Optimization
    #######################################################################################################
    cost = mass + nozzle_mass + ml.dp_mass_equivalent*(dp + nozzle_dp) # kg
    cost[hoop_stress > ml.allowable_stress[:, None]] = np.inf
    if np.all(np.isinf(cost)):
        raise ValueError(f"No standard {p.fluid} mainline tube holds {(ml.p_out + dp.min())/1e5:.1f} bar")
    material_index, size_index = np.unravel_index(np.argmin(cost), cost.shape)

    ml.material = ml.tube_material_options[material_index]
    ml.outer_diameter = ml.tube_od_options[size_index] # m
    ml.wall = ml.tube_wall_options[size_index] # m
    ml.tube_id = tube_id_options[size_index] # m
    ml.inner_area = inner_area[size_index] # m2
    ml.bulk_velocity = bulk_velocity[size_index] # m/s
    ml.reynolds = reynolds[size_index] # unitless
    ml.f_darcy = f_darcy[size_index] # unitless
    ml.dp_major = dp_major[size_index] # Pa
    ml.dp_minor = dp_minor[size_index] # Pa
    ml.dp = dp[size_index] # Pa
    ml.mass = mass[material_index, size_index] # kg
    ml.hoop_stress = hoop_stress[size_index] # Pa
    ml.nozzle_mass = nozzle_mass[material_index, size_index] # kg
    ml.nozzle_dp = nozzle_dp[size_index] # Pa

    p.d_mainline = ml.tube_id # m - discharge diffuser exit
    return ml.p_out + ml.dp # Pa - pump discharge pressure

def nozzle_cp(L_R1: float, A_R: float) -> float:
    if A_R < 1:
        return np.nan # contracting - rejected by the caller
    try:
        return discharge_cp(L_R1, A_R)
    except ValueError:
        return 0.0 # too short or too little area change to recover anything off the chart
