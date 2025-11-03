## Objective: finds the temperature distributions of the combustion chamber and nozzle walls

import numpy as np

from engine_state import inputs, tca
from rocketcea.cea_obj import CEA_Obj as rcea
from rocketprops.rocket_prop import get_prop as rprop

def thermal_balance():
    n1 = np.floor(tca.x_combustor/tca.dx);
    n2 = np.floor(tca.x2_throat/tca.dx) - n1;  
    n3 = np.length(tca.x) - n2 - n1;

    # Pressure
    p_gas = cea.output.eql.pressure;
    segment1 = np.linspace(p_gas(1), p_gas(2), n1);
    segment2 = np.linspace(p_gas(2), p_gas(3), n2);
    segment3 = np.linspace(p_gas(3), p_gas(4), n3);
    p_gas = 1E5 * [segment1, segment2, segment3]; # Pa

    # Specific Heat, Constant Pressure
    cp_gas = cea.output.eql.cp;
    segment1 = np.linspace(cp_gas(1), cp_gas(2), n1);
    segment2 = np.linspace(cp_gas(2), cp_gas(3), n2);
    segment3 = np.linspace(cp_gas(3), cp_gas(4), n3);
    cp_gas = [segment1, segment2, segment3]; # J/(kg*K)

    # Ratio of Specific Heats
    gamma_gas = cea.output.eql.gamma;
    segment1 = np.linspace(gamma_gas(1), gamma_gas(2), n1);
    segment2 = np.linspace(gamma_gas(2), gamma_gas(3), n2);
    segment3 = np.linspace(gamma_gas(3), gamma_gas(4), n3);
    gamma_gas = [segment1, segment2, segment3]; # unitless
    
	# Mach
    # cs_gas = sqrt(gamma_gas.*p_gas./dens_gas); # m/s - speed of sound
    M1 = cea.output.eql.mach(1) + inputs.v_injection/cs_gas(1); # Mach at injector
    M2 = cea.output.eql.mach(2); # Mach at start of convergence (combustor?)
    M3 = 1 # Mach at throat, always 1
    M4 = cea.output.eql.mach(4); # Mach at exit, should be supersonic

    M_gas = zeros(1,np.length(x)); # Mach Number


def thermal_balance(inputs: ThermalInputs) -> ThermalResult:

    x = tca.x
    dx = np.gradient(x)

    # n1/n2/n3 segmentation based on special x positions
    # (MATLAB floor(...) with 1-based indexing → Python 0-based; we keep counts)
    n1 = int(np.floor(tca.x_combustor / dx[np.searchsorted(x, tca.x_combustor, side="right")-1]))
    n2 = int(np.floor(tca.x2_throat / np.mean(dx))) - n1
    n3 = len(x) - n2 - n1

    # -------------------- # Pressure
    p_gas = _segments_from_ref(cea.p_gas_ref, n1, n2, n3) * 1.0  # Pa

    # -------------------- # Specific Heat, Constant Pressure
    cp_gas_piecewise = _segments_from_ref(cea.cp_gas_ref, n1, n2, n3)  # J/(kg*K)

    # -------------------- # Ratio of Specific Heats
    gamma_gas = _segments_from_ref(cea.gamma_gas_ref, n1, n2, n3)  # unitless

    # -------------------- # Mach
    # M1, M2, M3, M4 correspond to cea.M_ref[0..3]
    M_gas = np.zeros(len(x))
    # NOTE: you had M1 ~ injector, M2 at start of convergence (combustor?), M3=1 at throat, M4>1 at exit
    # Here we just linearly morph between 4 refs like above; if you have a solver, plug it here.
    M_gas = _segments_from_ref(cea.M_ref, n1, n2, n3)

    ## Combustion Thermal Environment
    T_free = tca.Tc / (1 + (gamma_gas - 1) / 2 * M_gas**2)  # Free-Stream Temperature

    ## Wall Conduction
    fin_thickness = tca.w_pipe    # m (gap_pipe in your comment looked like width)
    fin_height = tca.h_pipe       # m
    fin_length = x.copy()       # m

    A_gas = 2.0 * np.pi * tca.r1 * dx / tca.n_pipe   # m2 - Area of dx per coolant channel
    t_wall = tca.r2 - tca.r1                          # m - Wall thickness (projected vertically)

    ## Film Flow
    mdot_film = inputs.film_fraction * tca.mdot_cc * np.ones(2)  # kept your [0.5,0.5]*... idea if split; tune as needed
    mdot_gas = (1 - inputs.film_fraction) * tca.mdot_cc

    T_injection = tca.T_amb  # K - PLACEHOLDER UNTIL REPLACED WITH REGEN TEMPERATURE AT THAT LOCATION

    # Film coolant properties (replace with rocketprops as in your code)
    T_sat_film = inputs.rprop_Tsat_film()  # K - saturation temp of film liquid (e.tca., ethanol/water mix)

    T_film = tca.T_amb * np.ones(len(x))
    # Convert film injection x-locations to index positions on x-grid
    film_inj_idx = np.clip(np.searchsorted(x, inputs.film_injection_x), 0, len(x)-1)

    ## Thermal Balance (preallocations)
    # products / fractions
    products = inputs.products
    fractions = inputs.fractions_at(x)  # shape: (n_species, len(x))

    T_wall_cold = np.zeros(len(x))  # K - wall temperature on coolant side
    Tab = np.zeros(len(x))          # K
    T_cool = tca.T_amb * np.ones(len(x))  # K
    T_ab = tca.T_amb * np.ones(len(x))    # K
    cp_gas = np.zeros(len(x))             # J/kg-K - specific heat capacity (at constant pressure) of combustion gases
    visc_gas = np.zeros(len(x))           # Pa-s - dynamic viscosity of combustion gases
    cond_gas = np.zeros(len(x))           # W/m-K - thermal conductivity of combustion gases
    h_gas = np.zeros(len(x))              # W / m2-K - convective heat transfer coefficient - MUST be initialized as all zeros for film_cooling() to work
    q_gas = np.zeros(len(x))
    Pr_gas = np.zeros(len(x))
    Nu_cool = np.zeros(len(x))
    h_cool = np.zeros(len(x))

    warned = False
    start_flag = True  # (kept but unused; left to match your original structure)

    # Reference temperature for gas transport properties
    T_ref = T_free * (1 + 0.032 * M_gas**2 + 0.58 * (T_wall_cold / np.maximum(T_free, 1e-9) - 1))  # T_wall_hot not known yet; start with cold side
    # Get mixture properties along x (loop because you had a per-x species mix call)
    for i in range(len(x)):
        # mixture(products, fractions[:, i], T_ref[i], p_gas[i]) -> (cp, mu, k)
        cp_i, mu_i, k_i = inputs.mixture_fn(products, fractions[:, i], float(T_ref[i]), float(p_gas[i]))
        cp_gas[i] = cp_i
        visc_gas[i] = mu_i
        cond_gas[i] = k_i
    Pr_gas = np.divide(cp_gas * visc_gas, cond_gas, out=np.zeros_like(cp_gas), where=cond_gas > 0)  # Prandtl Number

    # # Bartz relation
    sigma = 1.0 / (
        (0.5 * (T_wall_cold / tca.Tc.clip(min=1e-9)) * ((1 + (gamma_gas - 1) / 2.0 * M_gas**2) + 0.5)) ** 0.68
        * ((1 + (gamma_gas - 1) / 2.0 * M_gas**2) ** 0.12)
    )  # NOTE: This keeps your structure, but you'll likely replace with your validated Bartz form.

    h_gas = (
        (0.026 / (tca.d_throat ** 0.2))
        * (visc_gas ** 0.2)
        * (cp_gas ** 0.6)
        * ((p_gas / tca.c_star_eff2) ** 0.8)
        * ((tca.d_throat / tca.rc_throat) ** 0.1)
        * sigma
    )

    # K - Adiabatic Wall Temperature. 0.33 ==> 0.5 IF GAS ...
    T_ab = tca.Tc * (1 + Pr_gas ** 0.33 * ((gamma_gas - 1) / 2) * M_gas ** 2) / (1 + ((gamma_gas - 1) / 2) * M_gas ** 2)
    # K - recovery temperature of combustion gas (fluid state if it slowed to stagnation)
    T_recovery = tca.c_star_eff2**2 * (T_free + Pr_gas ** (1/3) * (tca.Tc - T_free))

    # Film phase map: 0 means film coolant is a gas, and 1 means it is a liquid at that location
    film_phase = np.zeros(len(x), dtype=int)

    # -------------------- Liquid Film Coolant Distribution
    # We step each injection individually; if multiple, we stack effects (keep your approach)
    for inj_i, inj_idx in enumerate(film_inj_idx):
        # energy required to boil a slug injected at inj_idx until it turns to gas
        # (uses your placeholders)
        heat_vap_film = 0.0  # TODO: wire latent heat
        cp_film_inj = inputs.rprop_cp_film(T_sat_film)  # J/kg-K (evaluated at Tsat as placeholder)
        boiling_energy = inputs.injection_efficiency * mdot_film[min(inj_i, len(mdot_film)-1)] * (
            cp_film_inj * (T_sat_film - T_injection) + heat_vap_film
        )

        # accumulate "recovery" energy to see when we’ve boiled off
        # integrate (T_recovery - T_film_prev) * ... along j; simple forward-Euler here
        accumulated_energy = 0.0
        liquid_end = inj_idx
        T_film[inj_idx] = tca.T_amb  # initialize at injection

        for j in range(inj_idx, len(x)):
            # Simple gas-side heat pickup proxy (your commented Bartz-based dt); keep your original line & constants:
            # dT step for film due to gas-side h_gas
            # NOTE: Units sanity check recommended here once wired to your property funcs
            dt_film = (
                dx[j]
                * 1.628 * (2 * np.pi * tca.r1[j]) * h_gas[j] / np.maximum(cp_film_inj, 1e-9)
                * ((visc_gas[j]) / np.maximum(tca.v_injection, 1e-9))
                * (1.0 / (2 * np.pi * tca.r1[j]))
                * ((cp_film_inj / max(inputs.rprop_k_film(T_film[max(j-1, inj_idx)]), 1e-9)) ** 0.125)
                * (T_recovery[j] - T_film[max(j-1, inj_idx)])
            )
            T_film[j] = T_film[max(j-1, inj_idx)] + dt_film

            # crude energy accounting to find boiling front (replace with your actual)
            accumulated_energy += max(T_recovery[j] - T_film[j], 0.0)
            if accumulated_energy >= boiling_energy:
                liquid_end = j
                break

        # Mark phases in this injected region
        film_phase[inj_idx:liquid_end + 1] = 1  # 1 for liquid and 0 for gas
        # After liquid_end, film coolant assumed at saturation temperature (gas after)
        if liquid_end + 1 < len(x):
            T_film[liquid_end + 1:] = np.maximum(T_film[liquid_end + 1:], T_sat_film)

    # -------------------- March coolant (regen) along flow direction
    flow_dir = inputs.flow_direction  # +1 or -1
    if flow_dir == 1:
        march_indices = range(0, len(x))
    else:
        march_indices = range(len(x) - 1, -1, -1)

    # Initialize coolant inlet temperature
    T_input = tca.T_amb

    for i in march_indices:
        ## Regenerative Coolant Properties (evaluate at current T_cool)
        Tci = T_cool[i] if i != (0 if flow_dir == 1 else len(x)-1) else T_input
        cp_cool = inputs.rprop_cp_cool(Tci)        # J/kg-K
        visc_cool = inputs.rprop_mu_cool(Tci)      # Pa-s
        k_cool = inputs.rprop_k_cool(Tci)          # W/m-K
        Pr_cool = cp_cool * visc_cool / max(k_cool, 1e-12)  # Prandtl Number

        # Nusselt number correlations (Gnielinski / laminar interpolation)
        if tca.Re_cool[i] < 3000:  # Really should be < 2300, but using laminar flow in transition region for conservatism
            lam_x = np.array([1, 1.43, 2, 3, 4, 8, 1e10, 3.61, 3.73, 4.12, 4.79, 5.33, 6.49, 8.23])
            lam_y = lam_x.copy()   # placeholder—MATLAB used two rows: provide your table here
            # approximate laminar Nu in rectangular duct by interpolation on h_pipe/w_pipe
            ratio = np.maximum(tca.h_pipe[i] / np.maximum(tca.w_pipe[i], 1e-9), 1e-6)
            Nu_cool[i] = np.interp(ratio, lam_x, lam_y)
            if not warned:
                # fprintf("Warning: regen coolant goes laminar at #g m from injector (#g m from exit plane)\n", x(i), x_exit-x(i))
                warned = True
        elif (Pr_cool >= 0.7 and Pr_cool <= 2000 and tca.Re_cool[i] >= 3000 and tca.Re_cool[i] <= 5e6):
            f_cole = 0.316 * tca.Re_cool[i] ** (-0.25)  # simple Blasius for smooth tubes; replace with your f_cool
            Nu_cool[i] = (f_cole / 8.0) * (tca.Re_cool[i] - 1000.0) * Pr_cool / (1 + 12.7 * (f_cole / 8.0) ** 0.5 * (Pr_cool ** (2/3) - 1))
        else:
            raise ValueError("No correlation for coolant Reynold's and/or Prandtl number")

        # Convective heat transfer coefficient
        h_cool[i] = Nu_cool[i] * k_cool / max(tca.d_hydraulic[i], 1e-12)  # W/m2-K

        # Fin efficiency
        m_coeff = np.sqrt(2.0 * h_cool[i] / max(m.k_wall * max(fin_thickness[i], 1e-12), 1e-12))  # coefficient for fin efficiency equation
        fin_eff = np.tanh(m_coeff * (tca.h_pipe[i] + fin_thickness[i]) / 2) / (m_coeff * (tca.h_pipe[i] + fin_thickness[i]) / 2)
        A_cool = (2.0 * fin_height[i] * fin_eff * tca.w_pipe[i]) * dx[i]  # m2 - adjusted contact area of coolant on channel walls

        ## Heat Balancing
        if film_phase[i] == 1:    # Liquid State
            # q_wall: conduction from film to coolant through wall and gas-side film convection
            q_wall = (T_film[i] - T_cool[i]) / (1.0 / (h_cool[i] * A_cool) + 1.0 / (m.k_wall * A_gas[i] / max(t_wall[i], 1e-12)))  # W
        else:                      # Gaseous State
            q_wall = (T_film[i] - T_cool[i]) / (
                1.0 / (h_gas[i] * A_gas[i]) + 1.0 / (m.k_wall * A_gas[i] / max(t_wall[i], 1e-12)) + 1.0 / (h_cool[i] * A_cool)
            )  # W

        # Split wall temperatures (hot/cold) across the wall thickness with 1-D conduction
        T_wall_cold[i] = q_wall / (h_cool[i] * A_cool) + T_cool[i]  # K
        T_wall_hot_i = q_wall / (m.k_wall * A_gas[i] / max(t_wall[i], 1e-12)) + T_wall_cold[i]  # K

        # March coolant temperature
        dT = q_wall * tca.n_pipe[i] / max(cp_cool * tca.mdot_fuel_cc, 1e-12)  # Coolant Temperature change at one channel
        nxt = i + flow_dir
        if 0 <= nxt < len(x):
            T_cool[nxt] = T_cool[i] + dT

        # SHOULD ADD HEAT OF VAPORIZATION - THIS WILL DELAY BOILING
        # Check if the regen coolant boils. Gaseous coolant will be extremely ineffective, as well as harder to analyze
        # boiling_cool = rprop('T','p', p_cool(i), 'Q', 0,'Ethanol') # K - Ethanol Saturation Temperature
        # if T_cool(i) > boiling_cool: error("Coolant starts boiling at #g m from injector (#g m from exit plane)", x(i), x_exit-x(i))

        # Save hot wall
        if i == 0:
            T_wall_hot = np.full(len(x), T_wall_hot_i)
        else: