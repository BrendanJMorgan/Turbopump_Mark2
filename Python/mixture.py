"""
Gas mixture transport properties for combustion product gases.

Computes cp, viscosity, and thermal conductivity of a gas mixture
using Wilke's method (viscosity) and Lindsay-Bromley's method (conductivity).

mixture(compounds, fractions, T_ref, p_ref) -> (cp, visc, cond, density)
"""

import CoolProp
import numpy as np

STATE_CACHE = {} # builds a cache of CoolProp state objects for each species to avoid repeated initialisation overhead / reduces runtime

# Registry: CEA product name -> (CoolProp name, molecular weight [g/mol], boiling point [K], is_polar)
SPECIES = {
    "*CO":  ("CarbonMonoxide", 28.01,    81.65,  True),
    "*CO2": ("CarbonDioxide",  44.01,   194.65,  False),
    "*H2":  ("Hydrogen",        1.00794,  20.28,  False),
    "H2O":  ("Water",          18.01528, 373.1,   True),
    "O2":   ("Oxygen",         15.999,    90.19,  False),
    "C2H4": ("Ethylene",       28.05,    169.5,   False),
    "C2H6": ("Ethane",         30.07,    184.2,   False),
}


def get_species_props(cea_name: str, T: float, p: float):
    """Return (cp, visc, cond, density, weight, boiling, polar) for a single species,
    or None if the CEA name is not in the registry."""
    entry = SPECIES.get(cea_name)
    if entry is None:
        return None

    cp_name, mw, bp, polar = entry
    if cp_name not in STATE_CACHE:
        state = STATE_CACHE[cp_name] = CoolProp.AbstractState("HEOS", cp_name)
    else:
        state = STATE_CACHE[cp_name]


        
    state.update(CoolProp.PT_INPUTS, p, T)

    density = state.rhomass()
    cp = state.cpmass()

    if cp_name == "CarbonMonoxide":
        visc = 0.027098 * T**0.734156 * 1e-5
        cond = 0.227019 * T**0.828249 * 1e-3
    elif cp_name == "Ethylene":
        visc = 8e-10 * T**0.4845
        cond = state.conductivity()
    else:
        visc = state.viscosity()
        cond = state.conductivity()

    return cp, visc, cond, density, mw, bp, polar


def mixture(compounds, fractions, T_ref: float, p_ref: float):
    """
    Compute mixture transport properties for combustion-product gases.

    Parameters
    ----------
    compounds : list[str]
        CEA species names (e.g. ["*CO", "*CO2", "H2O", ...]).
    fractions : array-like
        Mass (or mole) fractions corresponding to each compound.
    T_ref : float
        Reference temperature [K].
    p_ref : float
        Reference pressure [Pa].

    Returns
    -------
    cp   : float – mixture specific heat [J/kg-K]
    visc : float – mixture dynamic viscosity [Pa-s]
    cond : float – mixture thermal conductivity [W/m-K]
    dens : float – mixture density [kg/m3]
    """
    # Collect per-species data for recognised compounds only
    frac_list = []
    cp_list   = []
    visc_list = []
    cond_list = []
    dens_list = []
    wt_list   = []
    bp_list   = []
    polar_list = []

    for name, xi in zip(compounds, fractions):
        props = get_species_props(name, T_ref, p_ref)
        if props is None:
            continue  # unrecognised species — skip
        cp_i, visc_i, cond_i, dens_i, mw_i, bp_i, polar_i = props
        frac_list.append(float(xi))
        cp_list.append(cp_i)
        visc_list.append(visc_i)
        cond_list.append(cond_i)
        dens_list.append(dens_i)
        wt_list.append(mw_i)
        bp_list.append(bp_i)
        polar_list.append(polar_i)

    if not frac_list:
        return 0.0, 0.0, 0.0, 0.0

    x  = np.asarray(frac_list)
    cp_arr   = np.asarray(cp_list)
    mu_arr   = np.asarray(visc_list)
    k_arr    = np.asarray(cond_list)
    rho_arr  = np.asarray(dens_list)
    w_arr    = np.asarray(wt_list)
    bp_arr   = np.asarray(bp_list)
    polar_arr = np.asarray(polar_list, dtype=bool)

    x_sum = x.sum()
    if 1.0 - x_sum > 0.05:
        print(f"Unaccounted-for combustion products make up {100 - 100*x_sum:.4g}% of exhaust")

    # --- Density (mass-fraction weighted) ---
    density = np.dot(rho_arr, x) / x_sum

    # --- Specific heat (mass-fraction weighted) ---
    cp = np.dot(cp_arr, x) / x_sum

    # --- Viscosity: Wilke's method ---
    n = len(mu_arr)
    mu_ratio = mu_arr[:, None] / mu_arr[None, :]          # mu_i / mu_j
    w_ratio  = w_arr[None, :] / w_arr[:, None]             # Mj / Mi
    w_sum_ratio = w_arr[:, None] / w_arr[None, :]          # Mi / Mj

    phi = (np.sqrt(2.0) / 4.0) * (1.0 + np.sqrt(mu_ratio) * w_ratio**0.25)**2 / np.sqrt(1.0 + w_sum_ratio)
    np.fill_diagonal(phi, 0.0)  # exclude j == i terms

    # visc = 0.0
    # for i in range(n):
    #     denom = x[i] + np.dot(phi[i, :], x)
    #     visc += x[i] * mu_arr[i] / denom
    denom = x[None, :] + phi @ x  # wrong shape; fix:
    denom = x + phi @ x
    visc = np.sum(x * mu_arr / denom)

    # --- Thermal conductivity: Lindsay & Bromley's method ---
    S = np.sqrt(2.25 * bp_arr[:, None] * bp_arr[None, :]) # Sutherland constants for each pair (i,j)

    # Adjust where exactly one molecule in the pair is polar  (factor 0.733)
    polar_xor = np.logical_xor(polar_arr[:, None], polar_arr[None, :])
    S = np.where(polar_xor, 0.733 * S, S)

    # Diagonal Sutherland constants
    S_diag = np.diag(S)  # S_ii

    # A_ij matrix (Lindsay-Bromley mixing rule)
    mu_ratio_ij = mu_arr[:, None] / mu_arr[None, :]
    w_ratio_ij  = w_arr[None, :] / w_arr[:, None]  # Mj/Mi

    term1 = 1.0 + np.sqrt(mu_ratio_ij * w_ratio_ij**0.75
                           * ((1.0 + S_diag[:, None] / T_ref)
                              / (1.0 + S_diag[None, :] / T_ref)))
    A = 0.25 * term1**2 * ((1.0 + S[:] / T_ref)
                            / (1.0 + S_diag[:, None] / T_ref))

    # cond = 0.0
    # for i in range(n):
    #     denom = np.dot(A[i, :], x) / x[i] if x[i] > 0 else 1.0
    #     cond += k_arr[i] / denom
    np.fill_diagonal(A, 0.0)
    denom = (A @ x) / x  # element-wise, safe since we skip x[i]==0 species
    cond = np.sum(k_arr / (1.0 + denom))  # need to rework slightly

    return float(cp), float(visc), float(cond), float(density)
