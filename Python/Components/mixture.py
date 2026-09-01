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
    mw = mw * 0.001  # convert g/mol to kg/mol
    if cp_name not in STATE_CACHE:
        state = STATE_CACHE[cp_name] = CoolProp.AbstractState("HEOS", cp_name)
    else:
        state = STATE_CACHE[cp_name]
        
    state.update(CoolProp.PT_INPUTS, p, T)

    density = state.rhomass()
    specific_heat= state.cpmass()

    if cp_name == "CarbonMonoxide":
        viscosity = 0.027098 * T**0.734156 * 1e-5
        cond = 0.227019 * T**0.828249 * 1e-3
    elif cp_name == "Ethylene":
        viscosity = 8e-10 * T**0.4845
        cond = state.conductivity()
    else:
        viscosity = state.viscosity()
        cond = state.conductivity()

    return specific_heat, viscosity, cond, density, mw, bp, polar


def mixture(cea_compounds, fractions, T_ref: float, p_ref: float):
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
    specific_heat  : float – mixture specific heat [J/kg-K]
    viscosity: float – mixture dynamic viscosity [Pa-s]
    cond : float – mixture thermal conductivity [W/m-K]
    dens : float – mixture density [kg/m3]
    """
    # Collect per-species data for recognised compounds only
    fractions_list,cp_list,visc_list,cond_list,dens_list,mol_weight_list,bp_list,polar_list = [],[],[],[],[],[],[],[]

    for cea_compound, fraction in zip(cea_compounds, fractions):
        properties = get_species_props(cea_compound, T_ref, p_ref) # pulling properties with coolprop and empirical fits
        if properties is None:
            continue  # unrecognised species — skip
        cp_i, visc_i, cond_i, dens_i, mw_i, bp_i, polar_i = properties
        fractions_list.append(float(fraction))
        cp_list.append(cp_i)
        visc_list.append(visc_i)
        cond_list.append(cond_i)
        dens_list.append(dens_i)
        mol_weight_list.append(mw_i)
        bp_list.append(bp_i)
        polar_list.append(polar_i)

    if not fractions_list:
        return 0.0, 0.0, 0.0, 0.0

    fractions_array  = np.asarray(fractions_list)
    density_array    = np.asarray(cp_list)
    viscosity_array  = np.asarray(visc_list)
    k_array          = np.asarray(cond_list)
    density_array    = np.asarray(dens_list)
    mol_weight_array = np.asarray(mol_weight_list)
    bp_array         = np.asarray(bp_list)
    polar_array      = np.asarray(polar_list, dtype=bool)

    recognized = fractions_array.sum()                                       # species found in SPECIES registry
    total = float(np.sum(np.asarray(fractions, dtype=float)))  # everything passed in
    if total > 0 and (total - recognized) / total > 0.10:
        print(f"Unaccounted-for combustion products make up {100*(total-recognized)/total:.4g}% of exhaust")

    molar_fractions = fractions_array / fractions_array.sum()                                    # mole fractions
    mass_fractions = molar_fractions * mol_weight_array / np.dot(molar_fractions, mol_weight_array)                   # mass fractions
    specific_heat     = np.dot(density_array, mass_fractions)                         # J/kg-K, mass-weighted
    density = p_ref * np.dot(molar_fractions, mol_weight_array) / (8.3145 * T_ref)   # ideal gas, kg/m3

    #######################################################################################################
    ### Viscosity: Wilke's method
    #######################################################################################################
    n = len(viscosity_array)
    viscosity_ratio = viscosity_array[:, None] / viscosity_array[None, :]          # viscosity_i / viscosity_j
    mol_weight_ratio  = mol_weight_array[None, :] / mol_weight_array[:, None]             # Mj / Mi
    w_sum_ratio = mol_weight_array[:, None] / mol_weight_array[None, :]          # Mi / Mj

    phi = (np.sqrt(2.0) / 4.0) * (1.0 + np.sqrt(viscosity_ratio) * mol_weight_ratio**0.25)**2 / np.sqrt(1.0 + w_sum_ratio)
    np.fill_diagonal(phi, 0.0)  # exclude j == i terms

    denom = fractions_array[None, :] + phi @ fractions_array  # wrong shape; fix:
    denom = fractions_array + phi @ fractions_array
    viscosity = np.sum(fractions_array * viscosity_array / denom)

    #######################################################################################################
    ### Thermal conductivity: Lindsay & Bromley's method
    #######################################################################################################
    S = np.sqrt(2.25 * bp_array[:, None] * bp_array[None, :]) # Sutherland constants for each pair (i,j)

    # Adjust where exactly one molecule in the pair is polar  (factor 0.733)
    polar_xor = np.logical_xor(polar_array[:, None], polar_array[None, :])
    S = np.where(polar_xor, 0.733 * S, S)

    # Diagonal Sutherland constants
    S_diag = np.diag(S)  # S_ii

    # A_ij matrix (Lindsay-Bromley mixing rule)
    viscosity_ratio_ij = viscosity_array[:, None] / viscosity_array[None, :]
    mol_weight_ratio_ij  = mol_weight_array[None, :] / mol_weight_array[:, None]  # Mj/Mi

    term1 = 1.0 + np.sqrt(viscosity_ratio_ij * mol_weight_ratio_ij**0.75
                           * ((1.0 + S_diag[:, None] / T_ref)
                              / (1.0 + S_diag[None, :] / T_ref)))
    A = 0.25 * term1**2 * ((1.0 + S[:] / T_ref)
                            / (1.0 + S_diag[:, None] / T_ref))

    np.fill_diagonal(A, 0.0)
    denom = (A @ fractions_array) / fractions_array  # element-wise, safe since we skip fractions_array[i]==0 species
    cond = np.sum(k_array / (1.0 + denom))  # need to rework slightly

    return float(specific_heat), float(viscosity), float(cond), float(density)
