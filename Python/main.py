# %%

import numpy as np
import CoolProp as cp
import rocketprops as rp
from dataclasses import dataclass, asdict

from thrust_chamber_assembly import thrust_chamber_assembly
from powerhead import powerhead
from plots import plots
from engine_state import engine, tca, pump, turbine, gg


def main():
    thrust_chamber_assembly() # contains submodules combustion_chamber, cc_geometry, cc_gas_flow, coolant_flow
    powerhead() # contains submodules pump, gas_generator, turbine
    isp_real = engine.thrust/(tca.mdot+gg.mdot)/engine.g
    plots()

    # # Results (match MATLAB names)
    # thrust = ph_results.get("thrust", np.nan) 
    # isp_ideal = ph_results.get("isp_ideal", np.nan)
    # isp_real = ph_results.get("isp_real", np.nan)
    # mdot_gg = ph_results.get("mdot_gg", np.nan)
    # mdot_total = ph_results.get("mdot_total", np.nan)

    # # Print results section 
    # print("\n== Results ==")
    # print(f"thrust [N]: {thrust}")
    # print(f"isp_ideal [s]: {isp_ideal}")
    # print(f"isp_real [s]: {isp_real}")
    # print(f"mdot_gg [kg/s]: {mdot_gg}")
    # print(f"mdot_total [kg/s]: {mdot_total}")

if __name__ == "__main__":
    main()
