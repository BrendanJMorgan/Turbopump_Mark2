import os, sys, subprocess
from pathlib import Path

# _vpy = Path(__file__).resolve().parent / ".venv" / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
# if _vpy and Path(sys.executable).resolve() != _vpy.resolve():
#     env = os.environ.copy()
#     env["PATH"] = str(_vpy.parent) + os.pathsep + env.get("PATH", "")
#     subprocess.Popen([str(_vpy), *sys.argv], env=env)  # spawn venv python
#     sys.exit(0)

import numpy as np
import CoolProp as cp
import rocketprops as rp
from dataclasses import dataclass, asdict

from thrust_chamber_assembly import thrust_chamber_assembly
from powerhead import powerhead
from engine_state import engine, tca, pump, turbine, gg
# from plots import plots

def main():

    thrust_chamber_assembly() # contains submodules combustion_chamber, cc_geometry, cc_gas_flow, coolant_flow
    powerhead() # contains submodules pump, gas_generator, turbine
    isp_real = engine.thrust/(tca.mdot+gg.mdot)/engine.g
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

    # # Plots
    # plots(inputs, tc_results, ph_results)

if __name__ == "__main__":
    main()
