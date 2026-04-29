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
    # powerhead() # contains submodules pump, gas_generator, turbine
    plots()

if __name__ == "__main__":
    main()
