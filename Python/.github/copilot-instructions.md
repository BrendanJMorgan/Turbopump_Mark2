# Copilot Instructions for Turbopump_Mark2/Python

## Project Overview
This codebase models a liquid rocket engine turbopump system, including combustion chamber, pumps, turbines, and supporting thermofluidic calculations. It is organized as a set of engineering scripts, each focused on a subsystem or analysis step. The main entry point is `main.py`, which defines the `Inputs` dataclass for all key parameters and orchestrates subsystem calculations.

## Major Components
- **main.py**: Central configuration and workflow. Defines all global parameters in `Inputs`.
- **combustion_chamber.py, cc_geometry.py, cc_gas_flow.py, thermal_balance.py**: Combustion chamber and nozzle geometry, gas flow, and heat transfer.
- **pump.py, impeller.py, inducer.py, blades.py, volute.py**: Pump and turbopump geometry, performance, and blade design.
- **powerhead.py, gas_generator.py, turbine.py, turbine_optimization.py**: Powerhead, gas generator, and turbine calculations.
- **plots.py**: Visualization utilities.

## Key Patterns & Conventions
- **Parameter Passing**: All major calculations use parameters from the `Inputs` dataclass in `main.py`. Update or extend this class to add new global parameters.
- **Script Structure**: Each subsystem script is self-contained, often with a MATLAB-like structure (section headers, vectorized calculations, and inline documentation).
- **Units**: SI units are used throughout unless otherwise noted in comments.
- **Data Files**: `.npz` files (e.g., `expansion_loss_coeff.npz`) store empirical or precomputed data for loss models and are loaded as needed.
- **External Dependencies**: Requires `numpy`, `CoolProp`, and `rocketprops` Python packages. Some scripts reference MATLAB/Octave syntax or legacy code—port as needed.

## Developer Workflows
- **Run Analysis**: Execute `main.py` to run a full system analysis. All parameters are set in the `Inputs` dataclass.
- **Add/Modify Subsystems**: Create or edit subsystem scripts. Reference and update `Inputs` for new parameters.
- **Data Flow**: Most scripts expect to be run in the context of `main.py` or with all required variables defined in the global scope.
- **Debugging**: Use print statements or add plotting code in `plots.py` for diagnostics. There is no formal test suite or build system.

## Integration & Cross-Component Patterns
- **Shared Variables**: Many scripts rely on variables defined in `main.py` or other scripts. Maintain consistency in naming and units.
- **Empirical Models**: Losses, efficiencies, and performance factors are often interpolated from `.npz` data or hardcoded arrays.
- **Output**: Some scripts write results to text files (e.g., `Output Parameters/`, `Curves/`).

## Example: Adding a New Pump Model
1. Add new parameters to `Inputs` in `main.py`.
2. Create a new script (e.g., `new_pump.py`) following the structure of `pump.py`.
3. Reference new parameters from `Inputs` and ensure units match existing conventions.
4. Update `main.py` to call the new script as needed.

## Tips for AI Agents
- Always check `main.py` for parameter definitions and workflow logic.
- Maintain SI units and vectorized calculations.
- When in doubt, follow the structure and conventions of existing subsystem scripts.
- Document any new empirical data or models added to `.npz` files.
