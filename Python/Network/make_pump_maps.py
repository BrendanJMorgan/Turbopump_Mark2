"""
Regenerate the placeholder pump maps into a FullFlow v2 HDF5 file.

The cycle YAML keeps scalars inline but pushes the large (Q, N) grids out to
``pump_maps.h5``.  This script builds that file from a single affinity-scaled
design point per pump - the same synthetic model that used to live in
``mark2_design_point.placeholder_map``.  It is a PLACEHOLDER: replace each group
with a real sweep of ``impellers.py`` over a (volumetric_flow, rotor_speed) grid.

    python make_pump_maps.py            # writes ./pump_maps.h5

Group layout (one per pump), readable by ``fullflow.Map.from_hdf5``:

    /<pump>/
        attrs["axis_order"] = ["volumetric_flow", "rotor_speed"]
        axes/volumetric_flow   (nq,)   m^3/s, ascending
        axes/rotor_speed       (nn,)   rpm,   ascending
        outputs/head_rise      (nq, nn) m
        outputs/efficiency     (nq, nn) -
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np

# ---------------------------------------------------------------------------
# design points (PLACEHOLDER - mirror mark2_design_point constants)
# ---------------------------------------------------------------------------
DESIGN_SPEED = 30_000.0          # rpm

PUMPS = {
    "oxidizer": dict(
        mass_flow=10.885,        # kg/s
        density=1142.0,          # kg/m^3
        head_rise=531.0,         # m
        efficiency=0.65,
    ),
    "fuel": dict(
        mass_flow=5.575,
        density=810.0,
        head_rise=995.0,
        efficiency=0.60,
    ),
}


def affinity_map(design_flow: float, design_head: float,
                 design_efficiency: float, design_speed: float) -> dict:
    """
    Affinity-scaled synthetic map.  PLACEHOLDER.

        H(Q, N)   = H_d * (N/N_d)^2 * psi(phi)
        eta(Q, N) = eta_d * (1 - 1.6 * (phi - 1)^2)
        phi       = (Q / Q_d) / (N / N_d)
    """
    flow = np.linspace(0.15 * design_flow, 1.9 * design_flow, 25)
    speed = np.linspace(0.25 * design_speed, 1.3 * design_speed, 15)

    phi = (flow[:, None] / design_flow) / (speed[None, :] / design_speed)
    speed_ratio = (speed[None, :] / design_speed) ** 2

    head = design_head * speed_ratio * (1.15 - 0.15 * phi - 0.20 * phi**2)
    efficiency = design_efficiency * (1.0 - 1.6 * (phi - 1.0) ** 2)
    efficiency = np.clip(efficiency, 0.05, 0.95)

    return dict(volumetric_flow=flow, rotor_speed=speed,
                head_rise=head, efficiency=efficiency)


def write_group(parent: h5py.Group, name: str, table: dict) -> None:
    group = parent.create_group(name)
    group.attrs["axis_order"] = json.dumps(["volumetric_flow", "rotor_speed"])

    axes = group.create_group("axes")
    axes.create_dataset("volumetric_flow", data=table["volumetric_flow"])
    axes.create_dataset("rotor_speed", data=table["rotor_speed"])
    for axis in axes.values():
        axis.attrs["spacing"] = "linear"

    outputs = group.create_group("outputs")
    outputs.create_dataset("head_rise", data=table["head_rise"])
    outputs.create_dataset("efficiency", data=table["efficiency"])


def main() -> None:
    out = Path(__file__).resolve().parent / "pump_maps.h5"
    with h5py.File(out, "w") as file:
        for pump, dp in PUMPS.items():
            table = affinity_map(
                dp["mass_flow"] / dp["density"],
                dp["head_rise"],
                dp["efficiency"],
                DESIGN_SPEED,
            )
            write_group(file, pump, table)
    print(f"wrote {out}  ({', '.join(PUMPS)})")


if __name__ == "__main__":
    main()
