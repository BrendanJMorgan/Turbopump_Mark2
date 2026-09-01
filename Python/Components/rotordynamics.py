import numpy as np
import ross as rs
#######################################################################################################
### Shaft
#######################################################################################################
shaft_length = 1.5     # m   - total shaft length (long shaft through bottom tank)
d_outer_shaft     = 0.0254*1.5     # m   - 10 mm shaft (matches Mark 1 stock)
d_inner_shaft     = 0.0254*1.4       # m   - solid shaft
shaft_material = rs.Material(name="SS303", rho=8000, E=1.93e11, Poisson=0.30)

#######################################################################################################
### Components
#######################################################################################################
#   name              pos[m]   o_d[m]   i_d[m]   width[m]
components = [
    ("ox_impeller",   0.060,   0.0678,  d_outer_shaft, 0.020),   # rexit ~33.9 mm -> od ~67.8 mm
    ("ox_inducer",    0.090,   0.0268,  d_outer_shaft, 0.015),   # rtip ~13.4 mm
    ("turbine_wheel", shaft_length-0.1,   0.1600,  d_outer_shaft, 0.012),   # r_pitchline 0.08 m -> od ~160 mm
    ("fuel_impeller", shaft_length-0.04,   0.0776,  d_outer_shaft, 0.020),   # rexit ~38.8 mm
    ("fuel_inducer",  shaft_length,   0.0330,  d_outer_shaft, 0.015),   # rtip ~16.5 mm
]

bearing_stiffness = 1.0e7    # N/m  - placeholder radial stiffness (TUNE to real bearings)
bearing_damping = 0.0      # N*s/m- placeholder damping

#   pos[m]   kxx   kyy   cxx  cyy
bearings = [
    (0.000,  bearing_stiffness,    bearing_stiffness,    bearing_damping, bearing_damping),   # top bearing, above ox impeller
    (0.300,  bearing_stiffness,    bearing_stiffness,    bearing_damping, bearing_damping),   # below ox section
    (0.600,  bearing_stiffness,    bearing_stiffness,    bearing_damping, bearing_damping),   # mid
    (0.900,  bearing_stiffness,    bearing_stiffness,    bearing_damping, bearing_damping),   # above turbine
    (shaft_length-0.08,  bearing_stiffness,    bearing_stiffness,    bearing_damping, bearing_damping),   # below turbine
]

#######################################################################################################
### Meshing + Construction
#######################################################################################################
positions = sorted(
    {0.0, shaft_length}
    | {p for _, p, *_ in components}
    | {b[0] for b in bearings}
)

# Check nothing hanging off the end of the shaft
if positions[0] < -1e-12 or positions[-1] > shaft_length + 1e-12:
    raise ValueError("A component or bearing lies outside the shaft length")

# One shaft element per gap between consecutive nodes
shaft_elements = []
for i in range(len(positions) - 1):
    L = positions[i + 1] - positions[i]
    if L <= 0:
        raise ValueError(f"Zero/negative shaft element length at {positions[i]} m "
                         "(two components share a position - merge or nudge them)")
    shaft_elements.append(
        rs.ShaftElement(
            L=L,
            idl=d_inner_shaft,
            odl=d_outer_shaft,
            material=shaft_material,
            shear_effects=True,
            rotary_inertia=True,
            gyroscopic=True,
        )
    )

def node_at(pos):
    """Return the node index sitting exactly at axial position `pos`."""
    idx = np.argmin(np.abs(np.array(positions) - pos))
    if abs(positions[idx] - pos) > 1e-9:
        raise ValueError(f"No node at {pos} m - mesh build failed")
    return int(idx)

disk_elements = [
    rs.DiskElement.from_geometry(
        n=node_at(pos), material=shaft_material, width=w, i_d=i_d, o_d=o_d
    )
    for name, pos, o_d, i_d, w in components
]

bearing_elements = [
    rs.BearingElement(n=node_at(pos), kxx=kxx, kyy=kyy, cxx=cxx, cyy=cyy)
    for pos, kxx, kyy, cxx, cyy in bearings
]

rotor = rs.Rotor(shaft_elements, disk_elements, bearing_elements)

print(f"Nodes: {len(positions)}   Shaft elements: {len(shaft_elements)}")
print(f"Total rotor mass:  {rotor.m:.3f} kg")
print(f"Center of gravity: {rotor.CG:.4f} m from top end")
for name, pos, *_ in components:
    print(f"  {name:<14} node {node_at(pos):>2}  @ {pos*1000:6.1f} mm")
for pos, *_ in bearings:
    print(f"  bearing        node {node_at(pos):>2}  @ {pos*1000:6.1f} mm")

#######################################################################################################
### Analysis
#######################################################################################################
shaft_speed = np.pi/30 * 20000

# Geometry plot
rotor.plot_rotor().show()

# Undamped critical speed map - where do bending criticals sit vs 20k rpm?
ucs = rotor.run_ucs(stiffness_range=(5, 9), num=30, num_modes=8)
ucs.plot().show()

# Campbell diagram across a speed sweep
speeds = np.linspace(0, 1.3 * shaft_speed, 40)
campbell = rotor.run_campbell(speeds)
campbell.plot().show()

# Modal analysis at operating speed -> natural frequencies & damping
modal = rotor.run_modal(speed=shaft_speed)
modal.plot_mode_2d(0).show()
print("\nNatural frequencies at "
      f"{shaft_speed*30/np.pi:.0f} rpm (Hz):")
for i, wn in enumerate(modal.wn[:6]):
    print(f"  mode {i}: {wn/(2*np.pi):8.1f} Hz   "
          f"({wn/(2*np.pi)*60:8.0f} rpm)   "
          f"zeta = {modal.damping_ratio[i]:.4f}")

print(f"\nOperating speed: {shaft_speed} rpm = {shaft_speed/(2*np.pi):.1f} Hz")