"""STEP export of pump flow-passage geometry for SolidWorks.

Each volute is exported as one solid: the fluid "core" running from just past the tongue, around the wrap, and out
through the discharge cone. It is lofted through closed section loops sampled straight from vol.contour and
vol.nozzle_contour, so any section shape exports the same way. In SolidWorks, subtract it from a casing body
(Combine > Subtract) or shell it. Insert it with 3D Interconnect and right-click > Update after re-running Python.
"""
from pathlib import Path

import numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCP.GeomAPI import GeomAPI_Interpolate
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_ReturnStatus
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer
from OCP.collections import HArray1_gp_Pnt
from OCP.gp import gp_Pnt

CAD_EXPORT_DIR = Path(__file__).resolve().parents[2] / "Solidworks" / "Volute"  # STEP files land here
VOLUTE_EXPORT_THETA_START = np.radians(10.0)  # rad - solid starts past the tongue (section collapses to a point at theta = 0)
VOLUTE_LOFT_SECTIONS = 24  # section loops along the wrap - more follows the area schedule closer but lofts slower
NOZZLE_LOFT_SECTIONS = 4   # section loops along the discharge cone, after the throat


def export_volute_step(pumps):
    """Write <pump>_volute.step for every pump that has a volute. pumps: {key: pump}."""
    CAD_EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    for key, pump in pumps.items():
        if 'volute' not in vars(pump):
            continue
        vol = pump.volute
        solid = volute_core_solid(vol)
        path = CAD_EXPORT_DIR / f"{key}_volute.step"
        write_step(solid, path)
        print(f"{key} volute core exported to {path} ({solid_volume(solid)/1e3:0.2f} cm3)")


def volute_core_solid(vol):
    """Fluid-passage solid in mm: the volute wrap and the discharge cone, lofted separately and fused at the throat.
    One loft through both overshoots where the curved wrap meets the straight cone (+4 % volume); the cone
    leaves the throat tangent to the wrap, so the joint is smooth without it."""
    volute_grid = vol.contour[vol.theta >= VOLUTE_EXPORT_THETA_START]  # [m,m,m] - (n_theta, n_phi, 3)
    i_volute = np.linspace(0, len(volute_grid) - 1, VOLUTE_LOFT_SECTIONS).round().astype(int)
    i_nozzle = np.linspace(0, len(vol.nozzle_contour) - 1, NOZZLE_LOFT_SECTIONS + 1).round().astype(int)

    wrap = loft_solid([volute_grid[i] for i in i_volute], ruled=False)  # smooth B-spline skin
    cone = loft_solid([vol.nozzle_contour[i] for i in i_nozzle], ruled=True)  # straight walls, exact for a cone

    fuse = BRepAlgoAPI_Fuse(wrap, cone)
    fuse.Build()
    unify = ShapeUpgrade_UnifySameDomain(fuse.Shape())  # merge the shared throat faces away
    unify.Build()
    solid = unify.Shape()
    if not BRepCheck_Analyzer(solid).IsValid():
        raise RuntimeError("Volute core solid is invalid")
    return solid


def loft_solid(loops, ruled):
    """Capped solid lofted through closed section loops (each (n, 3) in m)."""
    loft = BRepOffsetAPI_ThruSections(True, ruled)
    loft.CheckCompatibility(False)  # every loop starts at phi = 0 and runs the same way - don't let OCC re-clock them
    for loop in loops:
        loft.AddWire(closed_spline_wire(loop * 1000))  # mm
    loft.Build()
    if not loft.IsDone():
        raise RuntimeError("Loft failed")
    return loft.Shape()


def closed_spline_wire(points):
    """Periodic B-spline wire through a closed loop of points (n, 3). A repeated closing point is dropped."""
    if np.allclose(points[0], points[-1]):
        points = points[:-1]
    array = HArray1_gp_Pnt(1, len(points))
    for i, (x, y, z) in enumerate(points, start=1):
        array.SetValue(i, gp_Pnt(float(x), float(y), float(z)))
    spline = GeomAPI_Interpolate(array, True, 1e-6)  # periodic
    spline.Perform()
    edge = BRepBuilderAPI_MakeEdge(spline.Curve()).Edge()
    return BRepBuilderAPI_MakeWire(edge).Wire()


def write_step(shape, path):
    writer = STEPControl_Writer()
    Interface_Static.SetCVal_s("write.step.unit", "MM")
    writer.Transfer(shape, STEPControl_AsIs)
    if writer.Write(str(path)) != IFSelect_ReturnStatus.IFSelect_RetDone:
        raise RuntimeError(f"STEP write failed: {path}")


def solid_volume(shape):
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, props)
    return props.Mass()  # mm3
