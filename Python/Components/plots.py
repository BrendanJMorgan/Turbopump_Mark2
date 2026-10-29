import webbrowser
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from engine_state import engine, tca, gg
from powerhead import powerhead

PUMPS = (("ox", "Ox", "blue"), ("fuel", "Fuel", "red"))  # (powerhead attribute prefix, label, color)
ROW_HEIGHT = 400  # px per subplot row in the dashboard
INFO_MARGIN = 300  # px - right margin reserved for the info panel
NOZZLE_AXIAL_GAP = 0.001  # m - nozzle exit to blade leading edge; DRAWING ONLY, not in the turbine model yet


def plots():
    pumps = {key: getattr(powerhead, f"{key}_pump") for key, _, _ in PUMPS}
    turbine = powerhead.turbine

    figures = [
        tca_figure(),
        meridional_figure(pumps),
        overhead_figure(pumps, turbine),
        turbine_blade_figure(turbine),
        inducer_figure(pumps),
        inducer_unrolled_figure(pumps),
        volute_figure(pumps),
    ]

    for key, label, _ in PUMPS:
        pump = pumps[key]
        print(f"{label} pump: shaft power {pump.shaft_power/1e3:0.1f} kW, shaft speed {pump.shaft_speed*30/np.pi:0.0f} rpm")

    show_dashboard(figures, title="Mark 2 Engine")


#######################################################################################################
### Helpers
#######################################################################################################

def info_box(fig, lines, heading=None):
    """Queue lines for the figure's info panel. show_dashboard draws all of a figure's blocks as one panel in
    the right margin, under the legend, so nothing sits on top of the plot."""
    block = ([f"<b>{heading}</b>"] if heading else []) + list(lines)
    meta = dict(fig.layout.meta or {})
    meta["info"] = meta.get("info", []) + [block]
    fig.layout.meta = meta


def draw_info_panel(fig):
    """Render the queued info blocks (see info_box) as one annotation in the right margin."""
    blocks = (fig.layout.meta or {}).get("info", [])
    if not blocks:
        return
    fig.add_annotation(x=1.02, xanchor="left", y=0, yanchor="bottom", xref="paper", yref="paper",
                       showarrow=False, align="left", text="<br><br>".join("<br>".join(b) for b in blocks),
                       bgcolor="white", bordercolor="black", borderpad=4)
    fig.update_layout(margin_r=INFO_MARGIN)


def split_legends(fig):
    """Give each subplot row its own legend, placed beside that row (plotly >= 5.15 multi-legend).

    Traces are grouped by the vertical domain of their y axis; secondary y axes share their primary's row.
    3D traces (scenes) have no y axis and keep the default legend.
    """
    def domain_top(axis_ref):
        axis = fig.layout["yaxis" + axis_ref[1:]]
        if axis.overlaying: # secondary y axis: use the axis it sits on
            axis = fig.layout["yaxis" + axis.overlaying[1:]]
        return axis.domain[1] if axis.domain else 1.0

    rows = {}
    for trace in fig.data:
        if getattr(trace, "yaxis", None):
            rows.setdefault(domain_top(trace.yaxis), []).append(trace)
    if len(rows) < 2:
        return

    for k, top in enumerate(sorted(rows, reverse=True), start=1):
        name = "legend" if k == 1 else f"legend{k}"
        for trace in rows[top]:
            trace.legend = name
        fig.update_layout({name: dict(x=1.02, xanchor="left", y=top, yanchor="top")})


def equal_aspect(fig, row=None, col=None):
    """Lock 1:1 scale on a subplot (or the whole figure) and shrink its plot area to fit the data, so the
    locked box-zoom shape matches the geometry instead of a wide strip."""
    x_ref = fig.get_subplot(row, col).yaxis.anchor if row is not None else "x"  # the x axis this y axis sits on
    fig.update_yaxes(scaleanchor=x_ref, scaleratio=1, constrain="domain", row=row, col=col)
    fig.update_xaxes(constrain="domain", row=row, col=col)


def is_geometry(fig):
    """True if any axis in the figure is aspect-locked (see equal_aspect)."""
    return any(axis.scaleanchor for axis in fig.select_yaxes())


def rotate(x, y, angle):
    """Rotate points (x, y) about the origin by angle [rad]."""
    return x*np.cos(angle) - y*np.sin(angle), x*np.sin(angle) + y*np.cos(angle)


def line(fig, x, y, color, name, dash="solid", width=2, **subplot):
    """Line trace. Every trace is named; traces sharing a name share one legend entry (see dedupe_legend)."""
    fig.add_trace(go.Scatter(x=x, y=y, name=name, line=dict(color=color, dash=dash, width=width)), **subplot)


def dedupe_legend(fig):
    """One legend entry per name per legend: traces with the same name are grouped, so clicking the entry
    toggles all of them, and only the first shows in the legend."""
    seen = set()
    for trace in fig.data:
        if not trace.name:
            continue
        trace.legendgroup = trace.legendgroup or trace.name
        key = (trace.legend or "legend", trace.name)
        trace.showlegend = key not in seen
        seen.add(key)


#######################################################################################################
### TCA Contours + Coolant Properties
#######################################################################################################

def tca_figure():
    gas = tca.combustion_gas
    jacket = tca.regenerative_jacket
    cool = tca.regenerative_coolant
    film = tca.film_coolant

    rows = 4
    fig = make_subplots(rows=rows, cols=1, shared_xaxes=True,
                        specs=[[{}], [{"secondary_y": True}], [{}], [{}]], vertical_spacing=0.06)

    # Chamber wall: inner and outer contour, mirrored about the axis
    line(fig, tca.x, tca.r1, "brown", name="Thrust Chamber Contours", row=1, col=1)
    for y in (tca.r2, -tca.r1, -tca.r2):
        line(fig, tca.x, y, "brown", name="Thrust Chamber Contours", row=1, col=1)
    # No equal aspect here: an aspect lock on one row of a shared-x stack makes every row re-zoom together
    fig.update_yaxes(title_text="Radius (m)", row=1, col=1)

    # Regen coolant hydraulics
    line(fig, tca.x, cool.p / 1e5, "blue", name="Regen Pressure", row=2, col=1, secondary_y=False)
    line(fig, tca.x, cool.velocity, "red", name="Regen Velocity", row=2, col=1, secondary_y=True)
    fig.update_yaxes(title_text="Pressure (bar)", color="blue", row=2, col=1, secondary_y=False)
    fig.update_yaxes(title_text="Velocity (m/s)", color="red", row=2, col=1, secondary_y=True)

    # Temperatures through the wall stack
    temperatures = [(jacket.T_cold, "Cold Wall", "#1f77b4"), (jacket.T_hot, "Hot Wall", "#ff7f0e"),
                    (cool.T, "Regen Coolant", "#2ca02c"), (film.liquid.T, "Film Liquid", "#d62728"),
                    (film.gas.T, "Film Gas", "#9467bd"), (gas.T_free, "Free-stream Gas", "#8c564b"),
                    (gas.T_adiabatic_wall, "Adiabatic Wall Gas", "#e377c2"),
                    (gas.T_boundary, "Boundary Layer Reference", "#7f7f7f")]
    for y, name, color in temperatures:
        line(fig, tca.x, y, color, name=name, row=3, col=1)
    fig.update_yaxes(title_text="Temperature (K)", range=[0, 1.2 * tca.Tt], row=3, col=1)

    # Zoom on the wall and coolant temperatures, which the gas temperatures squash in the row above
    for y, name, color in temperatures[:5]:
        line(fig, tca.x, y, color, name=name, dash="dash" if name == "Film Gas" else "solid", row=4, col=1)
    fig.update_yaxes(title_text="Wall & Film T (K)", row=4, col=1)
    fig.update_xaxes(title_text="Distance from Injector (m)", row=4, col=1)

    fig.update_layout(title="Thrust Chamber Geometry and Regenerative Coolant Properties")
    return fig, rows


#######################################################################################################
### Meridional impeller contours + meanline parameters
#######################################################################################################

def meridional_figure(pumps):
    rows = 2
    fig = make_subplots(rows=rows, cols=1, specs=[[{}], [{"secondary_y": True}]],
                        subplot_titles=("Meridional Plane Impeller Contours (mm)", ""),
                        vertical_spacing=0.12)

    # Ox mirrored to negative r so both impellers share one plot
    for (key, label, color), side in zip(PUMPS, (-1.0, 1.0)):
        imp = pumps[key].impeller[0]
        for curve, part, dash in ((imp.shroud_curve, "Shroud", "solid"), (imp.hub_curve, "Hub", "dot"),
                                  (imp.blade.meanline_curve, "Blade Meanline", "dash")):
            line(fig, side * curve[:, 0] * 1000, curve[:, 1] * 1000, color, name=f"{label} {part}",
                 dash=dash, row=1, col=1)
    fig.add_vline(x=0, line=dict(color="green", width=1), name="Shaft Axis", showlegend=True, row=1, col=1)
    fig.add_hline(y=0, line=dict(color="green", width=1, dash="dot"), name="Inlet Plane", showlegend=True, row=1, col=1)
    fig.update_xaxes(title_text="r (mm)", row=1, col=1)
    fig.update_yaxes(title_text="z (mm)", row=1, col=1)
    equal_aspect(fig, row=1, col=1)

    for key, label, color in PUMPS:
        imp = pumps[key].impeller[0]
        line(fig, imp.arc_meanline * 1000, imp.area_meanline * 1e6, color, name=f"{label} Flow Area",
             row=2, col=1, secondary_y=False)
        line(fig, imp.arc_meanline * 1000, imp.blockage, color, name=f"{label} Blockage", dash="dash",
             row=2, col=1, secondary_y=True)
    fig.update_xaxes(title_text="Meanline Arc Length (mm)", row=2, col=1)
    fig.update_yaxes(title_text="Annular Flow Area (mm2)", row=2, col=1, secondary_y=False)
    fig.update_yaxes(title_text="Blockage (dashed)", row=2, col=1, secondary_y=True)

    fig.update_layout(title="Impeller Meridional & Meanline")
    return fig, rows


#######################################################################################################
### Impeller blades (plan view) + turbine circles
#######################################################################################################

def overhead_figure(pumps, turbine):
    fig = go.Figure()

    # Impeller blades and volute outline for each pump, overlaid on the shaft axis
    for key, label, color in PUMPS:
        pump = pumps[key]
        imp = pump.impeller[0]
        pitch_angle = pump.clocking * 2.0 * np.pi / imp.blade.count
        for i in range(imp.blade.count):
            x, y = rotate(imp.blade_curve[:, 0], imp.blade_curve[:, 1], i * pitch_angle)
            line(fig, x * 1000, y * 1000, color, name=f"{label} Impeller")

        # Volute in plan view: inner and outer edges of each cross section along the wrap
        if 'volute' in vars(pump):
            vol = pump.volute
            azimuth = pump.clocking * vol.theta
            for r in (vol.r_meanline - 0.5 * vol.d_hydraulic_section, vol.r_meanline + 0.5 * vol.d_hydraulic_section):
                line(fig, r * 1000 * np.cos(azimuth), r * 1000 * np.sin(azimuth), color, name=f"{label} Volute",
                     dash="dot", width=1.5)
        else:
            # Vaned diffuser: leading and trailing edge circles
            dif = pump.vaned_diffuser
            circle = np.linspace(0.0, 2.0 * np.pi, 200)
            for r in (dif.r_inlet, dif.r_outlet):
                line(fig, r * 1000 * np.cos(circle), r * 1000 * np.sin(circle), color, name=f"{label} Diffuser",
                     dash="dot", width=1.5)

    # Turbine rotor: pitchline (dashed), blade tip and base
    theta = np.linspace(0.0, 2.0 * np.pi, 1000)
    for r, name, dash in ((turbine.r_pitchline, "Turbine Pitchline", "dash"),
                          (turbine.r_tip, "Turbine Blade Tip/Base", "solid"),
                          (turbine.r_base, "Turbine Blade Tip/Base", "solid")):
        line(fig, r * 1000 * np.cos(theta), r * 1000 * np.sin(theta), "orange", name=name, dash=dash, width=1)

    info_box(fig, [
        f"Shaft Speed = {turbine.shaft_speed*30/np.pi:0.2f} rpm",
        f"GG = {gg.mdot:0.2g} kg/s, {gg.mdot*1790:0.0f} SCFM N2 ({gg.mdot/engine.mdot*100:.2g}%)",
    ])
    equal_aspect(fig)
    fig.update_layout(title="Impellers, Volutes, and Turbine Rotor (mm)")
    return fig, 1


#######################################################################################################
### Turbine blade profiles (2D cascade)
#######################################################################################################

def turbine_blade_figure(turbine):
    fig = go.Figure()
    chord, pitch, r_le = turbine.blade_chord, turbine.blade_pitch, turbine.radius_leading

    # Trailing (pressure-side) arc through the blade tips
    sagitta = chord / 2.0 * (1.0 - np.cos(turbine.blade_angle)) / np.sin(turbine.blade_angle)
    tip = np.array([-sagitta, chord / 2.0])
    blade_radius = sagitta / 2.0 + chord**2 / (8.0 * sagitta)
    center_arc = np.array([-blade_radius, 0.0])
    angle_trailing = np.arctan2(tip[1], tip[0] + blade_radius)
    arc_theta = np.linspace(-angle_trailing, angle_trailing, 120)
    x_trailing = center_arc[0] + blade_radius * np.cos(arc_theta)
    y_trailing = center_arc[1] + blade_radius * np.sin(arc_theta)

    # Leading-edge fillet, joined to the tips by tangent lines
    apex_leading = np.array([pitch - turbine.blade_opening, 0.0])
    center_leading = apex_leading + np.array([-r_le, 0.0])
    center_to_tip = tip - center_leading
    tangent_angle = np.arctan2(center_to_tip[1], center_to_tip[0]) - np.arccos(r_le / np.linalg.norm(center_to_tip))
    tangency_leading = center_leading + r_le * np.array([np.cos(tangent_angle), np.sin(tangent_angle)])
    theta_leading = np.linspace(-tangent_angle, tangent_angle, 60)
    x_leading = center_leading[0] + r_le * np.cos(theta_leading)
    y_leading = center_leading[1] + r_le * np.sin(theta_leading)

    # Three neighbouring blades
    for offset in (-pitch * 1000, 0.0, pitch * 1000):
        line(fig, x_trailing * 1000 + offset, y_trailing * 1000, "gold", name="Rotor Blades")
        line(fig, x_leading * 1000 + offset, y_leading * 1000, "gold", name="Rotor Blades")
        for side in (1.0, -1.0):
            line(fig, [tangency_leading[0] * 1000 + offset, tip[0] * 1000 + offset],
                 [side * tangency_leading[1] * 1000, side * tip[1] * 1000], "gold", name="Rotor Blades")

    # Nozzle plate above the blades: section through each nozzle axis in the cascade plane (x tangential, y axial).
    # Jets point down into the blades and along +x, onto the concave faces; tilt is angle_nozzle from axial.
    r_throat = turbine.d_throat_nozzle / 2
    r_exit = r_throat * np.sqrt(gg.Ae_At)                                  # same expansion ratio as the GG nozzle
    length_diverging = (r_exit - r_throat) / np.tan(turbine.diverge_angle_gg)
    length_converging = 2 * r_throat                                       # DRAWING ONLY - short inlet stub, not modeled
    along = np.array([np.sin(turbine.angle_nozzle), -np.cos(turbine.angle_nozzle)])   # unit vector along the jet
    across = np.array([along[1], -along[0]])                                           # unit normal to the jet

    n_nozzles = max(int(round(turbine.nozzle_number)), 1)
    nozzle_spacing = 2 * np.pi * turbine.r_pitchline * turbine.admission_fraction / n_nozzles  # m - along the pitchline
    plate_face_y = chord / 2 + NOZZLE_AXIAL_GAP  # m - bottom face of the nozzle plate, where the bores break out

    def to_plate_face(point, direction):
        """Extend a line from point along direction until it meets the plate's bottom face."""
        return point + (plate_face_y - point[1]) / direction[1] * direction

    for k in range(-1, 2):
        axis_exit = np.array([k * nozzle_spacing, plate_face_y])          # nozzle axis meets the plate face here
        throat = axis_exit - length_diverging * along                     # throat sits one cone length up the axis
        inlet = throat - length_converging * along
        for side in (1.0, -1.0):
            throat_wall = throat + side * r_throat * across
            cone = rotate(along[0], along[1], -side * turbine.diverge_angle_gg)  # wall direction, opening away from the axis
            wall = np.array([inlet + side * 2 * r_throat * across, throat_wall, to_plate_face(throat_wall, np.array(cone))])
            line(fig, wall[:, 0] * 1000, wall[:, 1] * 1000, "gray", name="Nozzles")
        axis_line = np.array([inlet, axis_exit])
        line(fig, axis_line[:, 0] * 1000, axis_line[:, 1] * 1000, "gray", name="Nozzle Axes", dash="dashdot", width=1)
    fig.add_hline(y=plate_face_y * 1000, line=dict(color="gray", width=1, dash="dot"), name="Nozzle Plate Face",
                  showlegend=True)

    info_box(fig, [f"Blade Count = {turbine.blade_count:g}",
                             f"Nozzles = {turbine.nozzle_number:0.1f}, spacing {nozzle_spacing*1000:0.1f} mm",
                             f"Nozzle tilt = {np.degrees(turbine.angle_nozzle):0.1f} deg from axial",
                             f"Axial gap = {NOZZLE_AXIAL_GAP*1000:0.1f} mm (assumed)"])
    equal_aspect(fig)
    fig.update_xaxes(title_text="Tangential (mm), blades move +x")
    fig.update_yaxes(title_text="Axial (mm)")
    fig.update_layout(title="Nozzles and Rotor Blade Geometry (mm)")
    return fig, 1


#######################################################################################################
### Inducer Geometry (3D)
#######################################################################################################

def inducer_figure(pumps):
    fig = make_subplots(rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "scene"}]],
                        subplot_titles=[f"{label} Inducer Geometry (mm)" for _, label, _ in PUMPS])

    for col, (key, label, color) in enumerate(PUMPS, start=1):
        pump = pumps[key]
        ind = pump.inducer
        for k in range(ind.blade_number):
            x, y = rotate(ind.x_bladeline, ind.y_bladeline, 2 * np.pi * k / ind.blade_number)
            z = ind.z_bladeline
            # Bladelines (hub to tip), then leading and trailing edges drawn heavier
            for n in range(x.shape[0]):
                fig.add_trace(go.Scatter3d(x=x[n] * 1000, y=y[n] * 1000, z=z[n] * 1000, mode="lines",
                                           name=f"{label} Bladelines", line=dict(color=color, width=2)), row=1, col=col)
            for edge, edge_name in ((0, "Leading Edges"), (-1, "Trailing Edges")):
                fig.add_trace(go.Scatter3d(x=x[:, edge] * 1000, y=y[:, edge] * 1000, z=z[:, edge] * 1000,
                                           mode="lines", name=f"{label} {edge_name}",
                                           line=dict(color=color, width=4, dash="solid" if edge == 0 else "dash")),
                              row=1, col=col)

        info_box(fig, heading=f"{label} Inducer", lines=[
            f"phi = {ind.flow_coeff:0.4g}",
            f"psi = {ind.head_coeff:0.3g}",
            f"Nss = {2733.0*ind.suction_specific_speed_available:0.0f}",  # Omega_ss -> US Nss (rpm, gpm, ft)
            f"Shaft = {30/np.pi*pump.shaft_speed:0.0f} RPM",
        ])

    fig.update_layout(scene=dict(aspectmode="data"), scene2=dict(aspectmode="data"), title="Inducer 3D Geometry")
    return fig, 1


#######################################################################################################
### Inducer unrolled bladelines + blade angle distribution
#######################################################################################################

def inducer_unrolled_figure(pumps):
    rows = 2
    fig = make_subplots(rows=rows, cols=2, vertical_spacing=0.12,
                        subplot_titles=[f"{label} Inducer - Unrolled Bladelines (light = hub, dark = tip)" for _, label, _ in PUMPS]
                                       + [f"{label} Inducer - Blade Angle" for _, label, _ in PUMPS])

    for col, (key, label, color) in enumerate(PUMPS, start=1):
        ind = pumps[key].inducer
        n_lines = ind.r_bladeline.shape[0]
        s = np.linspace(0, 1, ind.blade_angle.shape[1])
        for n in range(n_lines):
            shade = dict(color=color, width=1)
            opacity = 0.3 + 0.7 * n / (n_lines - 1)  # hub light -> tip dark
            fig.add_trace(go.Scatter(x=ind.r_bladeline[n] * ind.theta_bladeline[n] * 1000, y=ind.z_bladeline[n] * 1000,
                                     name=f"{label} Bladelines", line=shade, opacity=opacity,
                                     hovertext=f"bladeline {n}"), row=1, col=col)
            fig.add_trace(go.Scatter(x=s, y=np.rad2deg(ind.blade_angle[n]), name=f"{label} Blade Angle",
                                     line=shade, opacity=opacity, hovertext=f"bladeline {n}"), row=2, col=col)

        fig.update_xaxes(title_text="r*theta (mm)", row=1, col=col)
        fig.update_yaxes(title_text="z (mm)", row=1, col=col)
        equal_aspect(fig, row=1, col=col)
        fig.update_xaxes(title_text="s / L", row=2, col=col)
        fig.update_yaxes(title_text="beta (deg)", row=2, col=col)

    fig.update_layout(title="Inducer Blade Development")
    return fig, rows


#######################################################################################################
### Volute Geometry (3D)
#######################################################################################################

def volute_figure(pumps):
    fig = make_subplots(rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "scene"}]],
                        subplot_titles=[f"{label} Volute Geometry (mm)" for _, label, _ in PUMPS])

    def surface(x, y, z, color, opacity, col, name):
        fig.add_trace(go.Surface(x=x * 1000, y=y * 1000, z=z * 1000, opacity=opacity, showscale=False, name=name,
                                 showlegend=True, colorscale=[[0, color], [1, color]]), row=1, col=col)

    theta = np.linspace(0.0, 2.0 * np.pi, 64)
    for col, (key, label, color) in enumerate(PUMPS, start=1):
        pump = pumps[key]
        imp = pump.impeller[0]
        if 'volute' in vars(pump):
            vol = pump.volute
            surface(vol.contour[:, :, 0], vol.contour[:, :, 1], vol.contour[:, :, 2], color, 0.6, col, f"{label} Volute")

        # Impeller shroud and hub as surfaces of revolution, for reference
        for curve, shade, part in ((imp.shroud_curve, "gray", "Shroud"), (imp.hub_curve, "lightgray", "Hub")):
            r, z = curve[:, 0][:, None], curve[:, 1][:, None]
            x, y = r * np.cos(theta)[None, :], r * np.sin(theta)[None, :]
            surface(x, y, np.broadcast_to(z, x.shape), shade, 0.25, col, f"{label} Impeller {part}")

        if 'volute' in vars(pump):
            info_box(fig, heading=f"{label} Volute", lines=[
                f"method = {vol.design_method}",
                f"d_throat = {vol.d_throat*1000:0.2f} mm",
                f"c_throat = {vol.c_throat:0.1f} m/s",
                f"loss = {vol.total_loss:0.4g}",
            ])
        else:
            dif = pump.vaned_diffuser
            info_box(fig, heading=f"{label} Vaned Diffuser", lines=[
                f"vanes = {dif.vane_count}",
                f"r3 = {dif.r_inlet*1000:0.2f} mm, r4 = {dif.r_outlet*1000:0.2f} mm",
                f"a3 = {dif.throat_width_inlet*1000:0.2f} mm",
                f"blade angle = {np.degrees(dif.blade_angle_inlet):0.1f} deg",
            ])

    fig.update_layout(scene=dict(aspectmode="data"), scene2=dict(aspectmode="data"), title="Volute 3D Geometry")
    return fig, 1


#######################################################################################################
### Dashboard
#######################################################################################################

def show_dashboard(figures, title="Mark 2 Plots", open_browser=True):
    """Stitch (figure, n_rows) pairs into one scrollable HTML page where every plot stays fully
    interactive (3D rotate, hover, zoom). plotly.js is embedded once, so the page works offline.
    """
    sections = []
    for i, (fig, rows) in enumerate(figures):
        split_legends(fig)
        dedupe_legend(fig)
        heading = fig.layout.title.text or f"Figure {i + 1}"  # becomes the collapsible section's header
        fig.update_layout(title=None, template="plotly_white", autosize=True, margin=dict(l=10, r=10, t=30, b=10))
        fig.update_xaxes(showline=True, linecolor="#bbb", gridcolor="#e6e6e6", zerolinecolor="#cfcfcf")
        fig.update_yaxes(showline=True, linecolor="#bbb", gridcolor="#e6e6e6", zerolinecolor="#cfcfcf")
        draw_info_panel(fig)
        plot = fig.to_html(
            full_html=False,
            include_plotlyjs=(i == 0),  # embed plotly.js with the first figure only
            default_width="96%",
            default_height=ROW_HEIGHT * rows,
            config={"responsive": True, "scrollZoom": is_geometry(fig)},  # wheel zoom only on equal-aspect plots
        )
        sections.append(f"<details open><summary>{heading}</summary>{plot}</details>")

    body = "\n".join(sections)
    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>{title}</title>
<style>
  body {{ margin:0; background:#f3f3f3;
         font-family:-apple-system,Segoe UI,Roboto,sans-serif; }}
  h1 {{ padding:14px 20px; margin:0; font-size:18px;
        background:#1c1c1c; color:#eee; position:sticky; top:0; z-index:10; }}
  details {{ margin:14px auto; padding:8px; max-width:1200px;
             background:#fff; border-radius:8px;
             box-shadow:0 1px 4px rgba(0,0,0,.12); }}
  summary {{ cursor:pointer; padding:6px 10px; font-size:16px; color:#333; }}
  summary:hover {{ color:#000; }}
  .controls {{ max-width:1200px; margin:10px auto 0; text-align:right; }}
  .controls button {{ margin-left:6px; }}
</style></head><body>
<h1>{title}</h1>
<div class="controls"><button onclick="setAll(true)">Expand all</button><button onclick="setAll(false)">Collapse all</button></div>
{body}
<script>
  // Plots in a collapsed section have no size; re-fit them when their section opens.
  document.querySelectorAll("details").forEach(d => d.addEventListener("toggle", () => {{
    if (d.open) d.querySelectorAll(".plotly-graph-div").forEach(p => Plotly.Plots.resize(p));
  }}));
  function setAll(open) {{ document.querySelectorAll("details").forEach(d => d.open = open); }}
</script>
</body></html>"""

    out = Path(__file__).resolve().parent / "mark2_plots.html"  # next to this script; ignored by git
    out.write_text(html, encoding="utf-8")
    if open_browser:
        webbrowser.open(out.as_uri())
    return out


if __name__ == "__main__":
    plots()
