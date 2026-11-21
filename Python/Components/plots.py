import webbrowser
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from engine_state import engine, tca, gg
from powerhead import powerhead
from cad_export import export_volute_step

PUMPS = (("ox", "Ox", "blue"), ("fuel", "Fuel", "red"))  # (powerhead attribute prefix, label, color)
ROW_HEIGHT = 400  # px per subplot row in the dashboard
INFO_MARGIN = 300  # px - right margin reserved for the info panel
NOZZLE_AXIAL_GAP = 0.001  # m - nozzle exit to blade leading edge; DRAWING ONLY, not in the turbine model yet


def plots():
    pumps = {key: getattr(powerhead, f"{key}_pump") for key, _, _ in PUMPS}
    turbine = powerhead.turbine

    # Tabs of figures; a (name, [figures]) entry is a collapsible subsection, a string is raw HTML.
    # Geometry leads each section.
    tabs = {
        "System":  [summary_html(pumps, turbine, powerhead.link), overhead_figure(pumps, turbine)],
        "TCA":     [tca_figure()],
        "Pumps":   [("Impeller", [meridional_figure(pumps)]),
                    ("Inducer", [inducer_figure(pumps), inducer_unrolled_figure(pumps)]),
                    ("Volute / Diffuser", [volute_figure(pumps)])],
        "Turbine": [turbine_overhead_figure(turbine), turbine_blade_figure(turbine)],
    }

    for key, label, _ in PUMPS:
        pump = pumps[key]
        print(f"{label} pump: shaft power {pump.shaft_power/1e3:0.1f} kW, shaft speed {pump.shaft_speed*30/np.pi:0.0f} rpm")

    export_volute_step(pumps)
    show_dashboard(tabs, title="Mark 2 Engine")


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


def hide_geometry_grids(fig):
    """Gridlines are clutter on drawings: drop them on aspect-locked (geometry) axes, keeping axis lines and
    ticks for scale. Data subplots in the same figure keep their grid."""
    for y in fig.select_yaxes():
        if y.scaleanchor:
            y.showgrid = False
            fig.layout["xaxis" + (y.anchor or "x")[1:]].showgrid = False


def is_geometry(fig):
    """True if any axis in the figure is aspect-locked (see equal_aspect)."""
    return any(axis.scaleanchor for axis in fig.select_yaxes())


def rotate(x, y, angle):
    """Rotate points (x, y) about the origin by angle [rad]."""
    return x*np.cos(angle) - y*np.sin(angle), x*np.sin(angle) + y*np.cos(angle)


def line(fig, x, y, color, name, dash="solid", width=2, **subplot):
    """Line trace. Every trace is named; traces sharing a name share one legend entry (see dedupe_legend)."""
    fig.add_trace(go.Scatter(x=x, y=y, name=name, line=dict(color=color, dash=dash, width=width)), **subplot)


def shaded(fig, x, y, color, name, opacity=0.35, **subplot):
    """Closed outline filled translucent - for solid regions like the volute in plan view."""
    fig.add_trace(go.Scatter(x=x, y=y, name=name, fill="toself", opacity=opacity, mode="lines",
                             line=dict(color=color, width=1)), **subplot)


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

    rows = 2
    fig = make_subplots(rows=rows, cols=1, shared_xaxes=True, row_heights=[0.3, 0.7],
                        specs=[[{}], [{"secondary_y": True}]], vertical_spacing=0.06)

    # Chamber wall: inner and outer contour, mirrored about the axis
    line(fig, tca.x, tca.r1, "brown", name="Thrust Chamber Contours", row=1, col=1)
    for y in (tca.r2, -tca.r1, -tca.r2):
        line(fig, tca.x, y, "brown", name="Thrust Chamber Contours", row=1, col=1)
    # No equal aspect here: an aspect lock on one row of a shared-x stack makes every row re-zoom together
    fig.update_yaxes(title_text="Radius (m)", showgrid=False, row=1, col=1)  # not aspect-locked, but still a drawing
    fig.update_xaxes(showgrid=False, row=1, col=1)

    # Temperatures through the wall stack (left axis)
    temperatures = [(jacket.T_cold, "Cold Wall", "#1f77b4"), (jacket.T_hot, "Hot Wall", "#ff7f0e"),
                    (cool.T, "Regen Coolant", "#2ca02c"), (film.liquid.T, "Film Liquid", "#d62728"),
                    (film.gas.T, "Film Gas", "#9467bd"), (gas.T_free, "Free-stream Gas", "#8c564b"),
                    (gas.T_adiabatic_wall, "Adiabatic Wall Gas", "#e377c2"),
                    (gas.T_boundary, "Boundary Layer Reference", "#7f7f7f")]
    for y, name, color in temperatures:
        line(fig, tca.x, y, color, name=name, row=2, col=1)
    fig.update_yaxes(title_text="Temperature (K)", range=[0, 1.2 * tca.Tt], row=2, col=1, secondary_y=False)

    # Regen coolant pressure and velocity on their own axes, all stacked on the left; dotted to set them apart
    line(fig, tca.x, cool.p / 1e5, "blue", name="Regen Pressure", dash="dot", row=2, col=1, secondary_y=True)
    fig.update_yaxes(title_text="Regen Pressure (bar)", side="left", anchor="free", position=0.1, showgrid=False,
                     row=2, col=1, secondary_y=True)
    row_x, row_y = fig.data[-2].xaxis, fig.data[-2].yaxis  # the temperature (primary) axes of this row, e.g. 'x2', 'y2'
    fig.add_trace(go.Scatter(x=tca.x, y=cool.velocity, name="Regen Velocity", xaxis=row_x, yaxis="y4",
                             line=dict(color="red", dash="dot", width=2)))
    fig.update_layout(yaxis4=dict(title_text="Regen Velocity (m/s)", overlaying=row_y,
                                  side="left", anchor="free", position=0.0, showgrid=False))
    fig.update_xaxes(domain=[0.2, 1.0])  # room for the extra left axes
    fig.update_xaxes(title_text="Distance from Injector (m)", row=2, col=1)

    fig.update_layout(title="Thrust Chamber Geometry and Regenerative Coolant Properties")
    return fig, rows


#######################################################################################################
### Summary tables (HTML)
#######################################################################################################

def table(title, header, rows):
    """Small HTML table: header is the column labels, rows are [label, value, ...] lists."""
    head = "".join(f"<th>{h}</th>" for h in header)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in row) + "</tr>" for row in rows)
    return f"<div class='card'><h3>{title}</h3><table><tr>{head}</tr>{body}</table></div>"


def summary_html(pumps, turbine, link):
    ox, fuel = pumps["ox"], pumps["fuel"]
    rpm = lambda w: f"{w*30/np.pi:,.0f}"
    jacket = tca.regenerative_jacket

    # Mechanically linked parts share a shaft; the two shafts are tied only electrically (synchronous link)
    shafts = table("Shafts & Drive", ["Shaft", "Mechanically linked", "Speed", "Pole pairs"], [
        ["1", "Turbine + fuel pump + alternator", f"{rpm(fuel.shaft_speed)} rpm", link.alternator.pole_count // 2],
        ["2", "Ox pump + motor", f"{rpm(ox.shaft_speed)} rpm", link.motor.pole_count // 2],
        ["Electric link", "Alternator → cable → motor (synchronous)", f"{link.cable.frequency:,.0f} Hz",
         f"ratio {ox.shaft_speed/fuel.shaft_speed:0.3g} : 1"],
    ])

    def pump_rows(p):
        imp, ind = p.impeller[0], p.inducer
        discharge = "volute + cone" if "volute" in vars(p) else "vaned diffuser"
        return [f"{p.mdot:0.2f}", f"{p.vdot*1000:0.2f}", f"{p.p_in/1e5:0.2f}", f"{p.p_out/1e5:0.1f}",
                f"{p.head_rise:0.0f}", f"{imp.head_rise:0.0f}", f"{ind.head:0.0f}", f"{imp.specific_speed:0.3f}",
                f"{p.hydraulic_efficiency:0.3f}", f"{p.shaft_power/1e3:0.1f}", f"{imp.r_outlet*1000:0.1f}",
                imp.blade.count, f"{ind.r_tip[0]*1000:0.1f}", discharge]
    labels = ["Mass flow (kg/s)", "Volume flow (L/s)", "Inlet pressure (bar)", "Discharge pressure (bar)",
              "Pump head (m)", "Impeller head (m)", "Inducer head (m)", "Impeller ωs", "Hydraulic efficiency",
              "Shaft power (kW)", "Impeller r2 (mm)", "Impeller blades", "Inducer tip r (mm)", "Discharge"]
    pump_table = table("Pumps", ["", "Ox", "Fuel"], [[l, a, b] for l, a, b in zip(labels, pump_rows(ox), pump_rows(fuel))])

    turbine_table = table("Turbine & Gas Generator", ["", "Value"], [
        ["Shaft power (kW)", f"{turbine.shaft_power/1e3:0.1f}"],
        ["Efficiency", f"{turbine.power_efficiency:0.3f}"],
        ["Pitchline radius (mm)", f"{turbine.r_pitchline*1000:0.1f}"],
        ["Velocity ratio u/c0", f"{turbine.v_pitchline/turbine.stator.v_spouting_norm:0.3f}"],
        ["Rotor blades", turbine.rotor.blade_count],
        ["Nozzles", f"{turbine.stator.nozzle_number}, exit {turbine.stator.d_exit_nozzle*1000:0.1f} mm"],
        ["Admission arc", f"{np.degrees(turbine.stator.nozzle_arc/turbine.r_pitchline):0.0f} deg, area fraction {turbine.stator.admission_fraction:0.3f}"],
        ["GG flow (kg/s)", f"{gg.mdot:0.3f} ({gg.fraction*100:0.1f}% of engine)"],
        ["GG chamber", f"{gg.pc/1e5:0.1f} bar, {gg.Tt:0.0f} K, MR {gg.mixture_ratio:g}"],
    ])

    tca_table = table("Thrust Chamber", ["", "Value"], [
        ["Thrust (kN)", f"{engine.thrust/1e3:0.2f}"],
        ["Chamber pressure (bar)", f"{tca.pc/1e5:0.1f}"],
        ["Mass flow (kg/s)", f"{tca.mdot:0.2f}, MR {tca.mixture_ratio:g}"],
        ["c* (m/s)", f"{tca.c_star:0.0f}"],
        ["Barrel / throat / exit D (mm)", f"{2*tca.r1[0]*1000:0.1f} / {tca.d_throat*1000:0.1f} / {tca.d_exit*1000:0.1f}"],
        ["Expansion ratio", f"{tca.Ae_At:0.2f}"],
        ["Regen channels (barrel)", f"{int(jacket.n_pipe[0])} x {np.atleast_1d(jacket.h_pipe)[0]*1000:0.2f} mm deep"],
        ["Engine mass flow (kg/s)", f"{engine.mdot:0.2f}"],
    ])
    # Shafts and pumps stacked in the first column, turbine and TCA beside them
    return (f"<div class='cards'><div class='stack'>{shafts}{pump_table}</div>"
            f"<div class='stack'>{turbine_table}{tca_table}</div></div>")


#######################################################################################################
### TCA plan view: barrel wall + regen channels, throat and nozzle exit circles (drawn into the system view)
#######################################################################################################

def draw_tca_overhead(fig):
    """Barrel wall, regen channels, throat and nozzle exit circles in plan view."""
    jacket = tca.regenerative_jacket
    theta = np.linspace(0.0, 2.0 * np.pi, 400)

    # Barrel wall: hot-side liner (r1 to r2) plus the fin layer the channels are milled into (r2 to r2 + h).
    # The closeout over the channels isn't modelled, so it isn't drawn.
    r_inner, r_wall = tca.r1[0], tca.r2[0]
    n = int(jacket.n_pipe[0])
    h_pipe, gap_pipe = np.atleast_1d(jacket.h_pipe)[0], np.atleast_1d(jacket.gap_pipe)[0]  # m - at the barrel station
    r_out = r_wall + h_pipe
    shaded(fig, np.concatenate([r_out * np.cos(theta), r_inner * np.cos(theta[::-1])]) * 1000,
           np.concatenate([r_out * np.sin(theta), r_inner * np.sin(theta[::-1])]) * 1000, "brown", name="Barrel Wall")
    width = 2 * np.pi / n - gap_pipe / r_wall  # rad - channel angular width at the wall, fin removed
    for k in range(n):
        a = np.linspace(0.0, width, 8) + 2 * np.pi * k / n
        x = np.concatenate([r_wall * np.cos(a), r_out * np.cos(a[::-1]), [r_wall * np.cos(a[0])]])
        y = np.concatenate([r_wall * np.sin(a), r_out * np.sin(a[::-1]), [r_wall * np.sin(a[0])]])
        shaded(fig, x * 1000, y * 1000, "lightsteelblue", name=f"Regen Channels ({n})", opacity=1.0)

    for r, name, dash in ((tca.r_throat, "Throat", "dash"), (tca.r1[-1], "Nozzle Exit", "solid")):
        line(fig, r * np.cos(theta) * 1000, r * np.sin(theta) * 1000, "black", name=name, dash=dash, width=1.5)



#######################################################################################################
### Meridional impeller contours + meanline parameters
#######################################################################################################

def meridional_figure(pumps):
    rows = 2
    fig = make_subplots(rows=rows, cols=1, specs=[[{}], [{"secondary_y": True}]],
                        subplot_titles=("Meridional Plane Impeller Contours (mm)", ""),
                        vertical_spacing=0.12)

    # Ox mirrored to negative r so both impellers share one plot; the flow passage between hub and shroud is shaded
    for (key, label, color), side in zip(PUMPS, (-1.0, 1.0)):
        imp = pumps[key].impeller[0]
        passage = np.vstack([imp.shroud_curve, imp.hub_curve[::-1]])
        shaded(fig, side * passage[:, 0] * 1000, passage[:, 1] * 1000, color, name=f"{label} Passage", opacity=0.15,
               row=1, col=1)
        for curve, part, dash in ((imp.shroud_curve, "Shroud", "solid"), (imp.hub_curve, "Hub", "solid"),
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

        # Volute in plan view: shaded band between the inner and outer section edges, plus the discharge cone
        if 'volute' in vars(pump):
            vol = pump.volute
            azimuth = pump.clocking * vol.theta
            r_in = vol.r_meanline - 0.5 * vol.d_hydraulic_section
            r_out = vol.r_meanline + 0.5 * vol.d_hydraulic_section
            x = np.concatenate([r_out * np.cos(azimuth), (r_in * np.cos(azimuth))[::-1]])
            y = np.concatenate([r_out * np.sin(azimuth), (r_in * np.sin(azimuth))[::-1]])
            shaded(fig, x * 1000, y * 1000, color, name=f"{label} Volute")

            # Nozzle walls are the phi = 0 and phi = pi generators of the cone (in the volute midplane)
            n_phi = vol.nozzle_contour.shape[1]
            wall_out, wall_in = vol.nozzle_contour[:, 0, :], vol.nozzle_contour[:, (n_phi - 1) // 2, :]
            x = np.concatenate([wall_out[:, 0], wall_in[::-1, 0]])
            y = np.concatenate([wall_out[:, 1], wall_in[::-1, 1]])
            shaded(fig, x * 1000, y * 1000, color, name=f"{label} Discharge Cone")
        else:
            # Vaned diffuser: leading and trailing edge circles
            dif = pump.vaned_diffuser
            circle = np.linspace(0.0, 2.0 * np.pi, 200)
            for r in (dif.r_inlet, dif.r_outlet):
                line(fig, r * 1000 * np.cos(circle), r * 1000 * np.sin(circle), color, name=f"{label} Diffuser",
                     dash="dot", width=1.5)

    n_pumps = len(fig.data)

    # Turbine rotor: pitchline (dashed), blade tip and base
    theta = np.linspace(0.0, 2.0 * np.pi, 1000)
    for r, name, dash in ((turbine.r_pitchline, "Turbine Pitchline", "dash"),
                          (turbine.rotor.r_tip, "Turbine Blade Tip/Base", "solid"),
                          (turbine.rotor.r_base, "Turbine Blade Tip/Base", "solid")):
        line(fig, r * 1000 * np.cos(theta), r * 1000 * np.sin(theta), "orange", name=name, dash=dash, width=1)

    n_powerhead = len(fig.data)

    draw_tca_overhead(fig)

    # Crowded, so the legend is split into columns under the plot: pumps, turbine, thrust chamber
    for k, trace in enumerate(fig.data):
        trace.legend = "legend" if k < n_pumps else ("legend2" if k < n_powerhead else "legend3")
    for k, (name, title) in enumerate((("legend", "Pumps"), ("legend2", "Turbine"), ("legend3", "Thrust Chamber"))):
        fig.update_layout({name: dict(title_text=f"<b>{title}</b>", x=0.34 * k, xanchor="left", y=-0.08, yanchor="top")})

    equal_aspect(fig)
    fig.update_layout(title="System - Top View (mm)", margin_b=180)
    return fig, 2


#######################################################################################################
### Turbine top view: rotor blades + partial-admission nozzle arc
#######################################################################################################

def turbine_overhead_figure(turbine):
    fig = go.Figure()
    theta = np.linspace(0.0, 2.0 * np.pi, 400)
    r_pitch, r_base, r_tip = turbine.r_pitchline, turbine.rotor.r_base, turbine.rotor.r_tip

    # Blade row as a shaded annulus with one radial tick per blade
    shaded(fig, np.concatenate([r_tip * np.cos(theta), r_base * np.cos(theta[::-1])]) * 1000,
           np.concatenate([r_tip * np.sin(theta), r_base * np.sin(theta[::-1])]) * 1000, "orange", name="Blade Row",
           opacity=0.2)
    a = 2 * np.pi * np.arange(turbine.rotor.blade_count) / turbine.rotor.blade_count
    x = np.column_stack([r_base * np.cos(a), r_tip * np.cos(a), np.full_like(a, np.nan)]).ravel()  # NaN breaks the line
    y = np.column_stack([r_base * np.sin(a), r_tip * np.sin(a), np.full_like(a, np.nan)]).ravel()
    fig.add_trace(go.Scatter(x=x * 1000, y=y * 1000, mode="lines", name=f"Rotor Blades ({turbine.rotor.blade_count})",
                             line=dict(color="darkorange", width=1)))
    line(fig, r_pitch * np.cos(theta) * 1000, r_pitch * np.sin(theta) * 1000, "orange", name="Pitchline", dash="dash", width=1)

    # Active arc, starting at 0 deg: the jet sheet, then each nozzle exit as an ellipse on the pitchline.
    # Drawn symmetric; the true cone cut is lopsided (see stators.nozzle_footprint), but the arc length matches.
    arc = turbine.stator.nozzle_arc / r_pitch  # rad
    a = np.linspace(0.0, arc, 100)
    shaded(fig, np.concatenate([r_tip * np.cos(a), r_base * np.cos(a[::-1])]) * 1000,
           np.concatenate([r_tip * np.sin(a), r_base * np.sin(a[::-1])]) * 1000, "red", name="Active Arc", opacity=0.15)
    t = np.linspace(0.0, 2.0 * np.pi, 40)
    for k in range(turbine.stator.nozzle_number):
        center = (turbine.stator.nozzle_footprint / 2 + k * turbine.stator.nozzle_spacing) / r_pitch  # rad
        along, radial = turbine.stator.nozzle_footprint / 2 * np.cos(t), turbine.stator.d_exit_nozzle / 2 * np.sin(t)  # m - local ellipse
        x = (r_pitch + radial) * np.cos(center) - along * np.sin(center)
        y = (r_pitch + radial) * np.sin(center) + along * np.cos(center)
        line(fig, x * 1000, y * 1000, "gray", name=f"Nozzle Exits ({turbine.stator.nozzle_number})", width=1)

    info_box(fig, [f"Nozzles = {turbine.stator.nozzle_number}, exit D {turbine.stator.d_exit_nozzle*1000:0.1f} mm",
                   f"Footprint {turbine.stator.nozzle_footprint*1000:0.1f} mm, spacing {turbine.stator.nozzle_spacing*1000:0.1f} mm",
                   f"Arc = {np.degrees(arc):0.0f} deg ({turbine.stator.nozzle_arc*1000:0.0f} mm)",
                   f"Area fraction = {turbine.stator.admission_fraction:0.3f}"])
    equal_aspect(fig)
    fig.update_layout(title="Turbine - Top View with Nozzle Arc (mm)")
    return fig, 1


#######################################################################################################
### Turbine blade profiles (2D cascade)
#######################################################################################################

def turbine_blade_figure(turbine):
    fig = go.Figure()
    chord, pitch, r_le = turbine.rotor.blade_chord, turbine.rotor.blade_pitch, turbine.rotor.radius_leading

    # Trailing (pressure-side) arc through the blade tips
    sagitta = chord / 2.0 * (1.0 - np.cos(turbine.rotor.blade_angle)) / np.sin(turbine.rotor.blade_angle)
    tip = np.array([-sagitta, chord / 2.0])
    blade_radius = sagitta / 2.0 + chord**2 / (8.0 * sagitta)
    center_arc = np.array([-blade_radius, 0.0])
    angle_trailing = np.arctan2(tip[1], tip[0] + blade_radius)
    arc_theta = np.linspace(-angle_trailing, angle_trailing, 120)
    x_trailing = center_arc[0] + blade_radius * np.cos(arc_theta)
    y_trailing = center_arc[1] + blade_radius * np.sin(arc_theta)

    # Leading-edge fillet, joined to the tips by tangent lines
    apex_leading = np.array([pitch - turbine.rotor.blade_opening, 0.0])
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
    r_throat = turbine.stator.d_throat_nozzle / 2
    r_exit = r_throat * np.sqrt(gg.Ae_At)                                  # same expansion ratio as the GG nozzle
    length_diverging = (r_exit - r_throat) / np.tan(turbine.stator.diverge_angle)
    length_converging = 2 * r_throat                                       # DRAWING ONLY - short inlet stub, not modeled
    along = np.array([np.sin(turbine.stator.angle_nozzle), -np.cos(turbine.stator.angle_nozzle)])   # unit vector along the jet
    across = np.array([along[1], -along[0]])                                           # unit normal to the jet

    nozzle_spacing = turbine.stator.nozzle_spacing  # m - along the pitchline
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
            cone = rotate(along[0], along[1], -side * turbine.stator.diverge_angle)  # wall direction, opening away from the axis
            wall = np.array([inlet + side * 2 * r_throat * across, throat_wall, to_plate_face(throat_wall, np.array(cone))])
            line(fig, wall[:, 0] * 1000, wall[:, 1] * 1000, "gray", name="Nozzles")
        axis_line = np.array([inlet, axis_exit])
        line(fig, axis_line[:, 0] * 1000, axis_line[:, 1] * 1000, "gray", name="Nozzle Axes", dash="dashdot", width=1)
    fig.add_hline(y=plate_face_y * 1000, line=dict(color="gray", width=1, dash="dot"), name="Nozzle Plate Face",
                  showlegend=True)

    info_box(fig, [f"Blade Count = {turbine.rotor.blade_count:g}",
                             f"Nozzles = {turbine.stator.nozzle_number}, spacing {nozzle_spacing*1000:0.1f} mm",
                             f"Admission: arc {turbine.stator.nozzle_arc*1000:0.0f} mm, area fraction {turbine.stator.admission_fraction:0.3f}",
                             f"Nozzle tilt = {np.degrees(turbine.stator.angle_nozzle):0.1f} deg from axial",
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
            # Blade as a shaded surface (rows hub to tip), then leading and trailing edges drawn heavier
            fig.add_trace(go.Surface(x=x * 1000, y=y * 1000, z=z * 1000, name=f"{label} Blades", showlegend=k == 0,
                                     showscale=False, colorscale=[[0, color], [1, color]], opacity=0.85,
                                     lighting=dict(ambient=0.5, diffuse=0.8, specular=0.3)), row=1, col=col)
            for edge, edge_name in ((0, "Leading Edges"), (-1, "Trailing Edges")):
                fig.add_trace(go.Scatter3d(x=x[:, edge] * 1000, y=y[:, edge] * 1000, z=z[:, edge] * 1000,
                                           mode="lines", name=f"{label} {edge_name}",
                                           line=dict(color=color, width=4, dash="solid" if edge == 0 else "dash")),
                              row=1, col=col)

        # Hub as a surface of revolution along the blade's hub line
        theta = np.linspace(0.0, 2.0 * np.pi, 64)
        r_hub = np.hypot(ind.x_bladeline[0], ind.y_bladeline[0])[:, None]
        fig.add_trace(go.Surface(x=r_hub * np.cos(theta) * 1000, y=r_hub * np.sin(theta) * 1000,
                                 z=np.broadcast_to(ind.z_bladeline[0][:, None], r_hub.shape[:1] + theta.shape) * 1000,
                                 name=f"{label} Hub", showlegend=True, showscale=False,
                                 colorscale=[[0, "gray"], [1, "gray"]], opacity=0.6), row=1, col=col)

        info_box(fig, heading=f"{label} Inducer", lines=[
            f"phi = {ind.flow_coeff:0.4g}",
            f"psi = {ind.head_coeff:0.3g}",
            f"Nss = {2733.0*ind.suction_specific_speed_realized:0.0f}",  # Omega_ss -> US Nss (rpm, gpm, ft)
            f"Shaft = {30/np.pi*pump.shaft_speed:0.0f} RPM",
        ])

    camera = dict(eye=dict(x=1.8, y=1.8, z=1.2))
    fig.update_layout(scene=dict(aspectmode="data", camera=camera), scene2=dict(aspectmode="data", camera=camera),
                      title="Inducer 3D Geometry")
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
            noz = vol.nozzle_contour
            surface(noz[:, :, 0], noz[:, :, 1], noz[:, :, 2], color, 0.6, col, f"{label} Discharge Cone")

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

def show_dashboard(tabs, title="Mark 2 Plots", open_browser=True):
    """Stitch figures into one HTML page with a folder-tab row (one tab per dict key). Each tab holds figure
    sections; a (name, [figures]) entry becomes an expandable subsection. Every plot stays fully interactive
    (3D rotate, hover, zoom), and plotly.js is embedded once so the page works offline.
    tabs: {tab name: [(fig, n_rows) | (subsection name, [(fig, n_rows), ...]), ...]}
    """
    first = [True]  # embed plotly.js with the first figure only

    def figure_html(fig, rows):
        split_legends(fig)
        dedupe_legend(fig)
        heading = fig.layout.title.text or "Figure"  # becomes the collapsible section's header
        margin_b = max(10, fig.layout.margin.b or 0)  # keep room a figure reserved (e.g. legends under the plot)
        fig.update_layout(title=None, template="plotly_white", autosize=True, margin=dict(l=10, r=10, t=30, b=margin_b))
        fig.update_xaxes(showline=True, linecolor="#bbb", gridcolor="#e6e6e6", zerolinecolor="#cfcfcf")
        fig.update_yaxes(showline=True, linecolor="#bbb", gridcolor="#e6e6e6", zerolinecolor="#cfcfcf")
        hide_geometry_grids(fig)
        fig.update_scenes(xaxis_showgrid=False, yaxis_showgrid=False, zaxis_showgrid=False,
                          xaxis_showbackground=False, yaxis_showbackground=False, zaxis_showbackground=False)  # 3D plots are all geometry
        draw_info_panel(fig)
        plot = fig.to_html(full_html=False, include_plotlyjs=first[0], default_width="96%",
                           default_height=ROW_HEIGHT * rows,
                           config={"responsive": True, "scrollZoom": is_geometry(fig)})  # wheel zoom only on equal-aspect plots
        first[0] = False
        return f"<details open class='figure'><summary>{heading}</summary>{plot}</details>"

    buttons, panels = [], []
    for i, (name, entries) in enumerate(tabs.items()):
        parts = []
        for entry in entries:
            if isinstance(entry, str):  # raw HTML section (summary tables)
                parts.append(f"<details open><summary>Summary</summary>{entry}</details>")
            elif isinstance(entry[0], str):  # subsection
                inner = "\n".join(figure_html(*f) for f in entry[1])
                parts.append(f"<details open class='group'><summary>{entry[0]}</summary>{inner}</details>")
            else:
                parts.append(figure_html(*entry))
        active = " active" if i == 0 else ""
        buttons.append(f"<button class='tab{active}' onclick=\"showTab({i})\">{name}</button>")
        panels.append(f"<div class='panel{active}' id='tab{i}'>" + "\n".join(parts) + "</div>")

    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>{title}</title>
<style>
  body {{ margin:0; background:#f3f3f3;
         font-family:-apple-system,Segoe UI,Roboto,sans-serif; }}
  header {{ position:sticky; top:0; z-index:10; background:#1c1c1c; }}
  h1 {{ padding:14px 20px 8px; margin:0; font-size:18px; color:#eee; }}
  .tabs {{ display:flex; gap:4px; padding:0 20px; }}
  .tab {{ border:none; border-radius:8px 8px 0 0; padding:8px 22px; font-size:15px; cursor:pointer;
          background:#3a3a3a; color:#ccc; }}
  .tab:hover {{ color:#fff; }}
  .tab.active {{ background:#f3f3f3; color:#111; font-weight:600; }}
  .panel {{ display:none; }}
  .panel.active {{ display:block; }}
  details {{ margin:14px auto; padding:8px; max-width:1200px;
             background:#fff; border-radius:8px; box-shadow:0 1px 4px rgba(0,0,0,.12); }}
  details.group {{ background:#e9ecef; }}
  details.group > summary {{ font-weight:600; }}
  details details {{ margin:10px 0; }}
  summary {{ cursor:pointer; padding:6px 10px; font-size:16px; color:#333; }}
  summary:hover {{ color:#000; }}
  .cards {{ display:flex; flex-wrap:wrap; gap:14px; padding:6px 10px; align-items:flex-start; }}
  .stack {{ display:flex; flex-direction:column; gap:14px; }}
  .card h3 {{ margin:4px 0 6px; font-size:14px; color:#333; }}
  .card table {{ border-collapse:collapse; font-size:13px; }}
  .card th, .card td {{ padding:3px 10px; border-bottom:1px solid #e6e6e6; text-align:left; }}
  .card th {{ background:#f3f3f3; font-weight:600; }}
  .card td:not(:first-child) {{ font-variant-numeric:tabular-nums; }}
  .controls {{ max-width:1200px; margin:10px auto 0; text-align:right; }}
  .controls button {{ margin-left:6px; }}
</style></head><body>
<header><h1>{title}</h1><div class="tabs">{"".join(buttons)}</div></header>
<div class="controls"><button onclick="setAll(true)">Expand all</button><button onclick="setAll(false)">Collapse all</button></div>
{"".join(panels)}
<script>
  // Plots in a hidden tab or collapsed section have no size; re-fit them when they become visible.
  function refit(el) {{ el.querySelectorAll(".plotly-graph-div").forEach(p => Plotly.Plots.resize(p)); }}
  document.querySelectorAll("details").forEach(d => d.addEventListener("toggle", () => {{ if (d.open) refit(d); }}));
  function showTab(i) {{
    document.querySelectorAll(".tab").forEach((b, k) => b.classList.toggle("active", k === i));
    document.querySelectorAll(".panel").forEach((p, k) => p.classList.toggle("active", k === i));
    refit(document.getElementById("tab" + i));
    window.scrollTo(0, 0);
  }}
  function setAll(open) {{ document.querySelectorAll(".panel.active details").forEach(d => d.open = open); }}
</script>
</body></html>"""

    out = Path(__file__).resolve().parent / "mark2_plots.html"  # next to this script; ignored by git
    out.write_text(html, encoding="utf-8")
    if open_browser:
        webbrowser.open(out.as_uri())
    return out


if __name__ == "__main__":
    plots()
