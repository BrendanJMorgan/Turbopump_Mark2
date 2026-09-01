import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from engine_state import engine, tca, gg, pump, turbine
from powerhead import powerhead


def plots():
    ox_pump = powerhead.ox_pump
    fuel_pump = powerhead.fuel_pump
    turbine = powerhead.turbine

    gas = tca.combustion_gas
    jacket = tca.regenerative_jacket
    cool = tca.regenerative_coolant
    film = tca.film_coolant

    figures = []

    #######################################################################################################
    ### TCA Contours + Coolant Properties
    #######################################################################################################
    tca_fig_rows = 3
    tca_fig = make_subplots(rows=tca_fig_rows, cols=1, shared_xaxes=True,
        specs=[[{}], [{"secondary_y": True}], [{}]],
        vertical_spacing=0.08)


    tca_fig.add_trace(go.Scatter(x=tca.x, y=tca.r1, line=dict(color="brown"),
                            name="Thrust Chamber Contours"), row=1, col=1)
    for y in (tca.r2, -tca.r1, -tca.r2):
        tca_fig.add_trace(go.Scatter(x=tca.x, y=y, line=dict(color="brown"),
                                showlegend=False), row=1, col=1)

    tca_fig.add_trace(go.Scatter(x=tca.x, y=cool.p / 1e5, line=dict(color="blue"),
                            name="Regen Pressure"), row=2, col=1, secondary_y=False)
    tca_fig.add_trace(go.Scatter(x=tca.x, y=cool.velocity, line=dict(color="red"),
                            name="Regen Velocity"), row=2, col=1, secondary_y=True)

    tca_fig.update_yaxes(title_text="Radius (m)", scaleanchor="x", scaleratio=1,
                    row=1, col=1)
    tca_fig.update_yaxes(title_text="Pressure (bar)", color="blue",
                    row=2, col=1, secondary_y=False)
    tca_fig.update_yaxes(title_text="Velocity (m/s)", color="red",
                    row=2, col=1, secondary_y=True)

    temps = [(jacket.T_cold, "Cold Wall"), (jacket.T_hot, "Hot Wall"),
           (cool.T, "Regen Coolant"), (film.liquid.T, "Film Liquid"),
           (film.gas.T, "Film Gas"), (gas.T_free, "Free-blade Gas"),
           (gas.T_adiabatic_wall, "Adiabatic Wall Gas"),
           (gas.T_boundary, "Boundary Layer Reference")]
    for y, name in temps:
        tca_fig.add_trace(go.Scatter(x=tca.x, y=y, name=name), row=3, col=1)
    tca_fig.update_yaxes(title_text="Temperature (K)", range=[0, 1.2 * tca.Tt],
                    row=3, col=1)

    tca_fig.update_yaxes(title_text="Temperature (K)", row=3, col=1)
    tca_fig.update_xaxes(title_text="Distance from Injector (m)", row=3, col=1)
    tca_fig.update_layout(title="Thrust Chamber Geometry and Regenerative Coolant Properties")
    figures.append((tca_fig, 3))

    #######################################################################################################
    ### Meridional impeller contours + meanline parameters
    #######################################################################################################
    meridional_fig_rows = 2
    meridional_fig = make_subplots(rows=meridional_fig_rows, cols=1,
                       specs=[[{}], [{"secondary_y": True}]],
                       subplot_titles=("Meridional Plane Impeller Contours (mm)", ""),
                       vertical_spacing=0.12)

    def meridional(imp, sign, color, name):
        for curve, dash in [(imp.shroud_curve, "solid"),
                            (imp.hub_curve, "solid"),
                            (imp.blade.meanline_curve, "dash")]:
            meridional_fig.add_trace(go.Scatter(
                x=sign * curve[:, 0] * 1000.0, y=curve[:, 1] * 1000.0,
                line=dict(color=color, dash=dash), name=name,
                showlegend=False), row=1, col=1)

    meridional(ox_pump.impeller[0], -1.0, "blue", "Ox Pump")
    meridional(fuel_pump.impeller[0], 1.0, "red", "Fuel Pump")
    meridional_fig.update_yaxes(title_text="r (mm)", scaleanchor="x", scaleratio=1,
                    row=1, col=1)
    meridional_fig.update_xaxes(title_text="z (mm)", row=1, col=1)

    for imp, color, name in [(ox_pump.impeller[0], "blue", "Ox Pump"),
                             (fuel_pump.impeller[0], "red", "Fuel Pump")]:
        meridional_fig.add_trace(go.Scatter(x=imp.arc_meanline * 1000, y=imp.area_meanline * 1e6,
                                line=dict(color=color), name=name),
                     row=2, col=1, secondary_y=False)
        meridional_fig.add_trace(go.Scatter(x=imp.arc_meanline * 1000, y=imp.blockage,
                                line=dict(color=color, dash="dash"),
                                name=f"{name} Blockage", showlegend=False),
                     row=2, col=1, secondary_y=True)
    meridional_fig.update_xaxes(title_text="Meanline Arc Length (mm)", row=2, col=1)
    meridional_fig.update_yaxes(title_text="Annular Flow Area (mm2)", row=2, col=1, secondary_y=False)
    meridional_fig.update_yaxes(title_text="Blockage", row=2, col=1, secondary_y=True)
    meridional_fig.update_layout(title="Impeller Meridional & Meanline")
    figures.append((meridional_fig, meridional_fig_rows))

    #######################################################################################################
    ### Impellers, volutes, rotor (2D) + turbine contour circles
    #######################################################################################################
    overhead_fig_rows = 1
    overhead_fig = go.Figure()
    xs, ys = [], []

    def blades(pump, color, x_off):
        delta = pump.clocking * 2.0 * np.pi / pump.impeller[0].blade.count
        c = pump.impeller[0].blade_curve
        for i in range(pump.impeller[0].blade.count):
            rot = np.vstack([
                c[:, 0] * np.cos(i * delta) - c[:, 1] * np.sin(i * delta),
                c[:, 0] * np.sin(i * delta) + c[:, 1] * np.cos(i * delta)]).T
            X = rot[:, 0] * 1000.0 + x_off
            Y = rot[:, 1] * 1000.0
            xs.append(X); ys.append(Y)
            overhead_fig.add_trace(go.Scatter(x=X, y=Y, line=dict(color=color, width=2),
                                    showlegend=False))

    blades(ox_pump, "blue", 0.0)
    blades(fuel_pump, "red", -150.0)

    theta = np.linspace(0.0, 2.0 * np.pi, 1000)
    for r, dash in [(turbine.r_pitchline, "dash"),
                    (turbine.r_tip, "solid"),
                    (turbine.r_base, "solid")]:
        overhead_fig.add_trace(go.Scatter(x=r * 1000.0 * np.cos(theta),
                                y=r * 1000.0 * np.sin(theta),
                                line=dict(color="#ffA500", dash=dash),
                                showlegend=False))

    txt = (f"Shaft Speed = {turbine.shaft_speed*30/np.pi:0.2f} rpm<br>"
           f"GG = {gg.mdot:0.2g} kg/s, {gg.mdot*1790:0.0f} SCFM N2 "
           f"({gg.mdot/engine.mdot*100:.2g}%)")
    overhead_fig.add_annotation(x=0.2, y=0.1, xref="paper", yref="paper", showarrow=False,
                      text=txt, align="left", bgcolor="white", bordercolor="black",
                      borderpad=4, opacity=0.85)
    overhead_fig.update_yaxes(scaleanchor="x", scaleratio=1)
    overhead_fig.update_layout(title="Impellers, Volutes, and Rotor (mm)")
    figures.append((overhead_fig, overhead_fig_rows))

    #######################################################################################################
    ### Turbine Flowpath Geometry (2D)
    #######################################################################################################
    turbine_flow_fig_rows = 1
    turbine_flow_fig = go.Figure()

    sagitta = turbine.blade_chord / 2.0 * (1.0 - np.cos(turbine.blade_angle)) / np.sin(turbine.blade_angle)
    tip = np.array([-sagitta, turbine.blade_chord / 2.0])
    blade_radius = sagitta / 2.0 + turbine.blade_chord**2 / (8.0 * sagitta)
    center_arc = np.array([-blade_radius, 0.0])
    angle_trailing = np.arctan2(tip[1] - 0.0, tip[0] + blade_radius)
    arc_theta = np.linspace(-angle_trailing, angle_trailing, 120)
    x_trailing = center_arc[0] + blade_radius * np.cos(arc_theta)
    y_trailing = center_arc[1] + blade_radius * np.sin(arc_theta)

    apex_leading = np.array([turbine.blade_pitch - turbine.blade_opening, 0.0])
    center_leading = apex_leading + np.array([-turbine.radius_leading, 0.0])
    center_tip = tip - center_leading
    base_angle = np.arctan2(center_tip[1], center_tip[0])
    tan_angle = base_angle - np.arccos(turbine.radius_leading / np.linalg.norm(center_tip))
    tangency_leading = center_leading + turbine.radius_leading * np.array([np.cos(tan_angle), np.sin(tan_angle)])
    theta_leading = np.linspace(-tan_angle, tan_angle, 60)
    x_leading = center_leading[0] + turbine.radius_leading * np.cos(theta_leading)
    y_leading = center_leading[1] + turbine.radius_leading * np.sin(theta_leading)

    for off in [-turbine.blade_pitch * 1000.0, 0.0, turbine.blade_pitch * 1000.0]:
        for X, Y in [(x_trailing * 1000.0 + off, y_trailing * 1000.0),
                     (x_leading * 1000.0 + off, y_leading * 1000.0)]:
            turbine_flow_fig.add_trace(go.Scatter(x=X, y=Y, line=dict(color="gold", width=2),
                                    showlegend=False))
        # leading-edge tangent lines (both sides)
        for s in (1.0, -1.0):
            turbine_flow_fig.add_trace(go.Scatter(
                x=[tangency_leading[0] * 1000.0 + off, tip[0] * 1000.0 + off],
                y=[s * tangency_leading[1] * 1000.0, s * tip[1] * 1000.0],
                line=dict(color="gold", width=2), showlegend=False))

    turbine_flow_fig.add_annotation(x=0.2, y=0.5, xref="paper", yref="paper", showarrow=False,
                      text=f"Blade Count = {turbine.blade_count:g}",
                      bgcolor="white", bordercolor="black", borderpad=4, opacity=0.85)
    turbine_flow_fig.update_yaxes(scaleanchor="x", scaleratio=1)
    turbine_flow_fig.update_layout(title="Rotor Blade Geometry (mm)")
    figures.append((turbine_flow_fig, turbine_flow_fig_rows))

    #######################################################################################################
    ### Inducer Geometry
    #######################################################################################################
    inducer_fig_rows = 1
    inducer_fig = make_subplots(rows=1, cols=2,
                       specs=[[{"type": "scene"}, {"type": "scene"}]],
                       subplot_titles=("Ox Inducer Geometry (mm)",
                                       "Fuel Inducer Geometry (mm)"))

    def inducer(pump_obj, color, scene_col):
        ind = pump_obj.inducer
        n_sl = ind.x_bladeline.shape[0]
        for k in range(ind.blade_number):
            phi = 2 * np.pi * k / ind.blade_number
            xr = ind.x_bladeline * np.cos(phi) - ind.y_bladeline * np.sin(phi)
            yr = ind.x_bladeline * np.sin(phi) + ind.y_bladeline * np.cos(phi)
            zr = ind.z_bladeline
            for n in range(n_sl):
                inducer_fig.add_trace(go.Scatter3d(
                    x=xr[n] * 1000, y=yr[n] * 1000, z=zr[n] * 1000, mode="lines",
                    line=dict(color=color, width=2), showlegend=False),
                    row=1, col=scene_col)
            for edge in (0, -1):
                inducer_fig.add_trace(go.Scatter3d(
                    x=xr[:, edge] * 1000, y=yr[:, edge] * 1000, z=zr[:, edge] * 1000,
                    mode="lines", line=dict(color=color, width=4),
                    showlegend=False), row=1, col=scene_col)

    inducer(ox_pump, "blue", 1)
    inducer(fuel_pump, "red", 2)
    inducer_fig.update_layout(scene=dict(aspectmode="data"), scene2=dict(aspectmode="data"),
                     title="Inducer 3D Geometry")

    t1 = (f"phi = {ox_pump.inducer.flow_coeff:0.4g}<br>"
          f"psi = {ox_pump.inducer.head_coeff:0.3g}<br>"
          f"Nss = {2733.00*ox_pump.inducer.suction_specific_speed_available:0.0f}<br>"
          f"Shaft = {30/np.pi*ox_pump.shaft_speed:0.0f} RPM")
    t2 = (f"phi = {fuel_pump.inducer.flow_coeff:0.4g}<br>"
          f"psi = {fuel_pump.inducer.head_coeff:0.3g}<br>"
          f"Nss = {2733.00*fuel_pump.inducer.suction_specific_speed_available:0.0f}<br>"
          f"Shaft = {30/np.pi*fuel_pump.shaft_speed:0.0f} RPM")
    inducer_fig.add_annotation(x=0.0, y=0.0, xref="paper", yref="paper", text=t1,
                      showarrow=False, align="left", bgcolor="white", borderpad=4)
    inducer_fig.add_annotation(x=0.55, y=0.0, xref="paper", yref="paper", text=t2,
                      showarrow=False, align="left", bgcolor="white", borderpad=4)
    figures.append((inducer_fig, inducer_fig_rows))

    #######################################################################################################
    ### Volute Geometry
    #######################################################################################################
    volute_fig_rows = 1
    volute_fig = make_subplots(rows=1, cols=2,
                       specs=[[{"type": "scene"}, {"type": "scene"}]],
                       subplot_titles=("Ox Volute Geometry (mm)",
                                       "Fuel Volute Geometry (mm)"))

    def volute(pump_obj, color, scene_col):
        vol = pump_obj.volute
        imp = pump_obj.impeller[0]
        X = vol.contour[:, :, 0]; Y = vol.contour[:, :, 1]; Z = vol.contour[:, :, 2]
        volute_fig.add_trace(go.Surface(x=X * 1000, y=Y * 1000, z=Z * 1000, opacity=0.6,
                                showscale=False, colorscale=[[0, color], [1, color]]),
                     row=1, col=scene_col)
        theta_rev = np.linspace(0.0, 2.0 * np.pi, 64)
        for curve, col in [(imp.shroud_curve, "gray"), (imp.hub_curve, "lightgray")]:
            r = curve[:, 0][:, None]; z = curve[:, 1][:, None]
            Xr = r * np.cos(theta_rev)[None, :] * 1000
            Yr = r * np.sin(theta_rev)[None, :] * 1000
            Zr = np.broadcast_to(z, Xr.shape) * 1000
            volute_fig.add_trace(go.Surface(x=Xr, y=Yr, z=Zr, opacity=0.25, showscale=False,
                                    colorscale=[[0, col], [1, col]]),
                         row=1, col=scene_col)

    volute(ox_pump, "blue", 1)
    volute(fuel_pump, "red", 2)
    volute_fig.update_layout(scene=dict(aspectmode="data"), scene2=dict(aspectmode="data"),
                     title="Volute 3D Geometry")

    v1 = (f"method = {ox_pump.volute.design_method}<br>"
          f"d_throat = {ox_pump.volute.d_throat*1000:0.2f} mm<br>"
          f"c_throat = {ox_pump.volute.c_throat:0.1f} m/s<br>"
          f"loss = {ox_pump.volute.total_loss:0.4g}")
    v2 = (f"method = {fuel_pump.volute.design_method}<br>"
          f"d_throat = {fuel_pump.volute.d_throat*1000:0.2f} mm<br>"
          f"c_throat = {fuel_pump.volute.c_throat:0.1f} m/s<br>"
          f"loss = {fuel_pump.volute.total_loss:0.4g}")
    volute_fig.add_annotation(x=0.0, y=0.0, xref="paper", yref="paper", text=v1,
                      showarrow=False, align="left", bgcolor="white", borderpad=4)
    volute_fig.add_annotation(x=0.55, y=0.0, xref="paper", yref="paper", text=v2,
                      showarrow=False, align="left", bgcolor="white", borderpad=4)
    figures.append((volute_fig, volute_fig_rows))

    show_dashboard(figures, title="Mark 2 Engine")

import os
import tempfile
import webbrowser

def show_dashboard(figures, title="Mark 2 Plots", open_browser=True):
    """Stitch a list of Plotly figures into a single scrollable HTML page
    where every plot stays fully interactive (3D rotate, hover, zoom).

    plotly.js is embedded once (works offline). Each figure keeps its own
    title, which becomes its section heading.
    """
    fragments = []
    call_js = True
    for fig in figures:
        fig[0].update_layout(
            autosize=True,
            margin=dict(l=10, r=10, t=40, b=10),
            paper_bgcolor="white",
        )
        fragments.append(
            fig[0].to_html(
                full_html=False,
                include_plotlyjs=call_js,  # embed once    
                default_width="96%",
                default_height=400*fig[1],
                config={"responsive": True},
            )
        )
        call_js = False,

    body = "\n".join(f'<section>{f}</section>' for f in fragments)
    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>{title}</title>
<style>
  body {{ margin:0; background:#overhead_figoverhead_figturbine_flow_fig;
         font-family:-apple-system,Segoe UI,Roboto,sans-serif; }}
  h1 {{ padding:14px 20px; margin:0; font-size:18px;
        background:#1c1c1c; color:#eee; position:sticky; top:0; z-index:10; }}
  section {{ margin:14px auto; padding:8px; max-width:1200px;
             background:#fff; border-radius:8px;
             box-shadow:0 1px 4px rgba(0,0,0,.12); }}
</style></head><body>
<h1>{title}</h1>
{body}
</body></html>"""

    out = os.path.join(tempfile.gettempdir(), "mark2_plots.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    if open_browser:
        webbrowser.open("file://" + os.path.abspath(out))
    return out


if __name__ == "__main__":
    plots()
