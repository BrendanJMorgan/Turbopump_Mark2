import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
plt.ion()

from engine_state import engine, tca, gg, pump, turbine
from powerhead import powerhead


def plots():
    plt.close('all')

    ox_pump = powerhead.ox_pump
    fuel_pump = powerhead.fuel_pump
    turbine = powerhead.turbine

    # -------------------------------------------------------------------------
    # TCA (thrust chamber assembly) Contours
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 8))

    # Top subplot - Chamber Contours
    p1, p2, = ax1.plot(tca.x, tca.r1, tca.x, tca.r2, label="Thrust Chamber Contours",color="brown")
    p3, p4 = ax1.plot(tca.x, -tca.r1, tca.x, -tca.r2, color="brown")
    ax1.set_aspect('equal', adjustable='datalim')
    ax1.margins(y=0.2)
    ax1.set_ylabel("Radius (m)")
    ax1.yaxis.label.set_color("brown")
    ax1.legend(handles=[p1], labelcolor="linecolor")

    # Bottom subplot - Coolant Properties (twin axis)
    ax2_twin = ax2.twinx()

    p1, = ax2.plot(tca.x, tca.p_cool / 1e5, label="Regen Pressure", color="blue")
    p2, = ax2_twin.plot(tca.x, tca.v_cool, label="Regen Velocity", color="red")
    ax2.set_xlabel("Distance from Injector (m)")
    ax2.set_ylabel("Pressure (bar)")
    ax2_twin.set_ylabel("Velocity (m/s)")
    ax2.yaxis.label.set_color(p1.get_color())
    ax2_twin.yaxis.label.set_color(p2.get_color())
    ax2.legend(handles=[p1, p2], labelcolor="linecolor")

    plt.tight_layout()
    plt.show()

    # -------------------------------------------------------------------------
    # TCA Temperatures
    # -------------------------------------------------------------------------

    # plt.plot(
    #     tca.x, tca.T_wall_cold,
    #     tca.x, tca.T_wall_hot,
    #     tca.x, tca.T_cool,
    #     tca.x, tca.T_film,
    #     tca.x, tca.T_free,
    #     tca.x, tca.T_ab,
    #     tca.x, tca.T_recovery,
    #     tca.x, tca.T_ref,
    # )
    # plt.axhline(0.0)
    # plt.legend(
    #     [
    #         "Cold Wall",
    #         "Hot Wall",
    #         "Regen Coolant",
    #         "Film Coolant",
    #         "Free-Stream Gas",
    #         "Adiabatic",
    #         "Recovery",
    #         "Gas Property Reference",
    #     ],
    #     loc="northeast",
    # )
    plt.xlabel("Distance from Injector (m)")
    # plt.ylabel("Temperature (K)")
    # plt.title("Engine Steady-State Temperatures")

    # -------------------------------------------------------------------------
    # Pump Impeller and Blades
    # -------------------------------------------------------------------------
    plt.figure(3)
    plt.clf()
    plt.grid(True)
    plt.axis("equal")
    plt.title("Impeller and Shroud Contours (mm)")

    # --- Ox pump on the left (mirrored in x, cyan) ---
    plt.plot(
        -ox_pump.shroud_curve[:, 0] * 1000.0,
        ox_pump.shroud_curve[:, 1] * 1000.0,
        color="cyan",
    )
    plt.plot(
        -ox_pump.hub_curve[:, 0] * 1000.0,
        ox_pump.hub_curve[:, 1] * 1000.0,
        color="cyan",
    )
    plt.plot(
        -ox_pump.meanline_curve_bladed[:, 0] * 1000.0,
        ox_pump.meanline_curve_bladed[:, 1] * 1000.0,
        "--",
        color="cyan",
    )

    # --- Fuel pump on the right (red) ---
    plt.plot(
        fuel_pump.shroud_curve[:, 0] * 1000.0,
        fuel_pump.shroud_curve[:, 1] * 1000.0,
        color="red",
    )
    plt.plot(
        fuel_pump.hub_curve[:, 0] * 1000.0,
        fuel_pump.hub_curve[:, 1] * 1000.0,
        color="red",
    )
    plt.plot(
        fuel_pump.meanline_curve_bladed[:, 0] * 1000.0,
        fuel_pump.meanline_curve_bladed[:, 1] * 1000.0,
        "--",
        color="red",
    )


    ax = plt.gca()
    y_limits = ax.get_ylim()
    x_limits = ax.get_xlim()
    ax.plot([0.0, 0.0], y_limits, color="green")
    ax.plot(x_limits, [0.0, 0.0], color="green")

    

    # -------------------------------------------------------------------------
    # Impeller Blades / Volute
    # -------------------------------------------------------------------------
    plt.figure(4)
    plt.clf()
    plt.grid(True)
    plt.axis("equal")
    plt.title("Impellers, Volutes, and Rotor (mm)")

    # -------------------------------------------------------------------------
    # Ox Pump Contours
    # -------------------------------------------------------------------------
    delta_angle = ox_pump.clocking * 2.0 * np.pi / ox_pump.blade_count  # Calculate the angle to rotate each blade

    for i in range(ox_pump.blade_count):  # for i = 0:(blade_count(1)-1)
        c = ox_pump.blade_curve
        rotated_curve = np.vstack([c[:,0]*np.cos(i*delta_angle) - c[:,1]*np.sin(i*delta_angle), 
                                c[:,0]*np.sin(i*delta_angle) + c[:,1]*np.cos(i*delta_angle)]).T
        plt.plot(
            rotated_curve[:, 0] * 1000.0,
            rotated_curve[:, 1] * 1000.0,
            linewidth=2,
            color="cyan",
        )
        plt.plot(np.nan, np.nan)

    # plt.plot(
    #     ox_pump.volute_curve[:, 0, 0] * 1000.0 - 150.0,
    #     ox_pump.volute_curve[:, 1, 0] * 1000.0,
    #     color="cyan",
    # )

    # -------------------------------------------------------------------------
    # Fuel Pump Contours
    # -------------------------------------------------------------------------
    delta_angle = fuel_pump.clocking * 2.0 * np.pi / fuel_pump.blade_count  # Calculate the angle to rotate each blade

    for i in range(fuel_pump.blade_count):
        c = fuel_pump.blade_curve
        rotated_curve = np.vstack([c[:,0]*np.cos(i*delta_angle) - c[:,1]*np.sin(i*delta_angle), 
                                c[:,0]*np.sin(i*delta_angle) + c[:,1]*np.cos(i*delta_angle)]).T
        plt.plot(
            rotated_curve[:, 0] * 1000.0 - 150.0,
            rotated_curve[:, 1] * 1000.0,
            linewidth=2,
            color="red",
        )
        plt.plot(np.nan, np.nan)
    # plt.plot(
    #     fuel_pump.volute_curve[:, 0, 1] * 1000.0,
    #     fuel_pump.volute_curve[:, 1, 1] * 1000.0,
    #     color="red",
    # )

    # 

    # -------------------------------------------------------------------------
    # Turbine Contours
    # -------------------------------------------------------------------------
    ax = plt.gca()
    ax.set_aspect("equal", adjustable="box")

    theta = np.linspace(0.0, 2.0 * np.pi, 1000)  # rad

    ax.plot(
        turbine.r_pitchline * 1000.0 * np.cos(theta),
        turbine.r_pitchline * 1000.0 * np.sin(theta),
        linestyle="--",
        color="#ffA500",
    )
    ax.plot(
        turbine.r_tip * 1000.0 * np.cos(theta),
        turbine.r_tip * 1000.0 * np.sin(theta),
        color="#ffA500",
    )
    ax.plot(
        turbine.r_base * 1000.0 * np.cos(theta),
        turbine.r_base * 1000.0 * np.sin(theta),
        color="#ffA500",
    )

    text_str = f"GG = {gg.mdot:0.2g} kg/s, {gg.mdot*1790:0.0f} SCFM ({gg.mdot/engine.mdot*100:.2g}%)"

    ax.text(0.2,0.5,text_str,transform=ax.transAxes,bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    # -------------------------------------------------------------------------
    # Turbine Blades
    # -------------------------------------------------------------------------
    plt.figure(5)
    plt.clf()
    ax = plt.gca()
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True)
    ax.set_title("Rotor Blade Geometry (mm)")

    # Trailing edge
    sagitta = turbine.blade_chord / 2.0 * (1.0 - np.cos(turbine.blade_angle)) / np.sin(turbine.blade_angle)  # m
    tip = np.array([-sagitta, turbine.blade_chord / 2.0])  # [m,m] - chord point
    blade_radius = sagitta / 2.0 + turbine.blade_chord**2 / (8.0 * sagitta)  # m
    center_arc = np.array([-blade_radius, 0.0])  # y=0
    angle_trailing = np.arctan2(tip[1] - 0.0, tip[0] + blade_radius)  # rad
    arc_theta = np.linspace(-angle_trailing, angle_trailing, 120)  # rad
    x_trailing = center_arc[0] + blade_radius * np.cos(arc_theta)  # [m,m]
    y_trailing = center_arc[1] + blade_radius * np.sin(arc_theta)  # [m,m]

    # Leading edge
    apex_leading = np.array(
        [turbine.blade_pitch - turbine.blade_opening, 0.0]
    )  # [m,m] - location of blade apex point
    center_leading = apex_leading + np.array([-turbine.radius_leading, 0.0])  # [m,m] on y=0
    center_tip = tip - center_leading  # m
    base_angle = np.arctan2(center_tip[1], center_tip[0])  # rad
    tan_angle = base_angle - np.arccos(turbine.radius_leading / np.linalg.norm(center_tip))  # rad - upper tangent
    tangency_leading = center_leading + turbine.radius_leading * np.array(
        [np.cos(tan_angle), np.sin(tan_angle)]
    )  # [m,m]
    theta_leading = np.linspace(-tan_angle, tan_angle, 60)  # rad
    x_leading = center_leading[0] + turbine.radius_leading * np.cos(theta_leading)  # [m,m]
    y_leading = center_leading[1] + turbine.radius_leading * np.sin(theta_leading)  # [m,m]

    for off in [-turbine.blade_pitch * 1000.0, 0.0, turbine.blade_pitch * 1000.0]:
        ax.plot(
            x_trailing * 1000.0 + off,
            y_trailing * 1000.0,
            "y",
            linewidth=2,
        )  # trailing edge arc
        ax.plot(
            x_leading * 1000.0 + off,
            y_leading * 1000.0,
            "y",
            linewidth=2,
        )  # leading edge arc
        ax.plot(
            [tangency_leading[0] * 1000.0 + off, tip[0] * 1000.0 + off],
            [tangency_leading[1] * 1000.0, tip[1] * 1000.0],
            "y",
            linewidth=2,
        )  # leading edge tangent line
        ax.plot(
            [tangency_leading[0] * 1000.0 + off, tip[0] * 1000.0 + off],
            [-tangency_leading[1] * 1000.0, -tip[1] * 1000.0],
            "y",
            linewidth=2,
        )  # leading edge tangent line

    str_blades = f"Blade Count = {turbine.blade_count:g}"
    ax.text(
        0.2,
        0.5,
        str_blades,
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
    )

    plt.show(block=False)
    input("Press Enter to exit...")