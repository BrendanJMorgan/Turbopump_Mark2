import numpy as np
import matplotlib.pyplot as plt

from engine_state import engine, tca, gg, pump, turbine
from powerhead import powerhead


def plots():

    ox_pump = powerhead.ox_pump
    fuel_pump = powerhead.fuel_pump
    turbine = powerhead.turbine

    # -------------------------------------------------------------------------
    # TCA (thrust chamber assembly) Contours
    # -------------------------------------------------------------------------
    plt.figure(1)
    plt.clf()

    plt.plot(tca.x, tca.r1, tca.x, tca.r2, tca.x, -1 * tca.r1, tca.x, -1 * tca.r2, color="blue")
    plt.axis("equal")
    plt.xlabel("Distance from Injector (m)")
    plt.title("Combustion Chamber Contours")

    # -------------------------------------------------------------------------
    # TCA Temperatures
    # -------------------------------------------------------------------------

    # plt.figure(2)
    # plt.clf()
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
    # plt.xlabel("Distance from Injector (m)")
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
        -ox_pump.meanline_curve[:, 0] * 1000.0,
        ox_pump.meanline_curve[:, 1] * 1000.0,
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
        fuel_pump.meanline_curve[:, 0] * 1000.0,
        fuel_pump.meanline_curve[:, 1] * 1000.0,
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
    # delta_angle = clocking(1)*2*pi/blade_count(1);
    delta_angle = ox_pump.clocking * 2.0 * np.pi / ox_pump.blade_count  # Calculate the angle to rotate each blade

    for i in range(ox_pump.blade_count):  # for i = 0:(blade_count(1)-1)
        rotation_matrix = np.array(
            [
                [np.cos(i * delta_angle), -np.sin(i * delta_angle)],
                [np.sin(i * delta_angle), np.cos(i * delta_angle)],
            ]
        )
        # blade_curve(:,1:2,1) * rotation_matrix';
        rotated_curve = ox_pump.blade_curve[:, 0:2, 0] @ rotation_matrix.T

        plt.plot(
            rotated_curve[:, 0] * 1000.0 - 150.0,
            rotated_curve[:, 1] * 1000.0,
            linewidth=2,
            color="cyan",
        )
        # plot(NaN, NaN); % Prevent connection between different blades
        plt.plot(np.nan, np.nan)

    plt.plot(
        ox_pump.volute_curve[:, 0, 0] * 1000.0 - 150.0,
        ox_pump.volute_curve[:, 1, 0] * 1000.0,
        color="cyan",
    )

    # -------------------------------------------------------------------------
    # Fuel Pump Contours
    # -------------------------------------------------------------------------
    delta_angle = fuel_pump.clocking * 2.0 * np.pi / fuel_pump.blade_count  # Calculate the angle to rotate each blade

    for i in range(fuel_pump.blade_count[1]):  # for i = 0:(blade_count(2)-1)
        rotation_matrix = np.array(
            [
                [np.cos(i * delta_angle), -np.sin(i * delta_angle)],
                [np.sin(i * delta_angle), np.cos(i * delta_angle)],
            ]
        )
        rotated_curve = fuel_pump.blade_curve[:, 0:2, 1] @ rotation_matrix.T

        plt.plot(
            rotated_curve[:, 0] * 1000.0,
            rotated_curve[:, 1] * 1000.0,
            linewidth=2,
            color="red",
        )
        # Prevent connection between different blades
        plt.plot(np.nan, np.nan)

    plt.plot(
        fuel_pump.volute_curve[:, 0, 1] * 1000.0,
        fuel_pump.volute_curve[:, 1, 1] * 1000.0,
        color="red",
    )

    # -------------------------------------------------------------------------
    # Turbine Contours
    # -------------------------------------------------------------------------
    plt.figure(5)
    plt.clf()
    ax = plt.gca()
    ax.set_aspect("equal", adjustable="box")

    theta = np.linspace(0.0, 2.0 * np.pi, 1000)  # rad

    # plot(r_pitchline_rotor*1000*cos(theta), r_pitchline_rotor*1000*sin(theta), 'LineStyle','--','color','#ffA500')
    ax.plot(
        turbine.r_pitchline_rotor * 1000.0 * np.cos(theta),
        turbine.r_pitchline_rotor * 1000.0 * np.sin(theta),
        linestyle="--",
        color="#ffA500",
    )
    # plot(r_tip_rotor*1000*cos(theta), r_tip_rotor*1000*sin(theta),
    #      r_base_rotor*1000*cos(theta), r_base_rotor*1000*sin(theta), 'color','#ffA500')
    ax.plot(
        turbine.r_tip_rotor * 1000.0 * np.cos(theta),
        turbine.r_tip_rotor * 1000.0 * np.sin(theta),
        color="#ffA500",
    )
    ax.plot(
        turbine.r_base_rotor * 1000.0 * np.cos(theta),
        turbine.r_base_rotor * 1000.0 * np.sin(theta),
        color="#ffA500",
    )

    # str = sprintf('GG = %2.g kg/s, %0f SCFM (%.2g%%)', mdot_gg, mdot_gg*1790, mdot_gg/mdot_total*100);
    text_str = f"GG = {gg.mdot:0.2g} kg/s, {gg.mdot*1790:0.0f} SCFM ({gg.mdot/engine.mdot*100:.2g}%)"

    # annotation('textbox',[.2 .5 .3 .3],'String',str,'FitBoxToText','on');
    ax.text(
        0.2,
        0.5,
        text_str,
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
    )

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
    sagitta = turbine.blade_chord_rotor / 2.0 * (1.0 - np.cos(turbine.blade_angle_rotor)) / np.sin(
        turbine.blade_angle_rotor
    )  # % m
    tip = np.array([-sagitta, turbine.blade_chord_rotor / 2.0])  # % [m,m] - chord point
    blade_radius = sagitta / 2.0 + turbine.blade_chord_rotor**2 / (8.0 * sagitta)  # % m
    center_arc = np.array([-blade_radius, 0.0])  # % y=0
    angle_trailing = np.arctan2(tip[1] - 0.0, tip[0] + blade_radius)  # % rad
    arc_theta = np.linspace(-angle_trailing, angle_trailing, 120)  # %
    x_trailing = center_arc[0] + blade_radius * np.cos(arc_theta)  # % [m,m]
    y_trailing = center_arc[1] + blade_radius * np.sin(arc_theta)  # % [m,m]

    # Leading edge
    apex_leading = np.array(
        [turbine.blade_pitch_rotor - turbine.blade_opening_rotor, 0.0]
    )  # % [m,m] - location of blade apex point
    center_leading = apex_leading + np.array([-turbine.radius_leading, 0.0])  # % [m,m] on y=0
    center_tip = tip - center_leading  # % m
    base_angle = np.arctan2(center_tip[1], center_tip[0])  # % rad
    tan_angle = base_angle - np.arccos(turbine.radius_leading / np.linalg.norm(center_tip))  # % rad - upper tangent
    tangency_leading = center_leading + turbine.radius_leading * np.array(
        [np.cos(tan_angle), np.sin(tan_angle)]
    )  # % [m,m]
    theta_leading = np.linspace(-tan_angle, tan_angle, 60)  # % rad
    x_leading = center_leading[0] + turbine.radius_leading * np.cos(theta_leading)  # % [m,m]
    y_leading = center_leading[1] + turbine.radius_leading * np.sin(theta_leading)  # % [m,m]

    for off in [-turbine.blade_pitch_rotor * 1000.0, 0.0, turbine.blade_pitch_rotor * 1000.0]:
        ax.plot(
            x_trailing * 1000.0 + off,
            y_trailing * 1000.0,
            "y",
            linewidth=2,
        )  # % trailing edge arc
        ax.plot(
            x_leading * 1000.0 + off,
            y_leading * 1000.0,
            "y",
            linewidth=2,
        )  # % leading edge arc
        ax.plot(
            [tangency_leading[0] * 1000.0 + off, tip[0] * 1000.0 + off],
            [tangency_leading[1] * 1000.0, tip[1] * 1000.0],
            "y",
            linewidth=2,
        )  # % leading edge tangent line
        ax.plot(
            [tangency_leading[0] * 1000.0 + off, tip[0] * 1000.0 + off],
            [-tangency_leading[1] * 1000.0, -tip[1] * 1000.0],
            "y",
            linewidth=2,
        )  # % leading edge tangent line

    str_blades = f"Blade Count = {turbine.blade_count_rotor:g}"
    ax.text(
        0.2,
        0.5,
        str_blades,
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
    )