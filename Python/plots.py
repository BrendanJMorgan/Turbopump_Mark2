import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
plt.ion()

from engine_state import engine, tca, gg, pump, turbine
from powerhead import powerhead


def plots():
    plt.close('all')

    ox_pump = powerhead.ox_pump
    fuel_pump = powerhead.fuel_pump
    turbine = powerhead.turbine

    gas = tca.combustion_gas
    jacket = tca.regenerative_jacket
    cool = tca.regenerative_coolant
    film = tca.film_coolant

    #######################################################################################################
    ### TCA (thrust chamber assembly) Contours
    #######################################################################################################

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

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

    p1, = ax2.plot(tca.x, tca.regenerative_coolant.p / 1e5, label="Regen Pressure", color="blue")
    p2, = ax2_twin.plot(tca.x, tca.regenerative_coolant.velocity, label="Regen Velocity", color="red")
    ax2.set_xlabel("Distance from Injector (m)")
    ax2.set_ylabel("Pressure (bar)")
    ax2_twin.set_ylabel("Velocity (m/s)")
    ax2.yaxis.label.set_color(p1.get_color())
    ax2_twin.yaxis.label.set_color(p2.get_color())
    ax2.legend(handles=[p1, p2], labelcolor="linecolor")

    plt.tight_layout()
    plt.show()

    #######################################################################################################
    ### TCA Temperatures
    #######################################################################################################

    fig, (ax3, ax4) = plt.subplots(2, 1, sharex=True)

    # Top plot - all 8 temperatures
    ax3.plot(tca.x, jacket.T_cold)
    ax3.plot(tca.x, jacket.T_hot)
    ax3.plot(tca.x, cool.T)
    ax3.plot(tca.x, film.liquid.T)
    ax3.plot(tca.x, film.gas.T)
    ax3.plot(tca.x, gas.T_free)
    ax3.plot(tca.x, gas.T_adiabatic_wall)
    ax3.plot(tca.x, gas.T_boundary)
    ax3.legend([
        "Cold Wall",
        "Hot Wall",
        "Regen Coolant",
        "Film Liquid",
        "Film Gas",
        "Free-blade Gas",
        "Adiabatic Wall Gas",
        "Boundary Layer Reference",
    ])
    ax3.set_ylabel("Temperature (K)")
    ax3.set_title("Engine Steady-State Temperatures")
    ax3.set_ylim(0, 1.2*tca.Tt)

    # Bottom plot - first 5 temperatures only
    ax4.plot(tca.x, jacket.T_cold)
    ax4.plot(tca.x, jacket.T_hot)
    ax4.plot(tca.x, cool.T)
    ax4.plot(tca.x, film.liquid.T)
    ax4.plot(tca.x, film.gas.T, '--')
    ax4.legend([
        "Cold Wall",
        "Hot Wall",
        "Regen Coolant",
        "Film Liquid",
        "Film Gas",
    ])
    ax4.set_xlabel("Distance from Injector (m)")
    ax4.set_ylabel("Temperature (K)")
    ax4.set_title("Wall & Film Temperatures")

    plt.tight_layout()

    #######################################################################################################
    ### Pump Impeller and Blades
    #######################################################################################################

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False)
    ax1.axis("equal")
    ax1.set_title("Meridional Plane Impeller Contours (mm)")

    # Ox pump on the left (mirrored in x, blue)
    ax1.plot(-ox_pump.impeller[0].shroud_curve[:, 0] * 1000.0, ox_pump.impeller[0].shroud_curve[:, 1] * 1000.0, color="blue")
    ax1.plot(-ox_pump.impeller[0].hub_curve[:, 0] * 1000.0, ox_pump.impeller[0].hub_curve[:, 1] * 1000.0, color="blue")
    ax1.plot(-ox_pump.impeller[0].blade.meanline_curve[:, 0] * 1000.0, ox_pump.impeller[0].blade.meanline_curve[:, 1] * 1000.0, "--", color="blue")

    # Fuel pump on the right (red)
    ax1.plot(fuel_pump.impeller[0].shroud_curve[:, 0] * 1000.0, fuel_pump.impeller[0].shroud_curve[:, 1] * 1000.0, color="red")
    ax1.plot(fuel_pump.impeller[0].hub_curve[:, 0] * 1000.0, fuel_pump.impeller[0].hub_curve[:, 1] * 1000.0, color="red")
    ax1.plot(fuel_pump.impeller[0].blade.meanline_curve[:, 0] * 1000.0, fuel_pump.impeller[0].blade.meanline_curve[:, 1] * 1000.0, "--", color="red")

    ax1.plot([0.0, 0.0], ax1.get_ylim(), color="green")
    ax1.plot(ax1.get_xlim(), [0.0, 0.0], color="green")
    ax1.set_xlabel("z (mm)")
    ax1.set_ylabel("r (mm)")
    ax1.margins(y=0.2)
 
    # Meanline Parameters below
    ax2_twin = ax2.twinx()
    p1, = ax2.plot(ox_pump.impeller[0].arc_meanline*1000, ox_pump.impeller[0].area_meanline*1E6, label="Ox Pump", color="blue")
    p2, = ax2.plot(fuel_pump.impeller[0].arc_meanline*1000, fuel_pump.impeller[0].area_meanline*1E6, label="Fuel Pump", color="red")
    p3, = ax2_twin.plot(ox_pump.impeller[0].arc_meanline*1000, ox_pump.impeller[0].blockage, label="Ox Pump Blockage", color="blue", linestyle="--")
    p4, = ax2_twin.plot(fuel_pump.impeller[0].arc_meanline*1000, fuel_pump.impeller[0].blockage, label="Fuel Pump Blockage", color="red", linestyle="--")
    ax2.set_xlabel("Meanline Arc Length (mm)")
    ax2.set_ylabel("Annular Flow Area (mm2)")
    ax2_twin.set_ylabel("Blockage")
    ax2.legend(handles=[p1, p2], labelcolor="linecolor")
    plt.tight_layout()

    #######################################################################################################
    ### Impeller Blades / Volute
    #######################################################################################################

    plt.figure(4)
    plt.clf()
    plt.grid(True)
    plt.axis("equal")
    plt.title("Impellers, Volutes, and Rotor (mm)")

    # Ox Pump Contours
    delta_angle = ox_pump.clocking * 2.0 * np.pi / ox_pump.impeller[0].blade_count  # Calculate the angle to rotate each blade

    for i in range(ox_pump.impeller[0].blade_count):  # for i = 0:(blade_count(1)-1)
        c = ox_pump.impeller[0].blade_curve
        rotated_curve = np.vstack([c[:,0]*np.cos(i*delta_angle) - c[:,1]*np.sin(i*delta_angle), 
                                c[:,0]*np.sin(i*delta_angle) + c[:,1]*np.cos(i*delta_angle)]).T
        plt.plot(
            rotated_curve[:, 0] * 1000.0,
            rotated_curve[:, 1] * 1000.0,
            linewidth=2,
            color="blue",
        )
        plt.plot(np.nan, np.nan)

    # plt.plot(
    #     ox_pump.volute_curve[:, 0, 0] * 1000.0 - 150.0,
    #     ox_pump.volute_curve[:, 1, 0] * 1000.0,
    #     color="blue",
    # )

    # Fuel Pump Contours
    delta_angle = fuel_pump.clocking * 2.0 * np.pi / fuel_pump.impeller[0].blade_count  # Calculate the angle to rotate each blade

    for i in range(fuel_pump.impeller[0].blade_count):
        c = fuel_pump.impeller[0].blade_curve
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

    #######################################################################################################
    ### Turbine Contours
    #######################################################################################################

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

    text_str = f"Shaft Speed = {turbine.shaft_speed*30/np.pi:0.2f} rpm\nGG = {gg.mdot:0.2g} kg/s, {gg.mdot*1790:0.0f} SCFM N2 ({gg.mdot/engine.mdot*100:.2g}%)"

    ax.text(0.2,0.1,text_str,transform=ax.transAxes,bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    #######################################################################################################
    ### Turbine Blades
    #######################################################################################################

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

    #######################################################################################################
    ### Inducer 3D Blade Geometry
    #######################################################################################################

    fig6 = plt.figure(6)
    ax1 = fig6.add_subplot(121, projection='3d')
    ax2 = fig6.add_subplot(122, projection='3d')
    ax1.set_title("Ox Inducer Geometry (mm)")
    ax2.set_title("Fuel Inducer Geometry (mm)")
    ax1.set_aspect("equal")
    ax2.set_aspect("equal")

    def plot_inducer_blades(ax, pump_obj, color):
        ind = pump_obj.inducer
        n_sl = ind.x_bladeline.shape[0] # number of bladeline points
        for k in range(ind.blade_number):
            phi = 2 * np.pi * k / ind.blade_number # rad - azimuth
            xr = ind.x_bladeline * np.cos(phi) - ind.y_bladeline * np.sin(phi)
            yr = ind.x_bladeline * np.sin(phi) + ind.y_bladeline * np.cos(phi)
            zr = ind.z_bladeline
            for n in range(n_sl):
                ax.plot(xr[n] * 1000, yr[n] * 1000, zr[n] * 1000, color=color, linewidth=0.8) # plotting bladelines
            ax.plot(xr[:, 0] * 1000, yr[:, 0] * 1000, zr[:, 0] * 1000, color=color, linewidth=1.5) # plotting leading edges
            ax.plot(xr[:, -1] * 1000, yr[:, -1] * 1000, zr[:, -1] * 1000, color=color, linewidth=1.5) # plotting trailing edges

    plot_inducer_blades(ax1, ox_pump, color="blue")
    plot_inducer_blades(ax2, fuel_pump, color="red")
    ax1.set_box_aspect([np.ptp(ax1.get_xlim()), np.ptp(ax1.get_ylim()), np.ptp(ax1.get_zlim())])
    ax2.set_box_aspect([np.ptp(ax2.get_xlim()), np.ptp(ax2.get_ylim()), np.ptp(ax2.get_zlim())])
    
    text1 = f"phi = {ox_pump.inducer.flow_coeff:0.4g} \npsi = {ox_pump.inducer.head_coeff:0.3g}, \nωss = {ox_pump.inducer.suction_specific_speed_available:0.4g} (Nss = {2733.00*ox_pump.inducer.suction_specific_speed_available:0.0f})\nShaft Speed = {30/np.pi*ox_pump.shaft_speed:0.0f} RPM"
    ax1.text2D(0, -0.1, text1, transform=ax1.transAxes, bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    text2 = f"phi = {fuel_pump.inducer.flow_coeff:0.4g} \npsi = {fuel_pump.inducer.head_coeff:0.3g}, \nωss = {fuel_pump.inducer.suction_specific_speed_available:0.4g} (Nss = {2733.00*fuel_pump.inducer.suction_specific_speed_available:0.0f})\nShaft Speed = {30/np.pi*fuel_pump.shaft_speed:0.0f} RPM"
    ax2.text2D(0, -0.1, text2, transform=ax2.transAxes, bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    #######################################################################################################
    ### Volute Contours
    #######################################################################################################

    fig7 = plt.figure(7)
    ax1 = fig7.add_subplot(121, projection='3d')
    ax2 = fig7.add_subplot(122, projection='3d')
    ax1.set_title("Ox Volute Geometry (mm)")
    ax2.set_title("Fuel Volute Geometry (mm)")


    def plot_volute(ax, pump_obj, color):
        vol = pump_obj.volute
        imp = pump_obj.impeller[0]

        X = vol.contour[:, :, 0]
        Y = vol.contour[:, :, 1]
        Z = vol.contour[:, :, 2]

        # Volute tube surface
        ax.plot_surface(X * 1000, Y * 1000, Z * 1000, color=color, alpha=0.6, linewidth=0, antialiased=True)

        # Revolve impeller shroud + hub for context
        theta_rev = np.linspace(0.0, 2.0 * np.pi, 64)  # rad
        for curve, col in [(imp.shroud_curve, "0.5"), (imp.hub_curve, "0.7")]:
            r = curve[:, 0][:, None]
            z = curve[:, 1][:, None]
            Xr = r * np.cos(theta_rev)[None, :] * 1000
            Yr = r * np.sin(theta_rev)[None, :] * 1000
            Zr = np.broadcast_to(z, Xr.shape) * 1000
            ax.plot_surface(Xr, Yr, Zr, color=col, alpha=0.25, linewidth=0)

        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_zlabel("z (mm)")

    plot_volute(ax1, ox_pump, color="blue")
    plot_volute(ax2, fuel_pump, color="red")
    ax1.set_box_aspect([np.ptp(ax1.get_xlim()), np.ptp(ax1.get_ylim()), np.ptp(ax1.get_zlim())])
    ax2.set_box_aspect([np.ptp(ax2.get_xlim()), np.ptp(ax2.get_ylim()), np.ptp(ax2.get_zlim())])

    text1 = f"method = {ox_pump.volute.design_method}\nd_throat = {ox_pump.volute.d_throat*1000:0.2f} mm\nc_throat = {ox_pump.volute.c_throat:0.1f} m/s\nloss = {ox_pump.volute.total_loss:0.4g}"
    ax1.text2D(0, -0.1, text1, transform=ax1.transAxes, bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    text2 = f"method = {fuel_pump.volute.design_method}\nd_throat = {fuel_pump.volute.d_throat*1000:0.2f} mm\nc_throat = {fuel_pump.volute.c_throat:0.1f} m/s\nloss = {fuel_pump.volute.total_loss:0.4g}"
    ax2.text2D(0, -0.1, text2, transform=ax2.transAxes, bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))



    #######################################################################################################
    ### Inducer 2D Blade Geometry
    #######################################################################################################

    # fig7, (ax_ox, ax_fuel) = plt.subplots(2, 1, figsize=(10, 10))

    # # Ox inducer - all streamlines, blue gradient (light at hub → dark at tip)
    # n_sl = ox_pump.inducer.r_bladeline.shape[0]
    # cmap_ox = plt.cm.Blues
    # for n in range(n_sl):
    #     color_intensity = 0.3 + 0.7 * (n / (n_sl - 1))
    #     rtheta = ox_pump.inducer.r_bladeline[n, :] * ox_pump.inducer.theta_bladeline[n, :] * 1000
    #     z      = ox_pump.inducer.z_bladeline[n, :] * 1000
    #     ax_ox.plot(rtheta, z, color=cmap_ox(color_intensity), linewidth=0.8)

    # ax_ox.set_xlabel("r·θ (mm)")
    # ax_ox.set_ylabel("z (mm)")
    # ax_ox.set_aspect("equal", adjustable="datalim")
    # ax_ox.grid(True)
    # ax_ox.set_title(
    #     f"Ox Inducer – Unrolled Streamlines  "
    #     f"(β_LE: hub={np.rad2deg(ox_pump.inducer.blade_angle[0,0]):.1f}°, "
    #     f"tip={np.rad2deg(ox_pump.inducer.blade_angle[-1,0]):.1f}°)   "
    #     f"[light = hub, dark = tip]"
    # )

    # # Fuel inducer - all streamlines, red gradient (light at hub → dark at tip)
    # n_sl = fuel_pump.inducer.r_bladeline.shape[0]
    # cmap_fuel = plt.cm.Reds
    # for n in range(n_sl):
    #     color_intensity = 0.3 + 0.7 * (n / (n_sl - 1))
    #     rtheta = fuel_pump.inducer.r_bladeline[n, :] * fuel_pump.inducer.theta_bladeline[n, :] * 1000
    #     z      = fuel_pump.inducer.z_bladeline[n, :] * 1000
    #     ax_fuel.plot(rtheta, z, color=cmap_fuel(color_intensity), linewidth=0.8)

    # ax_fuel.set_xlabel("r·θ (mm)")
    # ax_fuel.set_ylabel("z (mm)")
    # ax_fuel.set_aspect("equal", adjustable="datalim")
    # ax_fuel.grid(True)
    # ax_fuel.set_title(
    #     f"Fuel Inducer – Unrolled Streamlines  "
    #     f"(β_LE: hub={np.rad2deg(fuel_pump.inducer.blade_angle[0,0]):.1f}°, "
    #     f"tip={np.rad2deg(fuel_pump.inducer.blade_angle[-1,0]):.1f}°)   "
    #     f"[light = hub, dark = tip]"
    # )

    # plt.tight_layout()

    # plt.figure()
    # n_pts = ox_pump.inducer.blade_angle.shape[1]
    # s = np.linspace(0, 1, n_pts)
    # for n in range(len(ox_pump.inducer.blade_angle)):
    #     plt.plot(s, np.rad2deg(ox_pump.inducer.blade_angle[n]), label=f'n={n}')
    # plt.xlabel('s/L')
    # plt.ylabel('β (deg)')
    # plt.grid(True)

    ###################################################################################################################

    print(f"ox_pump.shaft_power: {ox_pump.shaft_power}")
    print(f"fuel_pump.shaft_power: {fuel_pump.shaft_power}")
    print(f"ox_pump.shaft_speed: {ox_pump.shaft_speed}")
    print(f"fuel_pump.shaft_speed: {fuel_pump.shaft_speed}")
    
    plt.show(block=False)
    input("Press Enter to exit...")