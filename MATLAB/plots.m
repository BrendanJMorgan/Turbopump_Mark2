%% TCA (thrust chamber assembly) Contours
figure(1); clf
colororder('default')
plot(x,r1,x,r2,x,-1*r1,x,-1*r2, 'color','blue');
axis equal
xlabel("Distance from Injector (m)");
title("Combustion Chamber Contours")

%% TCA Temperatures 
figure(2); clf
plot(x,p_cool/1E5)
% plot(x,T_wall_cold,x,T_wall_hot,x,T_cool,x,T_film,x,T_free,x,T_ab,x,T_recovery,x,T_ref)
% yline(0)
% legend("Cold Wall","Hot Wall","Regen Coolant","Film Coolant","Free-Stream Gas","Adiabatic", "Recovery", "Gas Property Reference",'Location','northeast');
xlabel("Distance from Injector (m)");
% ylabel("Temperature (K)");
% title("Engine Steady-State Temperatures")

%% Pump Impeller and Blades
figure(3); clf; hold on; grid on; axis equal; title("Impeller and Shroud Contours (mm)")
plot(-shroud_curve(:,1,1)*1000, shroud_curve(:,2,1)*1000, 'color', 'cyan')
plot(-impeller_curve(:,1,1)*1000, impeller_curve(:,2,1)*1000, 'color', 'cyan')
plot(-meanline_curve(:,1,1)*1000, meanline_curve(:,2,1)*1000, '--', 'color', 'cyan');
plot(shroud_curve(:,1,2)*1000, shroud_curve(:,2,2)*1000, 'color', 'red')
plot(impeller_curve(:,1,2)*1000, impeller_curve(:,2,2)*1000, 'color', 'red')
plot(meanline_curve(:,1,2)*1000, meanline_curve(:,2,2)*1000, '--', 'color', 'red');
line([0 0], ylim, 'color', 'green');  % x-axis
line(xlim, [0 0], 'color', 'green');  % y-axis

%% Impeller Blades / Volute
figure(4); clf; hold on; grid on; axis equal; title('Impellers, Volutes, and Rotor (mm)');

% Ox Pump Contours
delta_angle=clocking(1)*2*pi/blade_count(1); % Calculate the angle to rotate each blade
for i = 0:(blade_count(1)-1) 
    rotation_matrix = [cos(i * delta_angle), -sin(i * delta_angle); sin(i * delta_angle), cos(i * delta_angle)];
    rotated_curve = blade_curve(:,1:2,1) * rotation_matrix';
    plot(rotated_curve(:, 1)*1000-150, rotated_curve(:, 2)*1000, 'LineWidth', 2, 'color', 'cyan'); 
    plot(NaN, NaN); % Prevent connection between different blades
end
plot(volute_curve(:,1,1)*1000-150, volute_curve(:,2,1)*1000, 'color', 'cyan');

% Fuel Pump Contours
delta_angle = clocking(2)*2*pi/blade_count(2); % Calculate the angle to rotate each blade
for i = 0:(blade_count(2)-1)
    rotation_matrix = [cos(i * delta_angle), -sin(i * delta_angle); sin(i * delta_angle), cos(i * delta_angle)];
    rotated_curve = blade_curve(:,1:2,2) * rotation_matrix';
    plot(rotated_curve(:, 1)*1000, rotated_curve(:, 2)*1000, 'LineWidth', 2, 'color', 'red'); 
    plot(NaN, NaN); % Prevent connection between different blades
end
plot(volute_curve(:,1,2)*1000, volute_curve(:,2,2)*1000, 'color', 'red');

% Turbine Contours
theta = linspace(0, 2*pi, 1000)'; % rad
plot(r_pitchline_rotor*1000*cos(theta), r_pitchline_rotor*1000*sin(theta), 'LineStyle', '--', 'color', '#ffA500')
plot(r_tip_rotor*1000*cos(theta), r_tip_rotor*1000*sin(theta), r_base_rotor*1000*cos(theta), r_base_rotor*1000*sin(theta), 'color', '#ffA500')

str = sprintf('GG = %.2g kg/s, %.0f SCFM (%.2g%%)',mdot_gg, mdot_gg*1790, mdot_gg/mdot_total*100);
annotation('textbox',[.2 .5 .3 .3],'String',str,'FitBoxToText','on');

%% Turbine Blades
figure(5), clf, hold on, axis equal, grid on, title('Rotor Blade Geometry (mm)');

% Trailing edge
sagitta        = blade_chord_rotor/2*(1-cos(blade_angle_rotor))/sin(blade_angle_rotor); % m
tip            = [-sagitta, blade_chord_rotor/2]; % [m,m] - chord point
blade_radius   = sagitta/2 + blade_chord_rotor^2/(8*sagitta); % m       
center_arc     = [-blade_radius, 0]; % y=0
angle_trailing = atan2(tip(2)-0, tip(1)+blade_radius); % rad             
arc_theta      = linspace(-angle_trailing, angle_trailing, 120)';                       
x_trailing     = center_arc(1)+blade_radius*cos(arc_theta); % [m,m]
y_trailing     = center_arc(2)+blade_radius*sin(arc_theta); % [m,m]

% Leading edge
apex_leading        = [blade_pitch_rotor-blade_opening_rotor, 0]; % [m,m] - location of blade apex point             
center_leading      = apex_leading + [-radius_leading, 0]; % [m,m] on y=0
center_tip          = tip - center_leading; % m                                  
base_angle          = atan2(center_tip(2), center_tip(1)); % rad                            
tan_angle           = base_angle - acos(radius_leading/norm(center_tip)); % rad - upper tangent
tangency_leading    = center_leading + radius_leading*[cos(tan_angle), sin(tan_angle)]; % [m,m]
theta_leading       = linspace(-tan_angle, tan_angle,60); % rad                      
x_leading           = center_leading(1)+radius_leading*cos(theta_leading); % [m,m]        
y_leading           = center_leading(2)+radius_leading*sin(theta_leading); % [m,m]


for off = [-blade_pitch_rotor*1000, 0, blade_pitch_rotor*1000]
    plot(x_trailing*1000+off, y_trailing*1000, 'y','LineWidth',2) % trailing edge arc
    plot(x_leading*1000+off, y_leading*1000, 'y','LineWidth',2); % leading edge arc
    plot([tangency_leading(1)*1000+off,tip(1)*1000+off],[tangency_leading(2)*1000,tip(2)*1000], 'y','LineWidth',2) % leading edge tangent line
    plot([tangency_leading(1)*1000+off,tip(1)*1000+off],[-tangency_leading(2)*1000,-tip(2)*1000], 'y','LineWidth',2) % leading edge tangent line
end

str = sprintf('Blade Count = %g ', blade_count_rotor);
annotation('textbox',[.2 .5 .3 .3],'String',str,'FitBoxToText','on');


%% Impeller Meanlines

figure(6), clf, hold on, grid on, title('Impeller Meanline'), xlabel('Radius (mm)'), ylabel('Velocity (m/s)')

p=1;
r_mm = 1000*meanline_curve(:,1,p);
plot(r_mm,u_blade(:,p), '--', 'color', 'cyan')
plot(r_mm,v_tangential(:,p), ':', 'color', 'cyan')
plot(r_mm,v_merid(:,p), '-.', 'color', 'cyan')
plot(r_mm,u_blade(:,p).*(1-slip_factor(:,p)) - v_tangential(:,p), 'color', 'cyan');

legend('Solid Body Velocity', 'Fluid Tangential Velocity', 'Meridional (Numerator) Velocity', 'Denominator Velocity')

figure(7), clf, hold on, grid on, title('Blade Meanline'), xlabel('Arc Length'), ylabel('Angle (deg)')
plot(blade_azimuth(:,1)*180/pi, 'color', 'cyan')
hold on
plot(blade_azimuth(:,2)*180/pi, 'color', 'red')









