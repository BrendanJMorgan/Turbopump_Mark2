%% Impeller Blades

%% Meanline Definition
r_inlet_blade(p) = r_inlet_impeller(p); % m - radial location of blade inlets - typically vertically aligned with eye radius - this assumption can be changed
meanline_curve_untrimmed(:,:,p) = 0.5 * (shroud_curve(:,:,p) + impeller_curve(:,:,p)); % [m,m] - halfway between shroud and impeller
inlet_index(p) = find(meanline_curve_untrimmed(:,1,p) > r_inlet_blade, 1); % unitless - point at which meanline curve lies directly beneath eye radius

r_meanline = linspace(r_inlet_blade(p), r_outlet_impeller(p), 1000)'; % m
z_meanline = pchip(meanline_curve_untrimmed(:,1,p), meanline_curve_untrimmed(:,2,p), r_meanline); % m
meanline_curve(:,:,p) = [r_meanline, z_meanline]; % [m,m] - blade meanline

% meanline_curve_blade =  meanline_curve(inlet_index(p):end,:,p); % [m,m] the portion of the meanline with a blade

%% Velocity Triangles
u_blade(:,p) = pump_shaft_speed(p)*meanline_curve(:,1,p); % m/s - blade tangential velocity
v_merid(:,p) = vdot_pump(p)/(pi*(r_inlet_impeller(p)^2-r_hub(p)^2)).*ones(1000,1); % m/s - NEEDS TO VARY ALONG MEANLINE
% v_merid(:,p) = vdot_pump(p) ./ (2*pi*meanline_curve(:,1,p).*vecnorm(shroud_curve(:,:,p)-impeller_curve(:,:,p),2,2)); % m/s - meridional velocity

hydraulic_efficiency(p) = 1 - 0.071 / vdot_pump(p)^0.25; % unitless - Jekat's Empirical Formula - valid for all specific speeds CONVERT TO ANDERSON?

v_tangential_inlet(p) = 0; % m/s - assuming no pre-swirl from the inducer
v_tangential_outlet(p) = g*head_pump(p) / (hydraulic_efficiency(p)*u_blade(end,p)) + v_tangential_inlet(p); 
    % m/s - this is Euler's pump equation rearranged (V_theta,2 in pump handbook). The tangential flow speed to achieve head requirement.                                                                                                                                                                                                                                                                                                                                                                                                                           
v_tangential(:,p) = linspace(v_tangential_inlet(p), v_tangential_outlet(p), length(meanline_curve(:,p)))'; % assuming linear ramp of fluid tangential velocity

fluid_azimuth(:,p) = atan(v_merid(:,p) ./ (u_blade(:,p) - v_tangential(:,p))); % rad - angle between fluid flow azimuth and local tangential azimuth (beta_f2 in pump handbook)

%% Camber Line Integration - Blade Curve Construction
blade_azimuth(:,p) = atan2(v_merid(:,p), (u_blade(:,p).*(1-slip_factor(:,p)) - v_tangential(:,p)));    % rad - angle between blade tip azimuth and local tangential azimuth (beta_2 in pump handbook)

dm = sqrt( diff(r_meanline).^2 + diff(z_meanline).^2 ); % meanline differential arc length
meanline_arc_length(:,p) = [0;cumsum(dm)]; % meanline cumulative arc length
theta_blade(:,p) = cumtrapz(r_meanline,cot(blade_azimuth(:,p))./r_meanline); % rad

blade_curve(:,:,p) = [meanline_curve(:,1,p).*[cos(theta_blade(:,p)),-clocking(p)*sin(theta_blade(:,p))],meanline_curve(:,2,p)]; % [m,m,m]

%% Number of Blades
blade_arc_length(p) = sum(sqrt(diff(blade_curve(:,1,p)).^2+diff(blade_curve(:,2,p)).^2+diff(blade_curve(:,3,p)).^2)); % m - arc length of one individual blade, measure at meanline line
solidity_ideal(p) = interp1([0, 0.4, 3], [1.8, 1.8, 1], specific_speed(p)); % solidity is the optimal ratio of blade chord to blade spacing. Pump handbook page 2.36 (sigma)
blade_count(p) = round(solidity_ideal(p) * 2*pi*r_outlet_impeller(p) / blade_arc_length(p)); % number of blades
solidity(p) = blade_count(p)*blade_arc_length(p)/(2*pi*r_outlet_impeller(p)); % unitless

%% Compute Blockage (pump handbook page 2.63)
blade_thickness(p) = 0.04*r_outlet_impeller(p); % m
boundary_layer_thickness(p) = 0.002*blade_arc_length(p); % m
meridional_length(p) = sum(sqrt(diff(blade_curve(:,1,p)).^2+diff(blade_curve(:,3,p)).^2)); % m - arc length of blade projected onto the meridional plane
blockage(p) = 1 - solidity(p)/meanline_arc_length(end,p)*(boundary_layer_thickness(p)+blade_thickness(p)); % unitless - proportion of area not blocked by blades

%% Velocity Slip - Pfliederer's Method - pump handbook pg 2.35
% Extending this to the full arc length of the meanline, instead of just the tip
% a = 0.75;  % unitless - 0.65 to 0.85 for volute
% r_ref = r_meanline(1); % m - scalar reference radius
% psi_prime = a * (1 + sin(blade_azimuth(:,p))); % unitless
% slip_factor(:,p) = psi_prime .* r_meanline .^ 2 ./  (blade_count(p) * 0.5 * (r_meanline.^2 - r_ref^2)); % unitless

% = fraction of impeller discharge meridional area (that is, the area
% normal to the velocity component Vm, 2) that is not blocked by the
% thickness of the blades and the boundary layer displacement
% thickness on blades and on hub and shroud surfaces

%% Cavitation Check
% eqns 25 and 26 in pump handbook







