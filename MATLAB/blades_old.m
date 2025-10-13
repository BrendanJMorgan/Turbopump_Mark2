%% Impeller Blades

%% Leading and Trailing Edges
r_inlet(p) = 0.5*r_eye_impeller(p); % m - radial location of blade inlets - typically vertically aligned with eye radius - this assumption can be changed

%% Surface of Revolution
rms_curve(:,:,p) = 0.5 * (shroud_curve(:,:,p) + impeller_curve(:,:,p)); % [m,m] - halfway between shroud and impeller

%% Inlet Blade Angles
inlet_index(p) = find(rms_curve(:,1) > r_inlet, 1); % unitless - point at which rms curve lies directly beneath eye radius
inlet_gap(p) = norm( shroud_curve(inlet_index(p), :) - impeller_curve(inlet_index, :) ); % m - inlet gap width in vicinity of eye radius
u_blade_inlet(p) = pump_shaft_speed(p)*r_inlet(p); % m/s - rotational speed at the inlet (U_1 in pump handbook)
v_merid_inlet(p) = vdot_pump(p)/(2*pi*r_inlet(p)*inlet_gap(p)); % m/s - fluid velocity at the inlet, assumed to be entirely radial. Can be found with CFD instead for more accuracy (V_1 in pump handbook)
v_tangential_inlet(p) = 0; % m/s - assuming no pre-swirl from the inducer
blade_angle_inlet(p) = atan2(v_merid_inlet(p),u_blade_inlet(p)); % rad - inlet angle of blades relative to tangential azimuth  (beta_1,f in pump handbook)

%% Outlet Blade Angles
slip_factor(p) = 0.15; % Buseman form "typically between 0.1 and 0.2 ... for frictionless flow ... with log-spiral blades ... and a 2D, radial-flow geometry with parallel hub and shroud"
u_blade_exit(p) = pump_shaft_speed(p)*r_exit_impeller(p);	% m/s - blade tangential velocity at exit (U_2 in pump handbook)
v_slip(p) = slip_factor(p)*u_blade_exit(p);			        % m/s - slip velocity (V_s in pump handbook)
v_merid_exit(p) = vdot_pump(p)/(2*pi*r_exit_impeller(p)*w_exit(p));	% m/s - meridional velocity (V_m,2 in pump handbook)                                                       

hydraulic_efficiency(p) = 1 - 0.071 / vdot_pump(p)^0.25; % unitless - Jekat's Empirical Formula - valid for all specific speeds

v_tangential_exit(p) = g*head_pump(p) / (hydraulic_efficiency(p)*u_blade_exit(p)) + v_tangential_inlet(p); % m/s - this is Euler's pump equation rearranged (V_theta,2 in pump handbook). The tangential flow speed to achieve head requirement.                                                                                                                                                                                                                                                                                                                                                                                                                           
blade_angle_outlet(p) = atan(v_merid_exit(p) / (u_blade_exit(p)*(1-slip_factor(p)) - v_tangential_exit(p)));    % rad - angle between blade tip azimuth and local tangential azimuth (beta_2 in pump handbook)
if blade_angle_outlet(p) <= 0
	blade_angle_outlet*180/pi
	error("Cannot converge to solution with valid blade outlet angle")
end

%% Individual Blade Curve
params_initial = [0.01,0.01,-clocking(p)]; % starting values; control1, control2
cost_function = @(params) cost_function_blade(params, r_inlet(p), r_exit_impeller(p), blade_angle_inlet(p), blade_angle_outlet(p));
blade_params_optimal(:,p) = fminsearch(cost_function, params_initial);

control1 = blade_params_optimal(1); % m - length between inlet point and inlet control point
control2 = blade_params_optimal(2); % m - length between outlet point and outlet control point
sweep(p) = blade_params_optimal(3); % rad
[blade_curve(:,:,p), blade_control_points(:,:,p)] = compute_blade_curve(control1, control2, sweep(p), r_inlet(p), r_exit_impeller(p), blade_angle_inlet(p), blade_angle_outlet(p));

%% Number of Blades
blade_z(:,p) = interp1(rms_curve(:,1,p), rms_curve(:,2,p), blade_curve(:,1,p)); % m - vertical position component of blade curve
blade_arc_length(p) = sum( vecnorm( diff([blade_curve(:,1,p), blade_curve(:,2,p), blade_z(:,p)])' ) ); % m - arc length of one individual blade, measure at rms line
solidity(p) = interp1([0, 0.4, 3], [1.8, 1.8, 1], specific_speed(p)); % solidity is the optimal ratio of blade chord to blade spacing. Pump handbook page 2.36 (sigma)
blade_number(p) = round(solidity(p) * 2*pi*r_exit_impeller(p) / blade_arc_length(p)); % number of blades, rounded to the nearest even number (n)

%% Compute Blockage (pump handbook page 2.63)
blade_thickness(p) = 0.04*r_exit_impeller(p); % m
boundary_layer_thickness(p) = 0.002*blade_arc_length(p); % m
meridional_length(p) = sum( vecnorm( diff([blade_curve(:,1), blade_z])' ) ); % m - arc length of blade projected onto the meridional plane
blockage(p) = 1 - (2*boundary_layer_thickness(p)+blade_thickness(p))/(0.0254*2*pi * meridional_length(p)/blade_arc_length(p)); 

%% Helper Functions
% Function to compute slopes and minimum curvature for blade curve in a pump
function [blade_curve, control_points] = compute_blade_curve(control1, control2, sweep, radius1, radius2, polar1, polar2)
	inlet_point = [radius1,0];																% [m,m] - blade inner edge starting location
	outlet_point = radius2*[cos(sweep), sin(sweep)];										% [m,m] - blade outer edge ending location
	inlet_azimuth = atan2(inlet_point(2), inlet_point(1)) - (pi/2 - polar1);			    % rad - blade inner edge absolute direction
	outlet_azimuth = atan2(outlet_point(2), outlet_point(1)) - (pi/2 - polar2);		        % rad - blade outer edge absolute direction
	inlet_control = control1*[cos(inlet_azimuth),sin(inlet_azimuth)] + inlet_point;			% [m,m] - control point for the inlet edge
	outlet_control = -control2*[cos(outlet_azimuth),sin(outlet_azimuth)] + outlet_point;	% [m,m] - control point for the outlet edge

	% Generate the shroud curve using four bezier control points
    control_points = [inlet_point;inlet_control;outlet_control;outlet_point];
    blade_curve = bezier(control_points); % [m,m]

	% control_error1 = atan2(blade_curve(2,2)-inlet_point(2), blade_curve(2,1)-inlet_point(1)) - inlet_azimuth;
	% control_error2 = atan2(blade_curve(end-1,2)-outlet_point(2), blade_curve(end-1,1)-outlet_point(1)) - outlet_azimuth;
end

function cost = cost_function_blade(params, radius1, radius2, polar1, polar2)
    control1 = params(1);   % m - length between inlet point and inlet control point
	control2 = params(2);   % m - length between outlet point and outlet control point
    sweep = params(3);      % rad
    [blade_curve, ~] = compute_blade_curve(control1, control2, sweep, radius1, radius2, polar1, polar2);

	% blade_arc_length = sum(sqrt(sum(diff(blade_curve).^2, 2))); % m - arc length of one individual blade
    s = [0; cumsum(sqrt(diff(blade_curve(:, 1)).^2 + diff(blade_curve(:, 2)).^2))];
    blade_curvature = abs(gradient(blade_curve(:, 1), s) .* gradient(gradient(blade_curve(:, 2), s), s) - gradient(blade_curve(:, 2), s) ...
		.* gradient(gradient(blade_curve(:, 1), s), s)) ./ (gradient(blade_curve(:, 1), s).^2 + gradient(blade_curve(:, 2), s).^2).^(3/2); % 1/m
	max_blade_curvature = max(blade_curvature); % 1/m
    
	cost = 1*abs(max_blade_curvature);
end







