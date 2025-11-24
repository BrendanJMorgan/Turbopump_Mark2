%% Impeller Contours - pg 2.29-2.32 in pump handbook

r_eye_margin(p) = 1.15 + 0.2116*(specific_speed(p)-0.2836); % unitless - margin correction factor, Gulich 7.1.4; for a normal impeller (not a suction impeller)
swirl_number(p) = 1; % NEEDS TO BE A FUNCTION OF INDUCER OUTLET CONDITIONS (ITERATIVE)
% r_inlet_impeller = 0.5 * r_eye_margin * sqrt(r_hub + 1.48E-3 * outlet_flow_coeff(p)/2* (52.9*specific_speed(p))^1.33 / swirl_number^0.67 );
r_inlet_impeller(p) = r_eye_margin(p) * sqrt(r_hub(p)^2 + 0.0026418 * outlet_flow_coeff(p)* specific_speed(p)^1.33 / swirl_number(p)^0.67 );

r_outlet_impeller(p) = 1 / abs(pump_shaft_speed(p)) * sqrt (g*head_pump(p) / head_coeff(p)); % m - impeller exit radius, Gulich eqn 7.1.3

%% The exit width b2
% w_exit(p) = vdot_pump(p) / (2*pi*pump_shaft_speed(p)*r_outlet_impeller(p)^2*outlet_flow_coeff(p)*blockage(p));
w_exit(p) = 2*r_outlet_impeller(p) * (0.017 + 0.1386*specific_speed(p) - 0.022387*specific_speed(p)^2 + 0.0013767*specific_speed(p)^3); % m - outlet width; empirical; Gulich eqn 7.1

% this needs boundary layer and leakage effect - pg 2.60?
% could parameterize a width at each station with a blockage vector, linear
% ramp of area?kk

%% Hub and shroud profiles

r_min(p) = 0.5*r_inlet_impeller(p); % m - minimum allowable radius of curvature
bezier_params_initial = [0.01, 0.01, 0.01]; % starting values
cost_function = @(bezier_params) cost_function_impeller(bezier_params, r_inlet_impeller(p), r_outlet_impeller(p), w_exit(p), impeller_thickness(p), r_hub(p), r_min(p));
bezier_params_optimal(:,p) = fminsearch(cost_function, bezier_params_initial);

bezier_vert_optimal(p) = bezier_params_optimal(1,p);
bezier_horiz_optimal(p) = bezier_params_optimal(2,p);
impeller_height_optimal(p) = bezier_params_optimal(3,p);

[~, ~, min_shroud_curvature(:,p), ~, shroud_curve(:,:,p), impeller_curve(:,:,p), control_points(:,:,p)] = ... 
    compute_curves(bezier_vert_optimal(p), bezier_horiz_optimal(p), r_inlet_impeller(p), r_outlet_impeller(p), w_exit(p), impeller_thickness(p), impeller_height_optimal(p), r_hub(p));

min_radius(p) = 1/max(min_shroud_curvature(:,p));
if min_radius(p) < 0.99*r_min(p)    
	fprintf("\nPump shroud curvature reaches %.2f mm radius. Need above ~%.2f mm (half of eye radius) to avoid cavitation.", min_radius*1000, r_min*1000)
end

%% Functions

% Function to compute slopes and minimum curvature for shroud and impeller curves in a pump
function [slope_shroud, slope_impeller, min_shroud_curvature, flatness_impeller, shroud_curve, impeller_curve, control_points] = compute_curves(bezier_vert, bezier_horiz, r_inlet_impeller, r_exit_impeller, w_exit, impeller_thickness, impeller_height, r_hub)

    % Generate the shroud curve using bezier control points
    control_points = [r_inlet_impeller, impeller_height; r_inlet_impeller, impeller_height - bezier_vert; r_exit_impeller - bezier_horiz, impeller_thickness + w_exit; r_exit_impeller, impeller_thickness + w_exit];
    shroud_curve = bezier(control_points, 1000); % [m,m]

    % Compute the flow areas at eye and exit
    A_pump_eye = pi * (r_inlet_impeller^2 - r_hub^2); % m2
    A_pump_exit = 2 * pi * r_exit_impeller * w_exit; % m2

    % Generate area array from eye to exit along the shroud curve
    % ARBITRARY ASSUMPTION TO MAKE IT LINEAR
    A_pump_flow = linspace(A_pump_eye, A_pump_exit, length(shroud_curve))'; % m2

    % Compute the angles of the normals at each point on the shroud curve
    normal_angles = atan2(gradient(shroud_curve(:,2)), gradient(shroud_curve(:,1))) + pi/2; % rad
    normal_angles(normal_angles > 2 * pi) = normal_angles(normal_angles > 2 * pi) - 2 * pi; % rad

    % Compute the crosswise gap between shroud and impeller curves at each point
    crosswise_gap = sec(normal_angles) .* (shroud_curve(:,1) - sqrt(pi^2 * shroud_curve(:,1).^2 - pi * A_pump_flow .* cos(normal_angles)) / pi); % m

    % Generate impeller curve based on the shroud curve and the crosswise gap
    impeller_curve = [shroud_curve(:,1) - crosswise_gap .* cos(normal_angles), shroud_curve(:,2) - crosswise_gap .* sin(normal_angles)]; % [m,m]

    % Calculate the slope of the shroud and impeller curves at their ends
    slope_shroud = (shroud_curve(end, 2) - shroud_curve(end - 1, 2)) / (shroud_curve(end, 1) - shroud_curve(end - 1, 1));
    slope_impeller = (impeller_curve(end, 2) - impeller_curve(end - 1, 2)) / (impeller_curve(end, 1) - impeller_curve(end - 1, 1));

    % Calculate flatness, i.e. the max height dip below the impeller exit plane
    flatness_impeller = min(impeller_curve(:,2)) - impeller_curve(end, 2); % m

    % Compute the minimum curvature of the shroud curve
    s = [0; cumsum(sqrt(diff(shroud_curve(:, 1)).^2 + diff(shroud_curve(:, 2)).^2))];
    min_shroud_curvature = abs(gradient(shroud_curve(:, 1), s) .* gradient(gradient(shroud_curve(:, 2), s), s) - gradient(shroud_curve(:, 2), s) .* gradient(gradient(shroud_curve(:, 1), s), s)) ./ (gradient(shroud_curve(:, 1), s).^2 + gradient(shroud_curve(:, 2), s).^2).^(3/2);
end

function cost = cost_function_impeller(bezier_params, r_inlet_impeller, r_exit_impeller, w_exit, impeller_thickness, r_hub, r_min)
    bezier_vert = bezier_params(1);
    bezier_horiz = bezier_params(2);
    impeller_height = bezier_params(3);
    [slope_shroud, slope_impeller, min_shroud_curvature, flatness_impeller] = compute_curves(bezier_vert, bezier_horiz, r_inlet_impeller, r_exit_impeller, w_exit, impeller_thickness, impeller_height, r_hub);
    
    if flatness_impeller > 0
        flatness_impeller = 0;
    end

    % The cost is the sum of the absolute slopes (want them to be 0 for a horizontal line) and the curvature
    cost = 100*abs(slope_shroud) + 100*abs(slope_impeller) + 1000*abs(1/max(min_shroud_curvature)-r_min) + 1000*abs(flatness_impeller);
end
