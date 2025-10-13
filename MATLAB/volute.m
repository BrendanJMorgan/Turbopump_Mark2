%% Pump Handbook page 2.68

r_tongue(p) = r_outlet_impeller(p)*(1.03+0.0025*specific_speed(p)+7E-8*density_pump(p)*head_pump(p)); % m - Gulich table 10.2, pg 729
t_tongue(p) = 0.04*r_outlet_impeller(p); % m - Gulich pg 493
A_throat(p) = ( vdot_pump(p)/(3.8*r_outlet_impeller(p)*v_tangential_outlet(p)) * ( 1 + sqrt(1+15.2*r_outlet_impeller(p)*v_tangential_outlet(p)*(r_tongue(p)+t_tongue(p))/vdot_pump(p)) ) )^2; % m2
% where does this come from ^ 
r_throat(p) = r_outlet_impeller(p) + 0.5*sqrt(A_throat(p)) + t_tongue(p); % m

theta_volute(:,p) = linspace(0, 2*pi, 1000)'; % rad
A_volute(:,p) = ( vdot_pump(p).*theta_volute(:,p)./(7.6*pi*r_outlet_impeller(p)*v_tangential_outlet(p)) .* ( 1 + sqrt(1+30.4*pi*r_outlet_impeller(p)*r_tongue(p)*v_tangential_outlet(p)./(vdot_pump(p)*theta_volute(:,p))) ) ).^2; % m2
A_volute(isnan(A_volute)) = 0;

h_volute(p) = sqrt(0.25*t_tongue(p)^2+A_volute(end,p)) - 0.5*t_tongue(p); % m - height of volute gap - set so the outlet is a square
r_volute(:,p) = r_tongue(p) + A_volute(:,p)/h_volute(p); % m - radius of volute wall
OD_volute_outlet(p) = 2*sqrt(A_volute(p,end)/pi); % m - diameter of volute outlet hole
r_volute_outlet(p) = r_volute(end,p) - h_volute(end,p)/2; % m

v_throat(p) = vdot_pump(p)./A_throat(p); % m/s - fluid outlet velocity into plumbing

volute_curve(:,:,p) = r_volute(:,p).*[cos(clocking(p)*theta_volute(:,p)), sin(clocking(p)*theta_volute(:,p))]; % [m, m]

%%

% "...considerable freedom to configure the cross section without risking major losses ... flat cross sections result in less intense secondary flow than circular cross
% sections, thereby generating fewer losses ... a ratio of width B to height H in the range of B/H = 2 to 3 may be considered the  optimum." - Gulich pg 492

