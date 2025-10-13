%% Pump Handbook page 2.68

r_tongue = r_exit_impeller*(1.03+0.0025*specific_speed+7E-8*density_pump*head_pump); % m - Gulich table 10.2, pg 729
t_tongue = 0.04*r_exit_impeller; % m - Gulich pg 493
A_throat = ( vdot_pump/(3.8*r_exit_impeller*v_tangential_exit) * ( 1 + sqrt(1+15.2*r_exit_impeller*v_tangential_exit*(r_tongue+t_tongue)/vdot_pump) ) )^2; % m2
% where does this come from ^ 
r_throat = r_exit_impeller + 0.5*sqrt(A_throat) + t_tongue; % m

theta = linspace(0, 2*pi, 1000)'; % rad
A_volute = ( vdot_pump.*theta./(7.6*pi*r_exit_impeller*v_tangential_exit) .* ( 1 + sqrt(1+30.4*pi*r_exit_impeller*r_tongue*v_tangential_exit./(vdot_pump*theta)) ) ).^2; % m2
A_volute(isnan(A_volute)) = 0;

h_volute = sqrt(0.25*t_tongue^2+A_volute(end)) - 0.5*t_tongue; % m - height of volute gap - set so the outlet is a square
r_volute = r_tongue + A_volute/h_volute; % m - radius of volute wall
OD_volute_outlet = 2*sqrt(A_volute(end)/pi); % m - diameter of volute outlet hole
r_volute_outlet = r_volute(end) - h_volute(end)/2;

v_throat = vdot_pump./A_throat; % m/s - fluid exit velocity into plumbing

volute_curve = r_volute.*[cos(clock*theta), sin(clock*theta)]; % [m, m]

%%

% "...considerable freedom to configure the cross section without risking major losses ... flat cross sections result in less intense secondary flow than circular cross
% sections, thereby generating fewer losses ... a ratio of width B to height H in the range of B/H = 2 to 3 may be considered the  optimum." - Gulich pg 492

