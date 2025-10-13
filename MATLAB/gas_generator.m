clear cea_gg

p_gg = min(p_out)/(1+gg_stiffness);

cea_gg = CEA('problem','rocket', ...
     'equilibrium', ...  % equilibrium flow, infinite chamber, 
     'o/f',OF_gg, ...          % Define OF ratio
     'case','CEAM-rocket', ...
     'p(psi)',p_gg/6894.76,'pi/p',p_gg/p_amb, ... % Define nozzle with pressure ratio
     'reactants', ... 
     'fuel','C2H5OH(L)','wt%',proof,'t(k)',T_amb,'H2O(L)','wt%',1-proof,'t(k)',T_amb, ...
     'oxid','O2(L)','wt%',100, ... % Liquid Oxygen for Oxidizer, at boiling point
     'output','short','massf','transport','mks', ...
     'end');

Ae_At_gg = cea_gg.output.eql.aeat(end); % Expansion / Area Ratio
Tt_gg = cea_gg.output.eql.temperature(1); % K - stagnation, gg chamber temperature
gamma_gg = cea_gg.output.eql.gamma(end); % unitless - ratio of specific heats of final exhaust

c_star_ideal_gg = cea_gg.output.eql.cstar(1); % m/s - Characteristic Velocity
c_star_gg = c_star_ideal_gg*c_star_eff_gg; % m/s

p_exit_gg = 1E5*cea_gg.output.eql.pressure(end); % Pa - exit pressure
c_tau_ideal_gg = sqrt( (2*gamma_gg^2/(gamma_gg-1) * (2/(gamma_gg+1))^((gamma_gg+1)/(gamma_gg-1)) * (1-(p_amb/p_gg)^((gamma_gg-1)/gamma_gg) ) ) ) + Ae_At*(p_exit-p_amb)/p_gg; % Ideal Thrust Coefficient
lambda_cone_gg = 0.5*(1+cos(diverge_angle_gg)); % Divergence correction factor for conical nozzle
c_tau_gg = c_tau_ideal_gg*lambda_cone_gg*c_tau_eff_gg; % Thrust Coefficient

v_exhaust_ideal_gg = cea_gg.output.eql.sonvel(end)*cea_gg.output.eql.mach(end); % m/s - ideal exhaust velocity
v_exhaust_gg = c_star_gg*c_tau_gg; % m/s - actual exhaust velocity

dens_gas_spouted = cea_gg.output.eql.density(3); % kg/m3

%% Injector Design
p_gg_ox_in = p_gg*(1+gg_stiffness);
p_gg_fuel_in = p_gg*(1+gg_stiffness);

% technically should use densities at pump outlet pressure
% but compressibility should be pretty negligible

Cd_gg_ox = 0.80; % guess
Cd_gg_fuel = 0.80; % guess

A_gg_ox = mdot_ox_gg / (Cd_gg_ox*sqrt(2*density_ox*(p_gg_ox_in-p_gg))); % m2
A_gg_fuel = mdot_fuel_gg / (Cd_gg_fuel*sqrt(2*density_fuel*(p_gg_fuel_in-p_gg))); % m2

A_gg_fuel_orifice = A_gg_fuel / orifice_number_fuel; % m2

d_gg_ox_orifice = sqrt(A_gg_ox*4/pi); % m
d_gg_fuel_orifice = sqrt(A_gg_fuel_orifice*4/pi); % m

v_gg_ox = mdot_ox_gg / (density_ox * A_gg_ox); % m/s
v_gg_fuel = mdot_fuel_gg / (density_fuel * A_gg_fuel); % m/s
