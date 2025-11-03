%% Establishing the Pump Configuration
dr = 0.001;
head_pump(p) = (p_out(p) - p_in(p)) / (density_pump(p)*g); % m
vdot_pump(p) = mdot_pump(p)/density_pump(p); % m3/s - Volumetric Flow Rate

%% Sizing the Pump

specific_speed(p) = abs(pump_shaft_speed(p))*sqrt(vdot_pump(p))/(head_pump(p)*g)^0.75;
specific_speed_imperial(p) = specific_speed(p)*2729.64; % pseudo unitless - Ns
outlet_flow_coeff(p) = 0.1715*sqrt(specific_speed(p)); % phi_i in pump handbook

ft = 1.0; % could increase up to 1.1 to be more aggressive on head, on low head pumps
head_coeff(p) = 0.605*ft*exp(-0.408*specific_speed(p)); % unitless - psi_opt; eqn 3.26 in Gulich
shutoff_coeff(p) = 1.25; % 1.25 for volute designs and 1.31 for diffuser designs
shutoff_head_coeff(p) = 0.625*exp(-0.159*specific_speed(p)); % unitless - psi_0; in the limit Q=0

blockage(p) = 0.85; % unitless - iterative initial guess - 1 is completely open, pump handbook says 0.85 is typical?
old_blockage = 0; % placeholder
blade_count(p) = 16; % unitless - iterative initial guess
old_blade_count = 0; % placeholder
slip_factor(:,p) = 0.15*ones(1000,1); % unitless - iterative initial guess


if specific_speed(p) < 0.1 || specific_speed(p) > 0.6
    fprintf("Pump Specific Speed is %g. Recommended is 0.1 to 0.6 for centrifugal geometries.\n", specific_speed)
end

%% Subfunctions
while (abs(blockage(p)-old_blockage) > 0.01)
    old_blockage = blockage(p)
    old_blade_count = blade_count(p)
    impeller % uses blockage parameter
    blades % computes a new blockage parameter
end
volute
inducer

%% Requirements for Turbine
pump_efficiency(p) = 0.94 - 6.983E-4*(vdot_pump(p)/(surface_roughness^2*abs(pump_shaft_speed(p))))^-0.2133 - 0.29*log10(0.83644/specific_speed(p))^2; % unitless - Anderson; pump handbook pg 2.25
power_pump(p) = vdot_pump(p) * head_pump(p) * density_pump(p) * g / pump_efficiency(p); % W
power_pump_theory(p) = vdot_pump(p) * head_pump(p) * density_pump(p) * g; % W

% Other
d_inlet_pump(p) = 3/8*0.0254; % m
v_inlet_pump(p) = vdot_pump(p) / (pi/4*d_inlet_pump(p)^2); % m/s
cavitation_inlet(p) = (p_amb - vapor_pressure_pump(p)) / (0.5*density_pump(p)*v_inlet_pump(p)^2); % unitless

d_outlet_pump(p) = 3/8*0.0254; % m
v_outlet_pump(p) = vdot_pump(p) / (pi/4*d_outlet_pump(p)^2); % m/s
cavitation_outlet(p) = (p_out(p) - vapor_pressure_pump(p)) / (0.5*density_pump(p)*v_outlet_pump(p)^2); % unitless

v_volute_pump(p) = sqrt(2*power_pump(p) / mdot_pump(p)); % m/s - average tangential velocity along the contour of the volute








