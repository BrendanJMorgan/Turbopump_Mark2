% NASA SP-8052 "LIQUID ROCKET ENGINE TURBOPUMP INDUCER" 1971
%% Net Positive Suction Head (NPSH) Requirements

k1 = 1.2; 
k2 = 0.2334 + (pump_shaft_speed(p)*r_eye_impeller(p)/128.3)^4;
cavitation_limit(p) = (pump_shaft_speed(p)*r_eye_impeller(p))^2/(2*g); % unitless
NPSH_SE(p) = ((k1+k2)*eye_flow_coeff(p)^2 + k2) * cavitation_limit(p); % m - NPSH, shockless entry. Pump handbook pg. 2.88
NPSH_required(p) = NPSH_SE(p); % m - assuming operation at best efficiency point (BEP) and shockless

NPSH_available(p) = (p_amb - vapor_pressure_pump(p))/(density_pump(p)*g); % m - Net positive suction head. Assuming zero plumbing water column and inlet velocity
NPSH_inducer(p) = (NPSH_required(p) - NPSH_available(p))*NPSH_margin(p); % m

% Should add some thought into NPSH upstream of inducer. Should cavitation
% on leading edge be prevented?
%% Tip Diameter
vdot_corrected(p) = vdot_pump(p) / (1-hub_tip_ratio_inducer(p)^2); % m3/s - volumetric flow rate through open portion of inducer

ss_speed(p) = pump_shaft_speed(p)*sqrt(vdot_pump(p))/(g * NPSH_inducer(p))^0.75; % unitless - suction specific speed
flow_coeff_inducer(p) = 1.3077*sqrt(1-hub_tip_ratio_inducer(p)^2)/ss_speed(p) / (1 + 0.5*sqrt( 1+6*(1.3077*sqrt(1-hub_tip_ratio_inducer(p)^2)/ss_speed(p))^2)); % unitless - phi
r_tip_inducer(p) = 0.5*0.3048*0.37843*(vdot_corrected(p)*15850.323 / (pump_shaft_speed(p)*30/pi*flow_coeff_inducer(p)))^(1/3); % m
r_hub_inducer(p) = r_tip_inducer(p)*hub_tip_ratio_inducer(p); % m

r_pitchline_inducer = 0.5*(1+hub_tip_ratio_inducer(p))*r_tip_inducer(p); % m

%% Blade Construction
v_inlet_inducer(:,p) = [0; 0];
v_inlet_inducer(1,p) = pump_shaft_speed(p)*r_tip_inducer(p); % m/s - tangential velocity at inducer pitchline tip
v_inlet_inducer(2,p) = vdot_corrected(p) / (pi*r_tip_inducer(p)^2); % m/s - axial velocity at inducer pitchline

flow_blade_ratio(p) = 0.35; % unitless

flow_incidence_angle(p) = atan(v_inlet_inducer(2) / v_inlet_inducer(1)); % rad - alpha
blade_angle_inducer(p) = flow_incidence_angle(p)/(1-flow_blade_ratio(p)); % rad - beta

blade_lead_inducer(p) = 2*pi*r_tip_inducer(p)*tan(blade_angle_inducer(p)); % m - blade lead - how much the blade advances per turn
solidity_inducer(p) = 2.5;
h_min_inducer(p) = blade_lead_inducer(p)*solidity_inducer(p)/blade_number_inducer(p)*sin(blade_angle_inducer(p)); % m - height of inducer

%% Axial Clearance
% Between inducer and impeller
clearance_axial_inducer(p) = blade_lead_inducer(p) / blade_number_inducer(p); % m

%% Head
v_tip_inducer(p) = pump_shaft_speed(p)*r_tip_inducer(p); % m/s
head_coeff_inducer(p) = NPSH_inducer(p)*g / v_tip_inducer(p)^2; % unitless - 0.15 is the cutoff between low and high head - camber needed if >0.075

%% Hub Construction
r_hub_base(p) = r_eye_inner(p); % m
r_hub_top(p) = 2*hub_tip_ratio_inducer(p)*r_tip_inducer(p) - r_hub_base(p); % m

%% Clearance Losses - not relevant with a proper shroud
r_cavity_inducer(p) = r_tip_inducer(p) + clearance_radial_inducer(p); % m
ss_speed_clearanced(p) = ss_speed(p)*(1-0.575*sqrt((r_cavity_inducer(p)-r_tip_inducer(p))/(r_tip_inducer(p)-r_hub_inducer(p)))); % unitless
head_coeff_inducer_clearanced(p) = head_coeff_inducer(p)*(1-1.0*sqrt((r_cavity_inducer(p)-r_tip_inducer(p))/(r_tip_inducer(p)-r_hub_inducer(p)))); % unitless

%% Other
inducer_torque(p) = mdot_pump(p)*pump_shaft_speed(p)*r_tip_inducer(p)^2; % N*m - this is a conservative maximum
