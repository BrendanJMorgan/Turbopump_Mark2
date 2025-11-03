clear all
format short g
PropsSI = @py.CoolProp.CoolProp.PropsSI;
rocketprop = @py.rocketprops.rocket_prop.get_prop;

%% Inputs 

% Constants
g = 9.81;				% m/s2
p_amb = 93900;          % Pa - ambient pressure at 2100 feet elevation
T_amb = 293;			% K - Ambient Temperature - lowest temperature for which 75/25 ethanol is solvable by CoolProp

% Overall Engine Performance Targets
thrust_target = 15000;  % N - Thrust
p_cc = 30E5;			% Pa - Chamber (Stagnation) Pressure
cc_stiffness = 0.25;	% Pa/Pa - guess   

% Combustion Chamber (CC)
fuel = 'Jet-A(L)';        % 'Jet-A(L)' 'C2H5OH(L)'
proof = 1;                % Proportion of fuel (remanining part is water)
OF = 2.0;                 % Oxidizer/Fuel Ratio (by mass)
c_star_eff = 0.75;        % Characteristic Vel Efficiency, experimental
c_tau_eff = 0.96;         % Thrust Coefficient Efficiency Factor
gamma_guess = 1.15;
c_tau_guess = 0.983*c_tau_eff*sqrt( (2*gamma_guess^2/(gamma_guess-1) * (2/(gamma_guess+1))^((gamma_guess+1)/(gamma_guess-1)) * (1-(p_amb/p_cc)^((gamma_guess-1)/gamma_guess) ) ) ); 
A_throat = thrust_target / (p_cc*c_tau_guess*c_star_eff); % m2 - Throat Area

% Gas Generator (GG)
gg_stiffness = 0.25;		% Pa/Pa - guess
gg_fraction_guess = 0.05;   % unitless - Fraction of total mass flow sent to the gas generator. Context: F1 = 0.030, J2 = 0.014 DO NOT CHANGE FROM 0.01 NOW
OF_gg = 0.4;				% unitless - OF Ratio - "[Most] operate at mixture ratios from 0.2 to 1.0, with hydrocarbons falling in the lower end, about 0.3" (NASA 1972)
c_star_eff_gg = 0.75;       % m/s - characteristic velocity efficiency, experimental ANY PAPERS ON THIS?
c_tau_eff_gg = 0.96;        % unitless - Thrust Coefficient Efficiency Factor
orifice_number_fuel = 12;   % 4 for impinging, 8 for showerhead film cooling on gg injector

% Pumps (first value is for ox pump, second for fuel pump)
clocking = [1,-1];                      % 1 for counterclockwise and -1 for clockwise (looking down at pump inlet)
pump_shaft_speed = 30000*pi/30*[1,1];   % rad/s - angular velocity of the pump shaft, impeller, and inducers
gear_efficiency = 1;                    % unitless - spur gears ~= 0.95; common shaft = 1
r_hub = [0.006,0.006];	        % m
r_shaft = 0.005;                        % m - portion of shaft that is stainless steel
impeller_thickness = [0.003,0.003];     % m - thickness of impeller at the exit point, not including blades
eye_flow_coeff = [0.25,0.25];           % unitless - higher means smaller impeller but larger inducer
                                        % phi_e in pump handbook; pg 2.29: 0.2-0.3 for impellers, ~0.1 or less for inducers
NPSH_margin = [1.5,1.5];                % unitless - Margin of extra net positive suction head (NPSH) to be provided by inducers to prevent cavitation
blade_number_inducer = [4,4];           % unitless - 3 or 4 is considered good
clearance_radial_inducer = [2E-4,2E-4]; % m - radial clearance between inducer blades and housing cavity
surface_roughness = 10E-6;              % m - surface roughness of the additive material of the impellers

% Turbine
r_pitchline_rotor = 0.070;      % m
gear_ratio = 1;                 % unitless - higher makes for a faster, smaller turbine. 1 = common shaft
d_throat_nozzle = 0.127*0.0254; % m - diameter of each nozzle leading off the manifold - mdot_gg is a direct function of this and nozzle_number
diverge_angle_gg = 15*pi/180;   % rad - half-cone divergence angle of nozzle plate
blade_width_rotor = 0.5*0.0254; % m - real turbines seems to be about 0.5 inches or so
admission_fraction = 20/360;    % unitless - fraction of the nozzle plate circle that actually has nozzles
blade_length_rotor = 1/4*0.0254;% m
tip_clearance = 0.0005;         % m
radius_leading = 0.002;         % m - blade leading edge fillet radius. NEEDS MORE RESEARCH

% Chamber/Nozzle Geometry
dx = 0.001;                 % m - position step 
converge_angle = 45*pi/180; % rad - half-cone convergence angle of combustion chamber end
diverge_angle = 15*pi/180;  % rad - half-cone divergence angle of nozzle
l_star = 1;                 % m - combustion length UPDATE FOR JET-A
rc_throat = 0.025;          % m - radius of curvature around the throat
d2_chamber = 0.125;         % m
thickness = 0.006;          % m % SLANT VS VERTICAL THICKNESS
d1_chamber = d2_chamber - 2*thickness; % m
r1_chamber = d1_chamber/2;  % m

% Coolant Channels
n_pipe1 = 16;   % number of channels along barrel
n_pipe2 = 16;   % number of channels near throat
n_pipe3 = 16;   % number of channels along lower nozzle section
gap_pipe = 1/4*0.0254;          % m - Gap between channels (fin thickness)
h_pipe = 1/16*0.0254;           % m - coolant channel height
merge_radius = 0.45*d1_chamber; % m - when contour is below this radius, transition to n_pipe2
flow_direction = -1;            % 1 = forward flow (injector to nozzle), -1 = counter flow (nozzle to injector)

% Film Cooling
film_fraction = 0.03;           % unitless - Fraction of the fuel mass flow dedicated to film cooling orifices - typically 3%-10% (Huzel and Huang)
v_injection = 10;               % m/s - combustion gas must have some initial velocity for injector film cooling to work mathematically
injection_efficiency = 1.0;     % I forgot what this is

%% Properties

k_wall = 253; % W/m-K - thermal conductivity of aluminum (at 571 C)

FSu_design = 2.0; % Minimum design factor of safety, ultimate

shaft_tensile_strength = 510E6; % Pa - stainless steel 304, ultimate
shaft_shear_strength = shaft_tensile_strength/sqrt(3);  % Pa - conservative Tresca criterion


%% Runs - these both have several (nested) subfunctions

thermals = false;
thrust_chamber_assembly
powerhead

%% Results

thrust
isp_ideal
isp_real
mdot_gg
mdot_total

plots

