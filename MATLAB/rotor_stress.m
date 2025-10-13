% Parameters
R = 0.16;                  % Disk radius [m]
rho = 7800;               % Density [kg/m^3]
omega = 30000*2*pi/60;    % Angular velocity [rad/s]
E = 200e9;                % Young's modulus [Pa]
nu = 0.29;                 % Poisson's ratio

% Radial discretization
N = 1000;
r = linspace(1e-6, R, N)';   % Avoid r=0 to prevent division by zero
dr = r(2) - r(1);

% Thickness profile f(r)
f = @(r) 0.01 - 0.001 * r;    % Example: linearly increasing thickness
f_vals = f(r);

% Construct matrices for finite difference
% Second derivative (central)
D2 = diag(-2*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
D2 = D2 / dr^2;

% First derivative (central)
D1 = diag(zeros(N,1));
for i = 2:N-1
    D1(i,i-1) = -1;
    D1(i,i+1) = 1;
end
D1 = D1 / (2*dr);

% Construct ODE system: A*u = b
A = E * (D2 + D1 .* diag(1./r) - diag(1./r.^2) * (1 - nu));
b = -rho * omega^2 * r;

% Apply boundary conditions
% u(0) = 0 (symmetry) → set u(1) = 0
% σ_r(R) = 0 → du/dr + nu * u / r = 0 at outer edge
A(1,:) = 0;
A(1,1:3) = [-1 0 1]/(2*dr);  % du/dr = 0 via central diff
b(1) = 0;
A(end,:) = 0;
A(end,end-1:end) = [-1, 1]/dr + nu * [1, 1]./(2*r(end));
b(end) = 0;

% Solve for u
u = A \ b;

% Compute strains
eps_r = [diff(u)/dr; 0];
eps_theta = u ./ r;

% Compute stresses
sigma_r = E ./ (1 - nu^2) .* (eps_r + nu * eps_theta);
sigma_theta = E ./ (1 - nu^2) .* (eps_theta + nu * eps_r);

% Find peak stress
peak_stress = max([abs(sigma_r); abs(sigma_theta)]);

% Plot
figure;
plot(r, sigma_r/1e6, r, sigma_theta/1e6);
legend('\sigma_r [MPa]', '\sigma_\theta [MPa]');
xlabel('Radius [m]'); ylabel('Stress [MPa]');
title('Stress Distribution in Rotating Disk');
grid on;
