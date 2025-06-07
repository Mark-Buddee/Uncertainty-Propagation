%% Setup
clear;clc; close all;
dllFolder = 'C:\Program Files\MATLAB\R2022a\bin\win64';
daceFolder = 'C:\Program Files (x86)\DACE\lib';
setenv('PATH', [getenv('PATH') ';' dllFolder ';' daceFolder]);

clear; clc; close all;

%% Constants and Units
mu = 3.986005e14; % Earth gravitational parameter [m^3/s^2]
Re = 6378.137e3;  % Earth radius [m]
J2 =  1.08262668e-3;
% J2 = 0;
g0 = 9.80665;     % [m/s^2]
m_spacecraft = 800;

% Normalised units
param.LU = Re;
param.VU = sqrt(mu/Re);
param.TU = param.LU/param.VU;
param.MU = m_spacecraft;

param.mu = mu / param.LU^3 * param.TU^2;
param.Re = Re / param.LU;
param.J2 = J2;

%% Chief initial orbit [a, e, i, RAAN, omega, nu]
kep_chief  = [500e3 / param.LU + param.Re, 0, deg2rad(97.4), 0, 0, 0];

[rrC, vvC] = CoordConv.po2pv(kep_chief, param.mu);
x0C = [rrC; vvC]

%% Deputy Initial Conditions
rho = 20 / param.LU; % m
omega = sqrt(param.mu / kep_chief(1)^3);
v_rel = omega * rho;

dr = [rho; 0; 0];
dv = [0; -v_rel; v_rel]; %% Parking orbit
% dv = [0; -1.2*v_rel; v_rel]; %% Final approach
% x0D = [rrC + dr; vvC + dv];

% digits(50); disp(vpa(dv(2)));
% digits(50); disp(vpa(x0D(5) - x0C(5)))

%% 

% Compute RTN frame from chief state
rhat = rrC / norm(rrC);
hhat = cross(rrC, vvC) / norm(cross(rrC, vvC));
RTN = [rhat'; cross(hhat, rhat)'; hhat'];

% Convert perturbation to ECI
rho_RTN = [dr; dv];
perturbation_ECI = [RTN' zeros(3); zeros(3) RTN'] * rho_RTN;

% Apply to chief to get deputy
x0D = x0C + perturbation_ECI

kep_chief_unnormalised = CoordConv.vec2orbElem(x0C(1:3)*param.LU - [Re; 0; 0], x0C(4:6)*param.VU, mu)
kep_deputy_unnormalised = CoordConv.vec2orbElem(x0D(1:3)*param.LU, x0D(4:6)*param.VU, mu)

%% Orbital Period
a = kep_chief(1);
T = 2*pi*sqrt(a^3 / param.mu);
tspan = linspace(0,2*T, 2000);

odeOpts = odeset('RelTol',1e-10,'AbsTol',1e-12);
twoBody = @(t,x) [x(4:6); -param.mu * x(1:3)/norm(x(1:3))^3];

[~, xC] = ode45(twoBody, tspan, x0C, odeOpts);
[~, xD] = ode45(twoBody, tspan, x0D, odeOpts);

%% Relative Position in LVLH Frame
r_rel_LVLH = zeros(3, length(tspan));
for k = 1:length(tspan)
    rC = xC(k,1:3).'; vC = xC(k,4:6).';
    rD = xD(k,1:3).';

    z_hat = -rC / norm(rC);
    y_hat = cross(rC, vC); y_hat = y_hat / norm(y_hat);
    x_hat = cross(y_hat, z_hat);
    R_eci2lvlh = [x_hat, y_hat, z_hat].';

    r_rel = rD - rC;
    r_rel_LVLH(:,k) = R_eci2lvlh * r_rel;
end

%% 3D Plot of Relative Motion
figure;

% Plot trajectory without legend entry
traj = plot3(r_rel_LVLH(1,:), -r_rel_LVLH(2,:), -r_rel_LVLH(3,:), 'b', ...
    'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on;

% Plot start marker for deputy
start_marker = plot3(r_rel_LVLH(1,end), -r_rel_LVLH(2,end), -r_rel_LVLH(3,end), 'bo', ...
    'MarkerSize', 8, 'MarkerFaceColor', [0, 0, 1]);

% Plot chief marker at origin
chief_marker = plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

projectionColor = [0.678, 0.847, 0.902];
lineWidth = 2;

% % Project to x–y plane at z = 20
% plot3(r_rel_LVLH(1,:), -r_rel_LVLH(2,:), -14600 * ones(1, size(r_rel_LVLH,2)), ...
%     'Color', projectionColor, 'LineWidth', lineWidth);
% 
% % Project to x–z plane at y = 20
% plot3(r_rel_LVLH(1,:), 1 * ones(1, size(r_rel_LVLH,2)), -r_rel_LVLH(3,:), ...
%     'Color', projectionColor, 'LineWidth', lineWidth);
% 
% % Project to y–z plane at x = 50
% plot3(-147200 * ones(1, size(r_rel_LVLH,2)), -r_rel_LVLH(2,:), -r_rel_LVLH(3,:), ...
%     'Color', projectionColor, 'LineWidth', lineWidth);
% 
% % Project to x–y plane at z = 20
% plot3(r_rel_LVLH(1,:), -r_rel_LVLH(2,:), -20 * ones(1, size(r_rel_LVLH,2)), ...
%     'Color', projectionColor, 'LineWidth', lineWidth);
% 
% % Project to x–z plane at y = 20
% plot3(-r_rel_LVLH(1,:), 21 * ones(1, size(r_rel_LVLH,2)), -r_rel_LVLH(3,:), ...
%     'Color', projectionColor, 'LineWidth', lineWidth);
% 
% % Project to y–z plane at x = 50
% plot3(-50 * ones(1, size(r_rel_LVLH,2)), -r_rel_LVLH(2,:), -r_rel_LVLH(3,:), ...
%     'Color', projectionColor, 'LineWidth', lineWidth);

% View angle
view(52.5, 30);
lgd = legend('Location', 'northeast'); 

% Use plot_latex for formatting
% axis equal;
plot_latex(traj, ...
    'Along-track [Re]', ...
    'Cross-track [Re]', ...
    'Radial [Re]', ...
    '', ...
    {'Deputy Satellite', 'Chief Satellite'});

%% Optional J2 (unused here)
function a_J2 = j2Acceleration(r, param)
    mu = param.mu;
    J2 = param.J2;
    Re = param.Re;
    x = r(1); y = r(2); z = r(3);
    r_norm = norm(r);
    
    factor = 1.5 * J2 * mu * Re^2 / r_norm^5;
    z2_r2 = 5 * (z^2) / r_norm^2;

    a_J2 = factor * [
        x * (1 - z2_r2);
        y * (1 - z2_r2);
        z * (3 - z2_r2 * 5)
    ];
end
