%% Setup
clear;clc; close all;
dllFolder = 'C:\Program Files\MATLAB\R2022a\bin\win64';
daceFolder = 'C:\Program Files (x86)\DACE\lib';
setenv('PATH', [getenv('PATH') ';' dllFolder ';' daceFolder]);

clear; clc; close all;

%% Constants and Units
mu = 3.986005e14; % [m^3/s^2]
Re = 6378.137e3;  % [m]
J2 =  1.08262668e-3;
m_spacecraft = 800;

param.mu = mu;
param.Re = Re;
param.J2 = J2;

%% Chief initial orbit [a, e, i, RAAN, omega, nu]
kep_chief = [500e3 + Re, 0.0, deg2rad(0), 0, 0, 0];
kep_dep = [6871.2657e3, 0.001, deg2rad(0), 0, deg2rad(180.09), deg2rad(180)];

% Convert to Cartesian
[rrC, vvC] = CoordConv.po2pv(kep_chief, param.mu);
[rrD, vvD] = CoordConv.po2pv(kep_dep, param.mu);

x0C = [rrC; vvC];
x0D = [rrD, vvD];

%% Deputy Initial Conditions
rho = 20; % m
omega = sqrt(param.mu / kep_chief(1)^3);
v_rel = omega * rho;

% dr = [rho; 0; 0];
% dv = [0; -v_rel; v_rel]; %% Parking orbit
% dv = [0; -1.2*v_rel; v_rel]; %% Final approach
% x0D = [rrC + dr; vvC + dv];

% digits(50); disp(vpa(dv(2)));

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
v_rel_LVLH = zeros(3, length(tspan));
for k = 1:length(tspan)
    rC = xC(k,1:3).'; vC = xC(k,4:6).';
    rD = xD(k,1:3).'; vD = xD(k,4:6).';

    % LVLH frame axes (ECI basis)
    z_hat = -rC / norm(rC);
    y_hat = cross(rC, vC); y_hat = y_hat / norm(y_hat);
    x_hat = cross(y_hat, z_hat);
    R_eci2lvlh = [x_hat, y_hat, z_hat].';

    % Relative position
    r_rel = rD - rC;
    r_rel_LVLH(:,k) = R_eci2lvlh * r_rel;

    % Angular velocity of LVLH frame (in ECI)
    hC = cross(rC, vC);
    omega_lvlh = hC / norm(rC)^2;

    % Relative velocity in inertial frame
    v_rel = vD - vC;

    % Transform to LVLH frame
    v_rel_LVLH(:,k) = R_eci2lvlh * (v_rel - cross(omega_lvlh, r_rel));
end

%% 3D Plot of Relative Motion
figure;

% Plot trajectory without legend entry
traj = plot3(r_rel_LVLH(1,:), -r_rel_LVLH(2,:), -r_rel_LVLH(3,:), 'b', ...
    'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on;

% Plot start marker for deputy
start_marker = plot3(r_rel_LVLH(1,end), -r_rel_LVLH(2,end), -r_rel_LVLH(3,end), 'bo', ...
    'MarkerSize', 8, 'MarkerFaceColor', [0, 0, 1]);

r_rel_LVLH(1,end), -r_rel_LVLH(2,end), -r_rel_LVLH(3,end)
v_rel_LVLH(1,end), -v_rel_LVLH(2,end), -v_rel_LVLH(3,end)

% Plot chief marker at origin
chief_marker = plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

projectionColor = [0.678, 0.847, 0.902];
lineWidth = 2;

% Project to x–y plane at z = 20
plot3(r_rel_LVLH(1,:), -r_rel_LVLH(2,:), -14600 * ones(1, size(r_rel_LVLH,2)), ...
    'Color', projectionColor, 'LineWidth', lineWidth);

% Project to x–z plane at y = 20
plot3(r_rel_LVLH(1,:), 1 * ones(1, size(r_rel_LVLH,2)), -r_rel_LVLH(3,:), ...
    'Color', projectionColor, 'LineWidth', lineWidth);

% Project to y–z plane at x = 50
plot3(-147200 * ones(1, size(r_rel_LVLH,2)), -r_rel_LVLH(2,:), -r_rel_LVLH(3,:), ...
    'Color', projectionColor, 'LineWidth', lineWidth);

% Final approach
% end_x = -r_rel_LVLH(1, end)
% 
% % Plot trajectory without legend entry
% traj = plot3(-r_rel_LVLH(1,:) - end_x, -r_rel_LVLH(2,:), -r_rel_LVLH(3,:), 'b', ...
%     'LineWidth', 1.5, 'HandleVisibility', 'off'); hold on;
% 
% % Plot start marker for deputy
% start_marker = plot3(-r_rel_LVLH(1,1)-end_x, -r_rel_LVLH(2,1), -r_rel_LVLH(3,1), 'bo', ...
%     'MarkerSize', 8, 'MarkerFaceColor', [0, 0, 1]);
% 
% % Plot chief marker at origin
% chief_marker = plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% 
% projectionColor = [0.678, 0.847, 0.902];
% lineWidth = 2;

% % Project to x–y plane at z = 20
% plot3(-r_rel_LVLH(1,:)-end_x, -r_rel_LVLH(2,:), -40 * ones(1, size(r_rel_LVLH,2)), ...
%     'Color', projectionColor, 'LineWidth', lineWidth);
% 
% % Project to x–z plane at y = 20
% plot3(-r_rel_LVLH(1,:)-end_x, 21 * ones(1, size(r_rel_LVLH,2)), -r_rel_LVLH(3,:), ...
%     'Color', projectionColor, 'LineWidth', lineWidth);
% 
% % Project to y–z plane at x = 50
% plot3(-190 * ones(1, size(r_rel_LVLH,2)), -r_rel_LVLH(2,:), -r_rel_LVLH(3,:), ...
%     'Color', projectionColor, 'LineWidth', lineWidth);
% % 416

% View angle
view(52.5, 30);
lgd = legend('Location', 'northeast'); 

% Use plot_latex for formatting
% axis equal;
plot_latex(traj, ...
    'Along-track [m]', ...
    'Cross-track [m]', ...
    'Radial [m]', ...
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
