% %% Setup
% clear,clc; close all;
% dllFolder = 'C:\Program Files\MATLAB\R2022a\bin\win64'; % libmex.dll libMatlabDataArray.dll
% daceFolder = 'C:\Program Files (x86)\DACE\lib'; % dace.dll
% setenv('PATH', [getenv('PATH') ';' dllFolder ';' daceFolder]);
% 
% clear; clc; close all;
% 
% %% Constants and Units
% mu = 3.986005e14; % Earth gravitational parameter [m^3/s^2]
% Re = 6378.137e3;  % Earth radius [m]
% J2 =  1.08262668e-3;
% g0 = 9.80665;     % [m/s^2]
% m_spacecraft = 800;
% 
% % Normalised units
% param.LU = Re;
% param.VU = sqrt(mu/Re);
% param.TU = param.LU/param.VU;
% param.MU = m_spacecraft;
% 
% % param.LU = 1;
% % param.VU = 1;
% % param.TU = 1;
% % param.MU = 1;
% 
% param.mu = mu / param.LU^3 * param.TU^2;
% param.Re = Re / param.LU;
% param.J2 = J2;
% 
% %% Define Chief initial state
% % Format: [a, e, i, RAAN, omega, nu]
% kep_chief  = [500e3 / param.LU + param.Re, 0.0, deg2rad(0), 0, 0, 0];
% 
% % Convert to Cartesian
% [rrC, vvC] = CoordConv.po2pv(kep_chief, param.mu);
% 
% x0C = [rrC; vvC] % Chief initial state
% 
% %% Define Deputy initial state
% % 
% % Relative motion parameters
% rho = 20 / param.LU % 20 m in normalised units
% omega = sqrt(param.mu / (kep_chief(1))^3); % mean motion
% v_rel = omega * rho; % tangential speed for circular motion
% % 
% % % Define relative position and velocity in Hill frame
% % r_rel_Hill = [rho; 0; 0];       % 20 m along radial
% % v_rel_Hill = [0; -1.5*v_rel; 0];     % tangential velocity
% % 
% % % LVLH frame rotation matrix at chief position
% % rC = rrC;
% % vC = vvC;
% % z_hat = -rC / norm(rC);
% % y_hat = cross(rC, vC); y_hat = y_hat / norm(y_hat);
% % x_hat = cross(y_hat, z_hat);
% % R_lvlh2eci = [x_hat, y_hat, z_hat];
% % 
% % % Convert deputy rel pos/vel to inertial frame
% % rD = rC + R_lvlh2eci * r_rel_Hill;
% % vD = vC + R_lvlh2eci * v_rel_Hill;
% % 
% % % Final deputy state
% % x0D = [rD; vD]
% % x0Ckep = CoordConv.vec2orbElem(x0C(1:3),x0C(4:6),param.mu)
% % x0Dkep = CoordConv.vec2orbElem(x0D(1:3),x0D(4:6),param.mu)
% dr = [rho; 0; 0];
% k = 1.00;
% dv = [0; -k*v_rel; k*v_rel];
% x0D = [rrC + dr; vvC + dv]
% 
% digits(50)
% x = vpa(dv(2));
% disp(x)
% 
% %% Period
% % Orbit period
% a = kep_chief(1);
% T = 2*pi*sqrt(a^3 / param.mu);  % normalised units
% tspan = linspace(0, 2*T, 2000); % 10 orbits
% 
% % Dynamics: two-body (or J2 if you want)
% odeOpts = odeset('RelTol',1e-10,'AbsTol',1e-12);
% twoBody = @(t,x) [x(4:6); -param.mu * x(1:3)/norm(x(1:3))^3];
% % twoBody = @(t, x) [ ...
% %     x(4:6); 
% %     -param.mu * x(1:3)/norm(x(1:3))^3 + ...
% %     j2Acceleration(x(1:3), param)
% % ];
% 
% [~, xC] = ode45(twoBody, tspan, x0C, odeOpts);
% [~, xD] = ode45(twoBody, tspan, x0D, odeOpts);
% 
% %% To LVLH
% r_rel_LVLH = zeros(3, length(tspan));
% for k = 1:length(tspan)
%     rC = xC(k,1:3).'; vC = xC(k,4:6).';
%     rD = xD(k,1:3).';
% 
%     % LVLH frame axes
%     z_hat = -rC / norm(rC);
%     y_hat = cross(rC, vC); y_hat = y_hat / norm(y_hat);
%     x_hat = cross(y_hat, z_hat);
%     R_eci2lvlh = [x_hat, y_hat, z_hat].';
% 
%     % Relative position in ECI
%     r_rel = rD - rC;
% 
%     % Express in LVLH
%     r_rel_LVLH(:,k) = R_eci2lvlh * r_rel;
% end
% 
% % %% Plot
% % figure; hold on;
% % 
% % % Main trajectory
% % p = plot(r_rel_LVLH(1,:), r_rel_LVLH(2,:), 'b');
% % 
% % % Start and end markers
% % plot(r_rel_LVLH(1,1), r_rel_LVLH(2,1), 'o', 'MarkerSize', 8, ...
% %      'MarkerFaceColor', [1, 0.5, 0], 'MarkerEdgeColor', 'k'); % orange start
% % 
% % plot(r_rel_LVLH(1,end), r_rel_LVLH(2,end), 'o', 'MarkerSize', 8, ...
% %      'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k'); % blue end
% % 
% % % Format using plot_latex
% % plot_latex(p, ...
% %     'Radial [Re]', ...
% %     'Along-track [Re]', ...
% %     '', ...
% %     '', ...
% %     {'', 'Epoch', 'Final'});
% % axis equal;
% 
% %% 3D Relative Motion Plot in LVLH
% figure;
% plot3(-r_rel_LVLH(1,:), -r_rel_LVLH(2,:), -r_rel_LVLH(3,:), 'b', 'LineWidth', 1.5); hold on;
% % plot(-r_rel_LVLH(1,:), -r_rel_LVLH(3,:), 'b', 'LineWidth', 1.5); hold on;
% 
% % Start marker
% plot3(-r_rel_LVLH(1,1), -r_rel_LVLH(2,1), -r_rel_LVLH(3,1), 'o', ...
%     'MarkerSize', 8, 'MarkerFaceColor', [1, 0.5, 0], 'MarkerEdgeColor', 'k'); % orange
% 
% % End marker
% plot3(-r_rel_LVLH(1,end), -r_rel_LVLH(2,end), -r_rel_LVLH(3,end), 'o', ...
%     'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k'); % blue
% 
% plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% 
% grid on;
% % axis equal;
% xlabel('Along-track [Re]', 'Interpreter', 'latex');
% ylabel('Cross-track [Re]', 'Interpreter', 'latex');
% zlabel('Radial [Re]', 'Interpreter', 'latex');
% title('3D Relative Motion in LVLH Frame', 'Interpreter', 'latex');
% legend({'Trajectory', 'Start', 'End'}, 'Interpreter', 'latex');
% 
% % Helper function to compute J2 perturbation acceleration
% function a_J2 = j2Acceleration(r, param)
%     mu = param.mu;
%     J2 = param.J2;
%     Re = param.Re;
% 
%     x = r(1); y = r(2); z = r(3);
%     r_norm = norm(r);
%     
%     factor = 1.5 * J2 * mu * Re^2 / r_norm^5;
%     z2_r2 = 5 * (z^2) / r_norm^2;
% 
%     a_J2 = factor * [
%         x * (1 - z2_r2);
%         y * (1 - z2_r2);
%         z * (3 - z2_r2 * 5)
%     ];
% end


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
