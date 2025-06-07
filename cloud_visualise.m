%% Setup
clear,clc; close all;
dllFolder = 'C:\Program Files\MATLAB\R2022a\bin\win64'; % libmex.dll libMatlabDataArray.dll
daceFolder = 'C:\Program Files (x86)\DACE\lib'; % dace.dll
setenv('PATH', [getenv('PATH') ';' dllFolder ';' daceFolder]);

clear; clc; close all;

%% Constants and Units
mu = 3.986005e14; % Earth gravitational parameter [m^3/s^2]
Re = 6378.137e3;  % Earth radius [m]
g0 = 9.80665;     % [m/s^2]
m_spacecraft = 800;

% Normalised units
param.LU = Re;
param.VU = sqrt(mu/Re);
param.TU = param.LU/param.VU;
param.MU = m_spacecraft;
param.mu = mu / param.LU^3 * param.TU^2;
param.Re = Re / param.LU;

%% Define Chief and Deputy (Keplerian Elements)
% Format: [a, e, i, RAAN, omega, nu]
kep_chief  = [500e3 + Re, 0, deg2rad(99), 0, 0, 0];
kep_deputy = [501e3 + Re, 0, deg2rad(98), 0, 0, deg2rad(.1)];

% Convert to Cartesian
[rrC, vvC] = CoordConv.po2pv(kep_chief, mu);
[rrD, vvD] = CoordConv.po2pv(kep_deputy, mu);

X0C = [rrC; vvC]; % Chief initial state
X0D = [rrD; vvD]; % Deputy initial state

%% Propagation Time
a0 = kep_chief(1);
Torb = 2*pi*sqrt(a0^3/mu); % Orbital period
tspan = linspace(0, 2*Torb, 1000); % Two orbits

% %% Propagate Chief and Deputy
opts = odeset('RelTol',1e-9,'AbsTol',1e-12);
% [~, Xchief]  = ode45(@(t,x) twoBodyODE(t,x,mu), tspan, X0C, opts);
% [~, Xdeputy] = ode45(@(t,x) twoBodyODE(t,x,mu), tspan, X0D, opts);

%% Setup and Perturbation Parameters
nPoints = 50; % number of perturbed states
r_var = 1e3 / param.LU; 
v_var = 1 / param.LU * param.TU;

% Save initial positions for plotting
perturbed_initials = zeros(nPoints, 3);
perturbed_finals = zeros(nPoints, 3);

for i = 1:nPoints
    val(i,:) = [randn(1,3), 0,0,0]'.*[r_var*ones(3,1); zeros(3,1)] + ...
        [0,0,0,randn(1,3)]'.*[zeros(3,1); v_var*ones(3,1)];

    % Update the cartesian initial state to the perturbed initial state
    CartXnew = X0D'+val(i,:); 

    % Save initial position
    perturbed_initials(i,:) = CartXnew(1:3) * param.LU;

    % Propagate perturbed deputy
    [~, Xdeputy_perturbed] = ode45(@(t,x) twoBodyODE(t,x,mu), tspan, CartXnew, opts);

    % Save final propagated state
    perturbed_finals(i,:) = Xdeputy_perturbed(end,1:3) * param.LU;
end

%% Plot Initial State Cloud
figure;
scatter3(perturbed_initials(:,1)/1e3, perturbed_initials(:,2)/1e3, perturbed_initials(:,3)/1e3, 'r.');
hold on;
% plot3(X0C(1)*param.LU/1e3, X0C(2)*param.LU/1e3, X0C(3)*param.LU/1e3, 'ko', 'MarkerFaceColor', 'k');
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Initial State');
% legend('Deputy Perturbed', 'Chief');
axis equal; grid on;

%% Plot Propagated State Cloud
figure;
scatter3(perturbed_finals(:,1)/1e3, perturbed_finals(:,2)/1e3, perturbed_finals(:,3)/1e3, 'r.');
hold on;
plot3(Xchief(end,1)*param.LU/1e3, Xchief(end,2)*param.LU/1e3, Xchief(end,3)*param.LU/1e3, 'ko', 'MarkerFaceColor', 'k');
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Propagated State After 2 Orbits');
legend('Deputy Propagated', 'Chief Propagated');
axis equal; grid on;

%% Helper Function
function dx = twoBodyODE(~, x, mu)
    r = x(1:3);
    v = x(4:6);
    a = -mu * r / norm(r)^3;
    dx = [v; a];
end

%% Compute Relative State in LVLH Frame
% rel_LVLH = zeros(length(tspan), 6);
% for k = 1:length(tspan)
%     rC = Xchief(k,1:3)';
%     vC = Xchief(k,4:6)';
%     rD = Xdeputy(k,1:3)';
%     vD = Xdeputy(k,4:6)';
%     
%     R = rD - rC;
%     V = vD - vC;
% 
%     % Compute LVLH basis
%     z = -rC/norm(rC);
%     y = -cross(rC,vC)/norm(cross(rC,vC));
%     x = cross(y,z);
% 
%     Q = [x y z]; % DCM from ECI to LVLH
%     omega = cross(rC,vC)/norm(rC)^2; % Angular velocity of LVLH frame
% 
%     rel_pos_LVLH = Q' * R;
%     rel_vel_LVLH = Q' * (V - cross(omega, R));
% 
%     rel_LVLH(k,:) = [rel_pos_LVLH' rel_vel_LVLH'];
% end

