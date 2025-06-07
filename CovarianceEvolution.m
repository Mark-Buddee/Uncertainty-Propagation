%% Setup
clear,clc; close all;
dllFolder = 'C:\Program Files\MATLAB\R2022a\bin\win64'; % libmex.dll libMatlabDataArray.dll
daceFolder = 'C:\Program Files (x86)\DACE\lib'; % dace.dll
setenv('PATH', [getenv('PATH') ';' dllFolder ';' daceFolder]);

clear; clc; close all;

%% Constants and Units
mu = 3.986005e14; % Earth gravitational parameter [m^3/s^2]
Re = 6378.137e3;  % Earth radius [m]
% J2 =  1.08262668e-3;
J2 = 0;
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

%% Define Chief and Deputy (Keplerian Elements)
% Format: [a, e, i, RAAN, omega, nu]
kep_chief  = [500e3 / param.LU + param.Re, 0.01, deg2rad(0.1), 0, 0, 0.1];
kep_deputy = [501e3 / param.LU + param.Re, 0.01, deg2rad(0.1), 0, 0, 0.1];

% Convert to Cartesian
[rrC, vvC] = CoordConv.po2pv(kep_chief, param.mu);
[rrD, vvD] = CoordConv.po2pv(kep_deputy, param.mu);

x0C = [rrC; vvC]; % Chief initial state
x0D = [rrD; vvD]; % Deputy initial state

% Get the initial state in different coordinates 
[Relx0, ~] = Cart2Rel2(x0D, x0C);
Roex0 = kep2roe(kep_chief, kep_deputy);

%% Propagation Time
a0 = kep_chief(1);
Torb = 2*pi*sqrt(a0^3/param.mu); % Orbital period
param.t0 = 0; 
param.tf = Torb
% param.tf = 5;

paramArray = [param.mu, param.J2, param.Re, param.tf];

%% Propagate Chief and Deputy
% tspan = linspace(0, param.tf, 1000); % Two orbits
% opts = odeset('RelTol',1e-9,'AbsTol',1e-12);
% [~, Xchief]  = ode45(@(t,x) twoBodyODE(t,x,param.mu), tspan, X0C, opts);
% [~, Xdeputy] = ode45(@(t,x) twoBodyODE(t,x,mu), tspan, X0D, opts);
[fx0C, STM0C] = CloudProp(paramArray, x0C);
[fx0D, STM0D] = CloudProp(paramArray, x0D);

fx0Ckep   = CoordConv.vec2orbElem(fx0C(1:3),fx0C(4:6),param.mu);
fx0Dkep   = CoordConv.vec2orbElem(fx0D(1:3),fx0D(4:6),param.mu);

fx0Drel = Cart2Rel2(fx0D, fx0C);
fx0Droe = kep2roe(fx0Ckep, fx0Dkep)';

%% Function to compute Jacobians and predict Keplerian STM
STM0Dcart = STM0D % Perturb deputy only
STM0Drel = convertSTM(STM0Dcart, x0C, fx0C, x0D, fx0D, param.mu, 'hillframe');
STM0Droe = convertSTM(STM0Dcart, x0C, fx0C, x0D, fx0D, param.mu, 'roe');

%% Perturbed state propagation 

nPoints = 100; % number of pertubed states we are looking at. 

% The covariance of the pertubation introduced in cartesian coordinates. It
% is 1 km and 1 m/s here. 
r_var = 1e3/param.LU 
v_var = 1/param.LU*param.TU


for i = 1: nPoints

    i

    r_perturb = r_var * randn(1,3);  % position noise
    v_perturb = v_var * randn(1,3);  % velocity noise
    val(i,:) = [r_perturb, v_perturb]';

    % Update the cartesian initial state to the perturbed initial state
    xeC = x0C;
    xeD = x0D+val(i,:)';
    
    KepxeC  = CoordConv.vec2orbElem(xeC(1:3),xeC(4:6),param.mu);
    KepxeD  = CoordConv.vec2orbElem(xeD(1:3),xeD(4:6),param.mu);

    [Relxe, ~] = Cart2Rel2(xeD, xeC);
    Roexe = kep2roe(KepxeC, KepxeD);

    % In each coordinate system, whats the difference between the perturbed
    % initial state and nominal initial state? 
    valRel(i,:)   = Relxe - Relx0;
    valRoe(i,:)   = Roexe - Roex0;

    for j = 1:4;
        param.tf = j*Torb/3 - Torb/3;
        paramArray = [param.mu, param.J2, param.Re, param.tf];
    
        % Like before, we propagate the now perturbed initial state to get the
        % perturbed state transition matricies and the perturbed final states.
        [fxeC, STMeC] = CloudProp(paramArray, xeC);
        [fxeD, STMeD] = CloudProp(paramArray, xeD);
    
        fxeCkep = CoordConv.vec2orbElem(fxeC(1:3),fxeC(4:6),param.mu);
        fxeDkep = CoordConv.vec2orbElem(fxeD(1:3),fxeD(4:6),param.mu);
    
        fxeDrel = Cart2Rel2(fxeD, fxeC);
        fxeDroe = kep2roe(fxeCkep, fxeDkep)';
    
        % Now we calculate how much the perturbed state differes from the
        % nominal state after propagation in each coordinate system.
        dCartXXnlc(j, i,:) = fxeD - fx0D; 
        drelXXnlc(j, i,:) = fxeDrel - fx0Drel;
        droeXXnlc(j, i,:) = fxeDroe - fx0Droe;

    end
    
end

%% Analysis 

figure;
subplot(2,2,1); hold on;
p = plot(dCartXXnlc(1, :,1), dCartXXnlc(1, :,2), 'r.');
plot(squeeze(dCartXXnlc(1, :, 1)'), squeeze(dCartXXnlc(1, :, 2)'), 'k.');
plot_latex(p, 'X[Re]', 'Y[Re]', '', 'Initial', {'Cloud prop', 'Initial'});

subplot(2,2,2); hold on;
p = plot(dCartXXnlc(2, :,1), dCartXXnlc(2, :,2), 'r.');
plot(squeeze(dCartXXnlc(1, :, 1)'), squeeze(dCartXXnlc(1, :, 2)'), 'k.');
plot_latex(p, 'X[Re]', 'Y[Re]', '', '1/3 Orbit', {'Cloud prop', 'Initial'});

subplot(2,2,3); hold on;
p = plot(dCartXXnlc(3, :,1), dCartXXnlc(3, :,2), 'r.');
plot(squeeze(dCartXXnlc(1, :, 1)'), squeeze(dCartXXnlc(1, :, 2)'), 'k.');
plot_latex(p, 'X[Re]', 'Y[Re]', '', '2/3 Orbit', {'Cloud prop', 'Initial'});

subplot(2,2,4); hold on;
p = plot(dCartXXnlc(4, :,1), dCartXXnlc(4, :,2), 'r.');
plot(squeeze(dCartXXnlc(1, :, 1)'), squeeze(dCartXXnlc(1, :, 2)'), 'k.');
plot_latex(p, 'X[Re]', 'Y[Re]', '', '1 Orbit', {'Cloud prop', 'Initial'});

sgtitle('Cartesian ECI Cloudprop', 'Interpreter', 'latex');

figure;
subplot(2,2,1); hold on;
p = plot(drelXXnlc(1, :,1), drelXXnlc(1, :,2), 'r.');
plot(squeeze(drelXXnlc(1, :, 1)'), squeeze(drelXXnlc(1, :, 2)'), 'k.');
plot_latex(p, 'X[Re]', 'Y[Re]', '', 'Initial', {'Cloud prop', 'Initial'});

subplot(2,2,2); hold on;
p = plot(drelXXnlc(2, :,1), drelXXnlc(2, :,2), 'r.');
plot(squeeze(drelXXnlc(1, :, 1)'), squeeze(drelXXnlc(1, :, 2)'), 'k.');
plot_latex(p, 'X[Re]', 'Y[Re]', '', '1/3 Orbit', {'Cloud prop', 'Initial'});

subplot(2,2,3); hold on;
p = plot(drelXXnlc(3, :,1), drelXXnlc(3, :,2), 'r.');
plot(squeeze(drelXXnlc(1, :, 1)'), squeeze(drelXXnlc(1, :, 2)'), 'k.');
plot_latex(p, 'X[Re]', 'Y[Re]', '', '2/3 Orbit', {'Cloud prop', 'Initial'});

subplot(2,2,4); hold on;
p = plot(drelXXnlc(4, :,1), drelXXnlc(4, :,2), 'r.');
plot(squeeze(drelXXnlc(1, :, 1)'), squeeze(drelXXnlc(1, :, 2)'), 'k.');
plot_latex(p, 'X[Re]', 'Y[Re]', '', '1 Orbit 20', {'Cloud prop', 'Initial'});

sgtitle('Cartesian Hillframe Cloudprop', 'Interpreter', 'latex');

figure;
subplot(2,2,1); hold on;
p = plot(droeXXnlc(1, :,1), droeXXnlc(1, :,2), 'r.');
plot(squeeze(droeXXnlc(1, :, 1)'), squeeze(droeXXnlc(1, :, 2)'), 'k.');
plot_latex(p, 'delta a[Re]', 'delta lambda[rad]', '', 'Initial', {'Cloud prop', 'Initial'});

subplot(2,2,2); hold on;
p = plot(droeXXnlc(2, :,1), droeXXnlc(2, :,2), 'r.');
plot(squeeze(droeXXnlc(1, :, 1)'), squeeze(droeXXnlc(1, :, 2)'), 'k.');
plot_latex(p, 'delta a[Re]', 'delta lambda[rad]', '', '1/3 Orbit', {'Cloud prop', 'Initial'});

subplot(2,2,3); hold on;
p = plot(droeXXnlc(3, :,1), droeXXnlc(3, :,2), 'r.');
plot(squeeze(droeXXnlc(1, :, 1)'), squeeze(droeXXnlc(1, :, 2)'), 'k.');
plot_latex(p, 'delta a[Re]', 'delta lambda[rad]', '', '2/3 Orbit', {'Cloud prop', 'Initial'});

subplot(2,2,4); hold on;
p = plot(droeXXnlc(4, :,1), droeXXnlc(4, :,2), 'r.');
plot(squeeze(droeXXnlc(1, :, 1)'), squeeze(droeXXnlc(1, :, 2)'), 'k.');
plot_latex(p, 'delta a[Re]', 'delta lambda[rad]', '', '1 Orbit', {'Cloud prop', 'Initial'});

sgtitle('Relative Orbital Elements Cloudprop', 'Interpreter', 'latex');

% %% Compute covariance matrices and overlay 3D ellipsoids on point clouds with plot_latex
% 
% % Calculate covariance matrices using first 3 elements (position)
% C5 = cov(squeeze(dCartXXnlc(5,:,1:2)));
% C10 = cov(squeeze(dCartXXnlc(10,:,1:2)));
% C15 = cov(squeeze(dCartXXnlc(15,:,1:2)));
% C20 = cov(squeeze(dCartXXnlc(20,:,1:2)));
% 
% % Calculate means for centering ellipsoids
% mu5 = mean(squeeze(dCartXXnlc(5,:,1:2)),1);
% mu10 = mean(squeeze(dCartXXnlc(10,:,1:2)),1);
% mu15 = mean(squeeze(dCartXXnlc(15,:,1:2)),1);
% mu20 = mean(squeeze(dCartXXnlc(20,:,1:2)),1);
% 
% % Plot point clouds with overlaid ellipsoids in 2x2 subplots
% figure;
% subplot(2,2,1); hold on;
% p = plot(dCartXXnlc(5,:,1), dCartXXnlc(5,:,2), 'r.');
% plot(val(:,1), val(:,2), 'k.');
% plot_gaussian_ellipsoid(mu5, C5, 1, 'r');
% plot_latex(p, 'X [Re]', 'Y [Re]',  '5Torb', {'Propagated', 'Initial'});
% 
% subplot(2,2,2); hold on;
% p = plot(dCartXXnlc(10,:,1), dCartXXnlc(10,:,2), 'g.');
% plot(val(:,1), val(:,2), 'k.');
% plot_gaussian_ellipsoid(mu10, C10, 1, 'g');
% plot_latex(p, 'X [Re]', 'Y [Re]',  '10Torb', {'Propagated', 'Initial'});
% 
% subplot(2,2,3); hold on;
% p = plot(dCartXXnlc(15,:,1), dCartXXnlc(15,:,2), 'b.');
% plot(val(:,1), val(:,2), 'k.');
% plot_gaussian_ellipsoid(mu15, C15, 1, 'b');
% plot_latex(p, 'X [Re]', 'Y [Re]',  '15Torb', {'Propagated', 'Initial'});
% 
% subplot(2,2,4); hold on;
% p = plot(dCartXXnlc(20,:,1), dCartXXnlc(20,:,2),  'm.');
% plot(val(:,1), val(:,2), 'k.');
% plot_gaussian_ellipsoid(mu20, C20, 1, 'm');
% plot_latex(p, 'X [Re]', 'Y [Re]',  '20Torb', {'Propagated', 'Initial'});

%% Clean implementation: covariance propagation and quadratic interpolation + compute mus

% Define initial covariance matrix as before
r_var = 1e3 / param.LU;
v_var = 1 / param.LU * param.TU;
P0Cart = diag([r_var^2, r_var^2, r_var^2, v_var^2, v_var^2, v_var^2]);

% Get Jacobians from Cartesian to each frame
Jrel = jacobianCart2Hill(x0C, x0D);             % for Hill frame
Jroe = jacobianCart2Roe(x0C, x0D, param.mu);    % for ROEs

% Transform Cartesian covariance to Hill frame and ROE
P0Rel = Jrel * P0Cart * Jrel';
P0Roe = Jroe * P0Cart * Jroe';

% Calculate propagated final covariance
PfCart = STM0Dcart * P0Cart * STM0Dcart';
PfRel = STM0Drel * P0Rel * STM0Drel';
PfRoe = STM0Droe * P0Roe * STM0Droe';

% disp('Eigenvalues of P0:'), disp(eig(P0));
% disp('Eigenvalues of Pf:'), disp(eig(Pf));

% Define interpolation fractions
fractions = [0, 0.333, 0.666, 1.0];
PsCart = cell(1,4);
musCart = cell(1,4);
PsRel = cell(1,4);
musRel = cell(1,4);
PsRoe = cell(1,4);
musRoe = cell(1,4);

points = squeeze(dCartXXnlc(1, :, :));  % Now size is [nPoints x 6]

% Compute the mean of each column (i.e., each state variable)
mu = mean(points, 1);  % 1 x 6 row vector

% Compute interpolated covariance matrices and means (quadratic interpolation)
for i = 1:4
    frac = fractions(i);
    PsCart{i} = (1 - frac)^2 * P0Cart + frac^2 * PfCart;
%     musCart{i} = (1 - frac)^2 * mu + frac^2 * (fx0C-fx0D);
    musCart{i} = [-0.00015, 0.0015];
    PsRel{i} = (1 - frac)^2 * P0Rel + frac^2 * PfRel;
%     musRel{i} = (1 - frac)^2 * zeros(6,1) + frac^2 * (STM0Drel * zeros(6,1)); % linear propagation of zero mean
    musRel{i} = [0, 0.00152];
    PsRoe{i} = (1 - frac)^2 * P0Roe + frac^2 * PfRoe;
%     musRoe{i} = (1 - frac)^2 * zeros(6,1) + frac^2 * (STM0Droe * zeros(6,1)); % linear propagation of zero mean
    musRoe{i} = [-0.00015, 0.00134];
end

% Plot x-y position covariance ellipses for each interpolation
colors = {'b', 'b', 'b', 'b'};

figure;
titles = {'Initial', '1/3 Orbit', '2/3 Orbit', '1 Orbit'};
for i = 1:4
    subplot(2,2,i); hold on;
    Cxy = PsCart{i}(1:2,1:2); % extract x-y covariance
    mu_xy = musCart{i}(1:2)';
    plot_gaussian_ellipsoid(mu_xy, Cxy, 1, colors{i});
    p = plot(dCartXXnlc(i,:,1), dCartXXnlc(i,:,2),  'r.');
    plot(squeeze(dCartXXnlc(1, :, 1)'), squeeze(dCartXXnlc(1, :, 2)'), 'k.');
    plot_latex([], 'X [Re]', 'Y [Re]', '', titles{i}, {});
    axis equal; grid on;
    sgtitle('Cartesian ECI Linearised Covariance', 'Interpreter', 'latex');
end

figure;
titles = {'Initial', '1/3 Orbit', '2/3 Orbit', '1 Orbit'};
for i = 1:4
    subplot(2,2,i); hold on;
    Cxy = PsRel{i}(1:2,1:2); % extract x-y covariance
    mu_xy = musRel{i}(1:2)';
    plot_gaussian_ellipsoid(mu_xy, Cxy, 1, colors{i});
    p = plot(drelXXnlc(i,:,1), drelXXnlc(i,:,2),  'r.');
    plot(squeeze(drelXXnlc(1, :, 1)'), squeeze(drelXXnlc(1, :, 2)'), 'k.');
    plot_latex([], 'X [Re]', 'Y [Re]', '', titles{i}, {});
    axis equal; grid on;
    sgtitle('Cartesian Hillframe Linearised Covariance', 'Interpreter', 'latex');
end

figure;
titles = {'Initial', '1/3 Orbit', '2/3 Orbit', '1 Orbit'};
for i = 1:4
    subplot(2,2,i); hold on;
    Cxy = PsRoe{i}(1:2,1:2); % extract x-y covariance
    mu_xy = musRoe{i}(1:2)';
    plot_gaussian_ellipsoid(mu_xy, Cxy, 1, colors{i});
    p = plot(droeXXnlc(i,:,1), droeXXnlc(i,:,2),  'r.');
    plot(squeeze(droeXXnlc(1, :, 1)'), squeeze(droeXXnlc(1, :, 2)'), 'k.');
    plot_latex([], 'delta a [Re]', 'delta lamba [rad]', '', titles{i}, {});
    axis equal; grid on;
    sgtitle('ROE Linearised Covariance', 'Interpreter', 'latex');
end

% Helper function to plot 2D covariance ellipse (filled)
function plot_gaussian_ellipsoid(mu, C, sdwidth, color)
    [eigvec, eigval] = eig(C);
    sqrt_eigval = sqrt(max(diag(eigval), 0));
    theta = linspace(0, 2*pi, 100);
    ellipse_pts = eigvec * diag(sqrt_eigval) * [cos(theta); sin(theta)] * sdwidth;
    fill(ellipse_pts(1,:) + mu(1), ellipse_pts(2,:) + mu(2), color, 'FaceAlpha', 0.3, 'EdgeColor', color, 'LineWidth', 1.5);
end

function J = jacobianCart2Kep(x, mu)
    % Numerically estimate Jacobian d(kep)/d(cartesian)
    delta = 1e-6;
    J = zeros(6,6);
    for i = 1:6
        dx = zeros(6,1);
        dx(i) = delta;
        x_plus = x + dx;
        x_minus = x - dx;
        kep_plus = CoordConv.vec2orbElem(x_plus(1:3), x_plus(4:6), mu);
        kep_minus = CoordConv.vec2orbElem(x_minus(1:3), x_minus(4:6), mu);
        J(:, i) = (kep_plus - kep_minus) / (2*delta);
    end
end

function J = jacobianCart2Hill(xC, xD)
    % Numerically estimate Jacobian d(hill)/d(cartesian)
    delta = 1e-6;
    J = zeros(6,6);
    for i = 1:6
        dx = zeros(6,1);
        dx(i) = delta;
        xD_plus = xD + dx;
        xD_minus = xD - dx;
        [hill_plus, ~]  = Cart2Rel2([xD_plus(1:3)', xD_plus(4:6)']', xC);
        [hill_minus, ~] = Cart2Rel2([xD_minus(1:3)', xD_minus(4:6)']', xC);
        J(:, i) = (hill_plus - hill_minus) / (2*delta);
    end
end

function J = jacobianCart2Roe(xC, xD, mu)
    kep_chief = CoordConv.vec2orbElem(xC(1:3), xC(4:6), mu);
    % Numerically estimate Jacobian d(hill)/d(cartesian)
    delta = 1e-6;
    J = zeros(6,6);
    for i = 1:6
        dx = zeros(6,1);
        dx(i) = delta;
        xD_plus = xD + dx;
        xD_minus = xD - dx;
        kep_plus = CoordConv.vec2orbElem(xD_plus(1:3), xD_plus(4:6), mu);
        kep_minus = CoordConv.vec2orbElem(xD_minus(1:3), xD_minus(4:6), mu);
        roe_plus = kep2roe(kep_chief, kep_plus);
        roe_minus = kep2roe(kep_chief, kep_minus);
        J(:, i) = (roe_plus - roe_minus) / (2*delta);
    end
end


