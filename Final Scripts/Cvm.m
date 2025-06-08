%% Setup
clear,clc; close all;
dllFolder = 'C:\Program Files\MATLAB\R2022a\bin\win64'; % libmex.dll libMatlabDataArray.dll
daceFolder = 'C:\Program Files (x86)\DACE\lib'; % dace.dll
setenv('PATH', [getenv('PATH') ';' dllFolder ';' daceFolder]);

clear; clc; close all;

addpath('C:\Users\mbudd\Documents\thesis\UncertaintyProp-compiled');

rng(1);

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

%% Define Chief (Keplerian Elements)
% Format: [a, e, i, RAAN, omega, nu]
kep_chief  = [500e3 / param.LU + param.Re, 0.001, deg2rad(97.4), deg2rad(0), deg2rad(45), deg2rad(-45)];
% Convert to Cartesian
[rrC, vvC] = CoordConv.po2pv(kep_chief, param.mu);

x0C = [rrC; vvC] % Chief initial state

%% Deputy Initial Conditions
rho = 20 / param.LU; % m
omega = sqrt(param.mu / kep_chief(1)^3);
v_rel = omega * rho;

% dr = [rho; 0; 0]; %% Parking orbit
% dr = [rho; -10e3/param.LU; 0]; %% Initial approach
dr = [-1.4307e3/param.LU; -1.4022e5/param.LU; 0];

% dv = [0; -v_rel; v_rel]; %% Parking orbit
% dv = [0; -1.2*v_rel; v_rel]; %% Initial approach
dv = [-0.0657/param.VU; 3.8066/param.VU; 0];

% Compute RTN frame from chief state
rhat = rrC / norm(rrC);
hhat = cross(rrC, vvC) / norm(cross(rrC, vvC));
RTN = [rhat'; cross(hhat, rhat)'; hhat'];

% Convert perturbation to ECI
rho_RTN = [dr; dv];
perturbation_ECI = [RTN' zeros(3); zeros(3) RTN'] * rho_RTN;

% Apply to chief to get deputy
x0D = x0C + perturbation_ECI

kep_deputy = CoordConv.vec2orbElem(x0D(1:3), x0D(4:6), param.mu);

%%

% Get the initial state in different coordinates 
[Relx0, ~] = Cart2Rel2(x0D, x0C);
Roex0 = kep2roe(kep_chief, kep_deputy);

%% Propagation Time
a0 = kep_chief(1);
Torb = 2*pi*sqrt(a0^3/param.mu); % Orbital period
% stm_time = 5; % 5s for the STM Generation
% param.tf = stm_time/param.TU;
param.tf = Torb;
stm_time = Torb*param.TU;

paramArray = [param.mu, param.J2, param.Re, param.tf];

%% Propagate Chief and Deputy
[fx0C, STM0C] = CloudProp(paramArray, x0C);
[fx0D, STM0D] = CloudProp(paramArray, x0D);

fx0Ckep   = CoordConv.vec2orbElem(fx0C(1:3),fx0C(4:6),param.mu);
fx0Dkep   = CoordConv.vec2orbElem(fx0D(1:3),fx0D(4:6),param.mu);

fx0Drel = Cart2Rel2(fx0D, fx0C);
fx0Droe = kep2roe(fx0Ckep, fx0Dkep)';

%% Function to compute Jacobians and predict Keplerian STM
STM0Dcart = STM0D; % Perturb deputy only
STM0Drel = convertSTM(STM0Dcart, x0C, fx0C, x0D, fx0D, param.mu, 'hillframe');
STM0Droe = convertSTM(STM0Dcart, x0C, fx0C, x0D, fx0D, param.mu, 'roe');

%% Get nominal path

nSteps = 8;
% nSteps = 13
% skips_per_fig = nSteps/5;
skips_per_fig = 1;

currXC = x0C;
currXD = x0D;
for i = 1:nSteps
    [fx0C, ~] = CloudProp(paramArray, currXC);
    [fx0D, ~] = CloudProp(paramArray, currXD);

    fx0Ckep = CoordConv.vec2orbElem(fx0C(1:3),fx0C(4:6),param.mu);
    fx0Dkep = CoordConv.vec2orbElem(fx0D(1:3),fx0D(4:6),param.mu);

    fx0Drel = Cart2Rel2(fx0D, fx0C);
    fx0Droe = kep2roe(fx0Ckep, fx0Dkep)';

    % Now we calculate how much the perturbed state differs from the
    % nominal state after propagation in each coordinate system.
    nomcartLocs(i,:) = fx0D - fx0C;
    nomrelLocs(i,:) = fx0Drel;
    nomroeLocs(i,:) = fx0Droe;

    currXC = fx0C;
    currXD = fx0D;
end

%% Perturbed state propagation 

nPoints = 100; % number of pertubed states we are looking at.

% The covariance of the pertubation introduced in cartesian coordinates. It
% is 1 km and 1 m/s here. 
% r_var = 1e3/param.LU;
% v_var = 1/param.LU*param.TU;
r_var = 10/param.LU;
v_var = 0.01/param.LU*param.TU;
% r_var = 5000/param.LU;
% v_var = 5/param.LU*param.TU;
% r_var = 50000/param.LU;
% v_var = 50/param.LU*param.TU;

% skips_per_fig = 1;
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

    currXC = xeC;
    currXD = xeD;
    for j = 1:nSteps

        % Like before, we propagate the now perturbed initial state to get the
        % perturbed state transition matricies and the perturbed final states.
        [fxeC, ~] = CloudProp(paramArray, currXC);
        [fxeD, ~] = CloudProp(paramArray, currXD);
    
        fxeCkep = CoordConv.vec2orbElem(fxeC(1:3),fxeC(4:6),param.mu);
        fxeDkep = CoordConv.vec2orbElem(fxeD(1:3),fxeD(4:6),param.mu);
    
        fxeDrel = Cart2Rel2(fxeD, fxeC);
        fxeDroe = kep2roe(fxeCkep, fxeDkep)';
    
        % Now we calculate how much the perturbed state differs from the
        % nominal state after propagation in each coordinate system.
        cartLocs(j, i,:) = fxeD - fxeC;
        relLocs(j, i,:) = fxeDrel;
        roeLocs(j, i,:) = fxeDroe;

        currXC = fxeC;
        currXD = fxeD;

    end
    
end

%% Analysis 

% Define initial covariance matrix as before
PCart = diag([r_var^2, r_var^2, r_var^2, v_var^2, v_var^2, v_var^2]);

% Get Jacobians from Cartesian to each frame
Jrel = jacobianCart2Hill(x0C, x0D);             % for Hill frame
Jroe = jacobianCart2Roe(x0C, x0D, param.mu);    % for ROEs

% Transform Cartesian covariance to Hill frame and ROE
PRel = Jrel * PCart * Jrel';
PRoe = Jroe * PCart * Jroe';

% dcart = fx0D - x0D;
% drel = fx0Drel - Relx0;
% droe = fx0Droe - Roex0;

% Define interpolation fractions
PsCart = cell(1,4);
musCart = cell(1,4);
PsRel = cell(1,4);
musRel = cell(1,4);
PsRoe = cell(1,4);
musRoe = cell(1,4);

% points = squeeze(dCartXXnlc(1, :, :));  % Now size is [nPoints x 6]
% 
% % Compute the mean of each column (i.e., each state variable)
% mu = mean(points, 1);  % 1 x 6 row vector

% for i = 1:4
%     PsCart{i} = PCart;
%     (i-1)*skips_per_fig + 1
%     musCart{i} = nomcartLocs((i-1)*skips_per_fig + 1, :);
% %     musCart{i} = [0, 0, 0, 0];
% 
%     PsRel{i} = PRel;
%     musRel{i} = nomrelLocs((i-1)*skips_per_fig + 1, :);
% %     musRel{i} = [0, 0, 0, 0];
% 
%     PsRoe{i} = PRoe;
%     musRoe{i} = nomroeLocs((i-1)*skips_per_fig + 1, :);
% %     musRoe{i} = [0, 0, 0, 0];
% 
%     % Move it along 5 seconds
%     PCart = STM0Dcart^skips_per_fig * PCart * STM0Dcart'^skips_per_fig;
%     PRel = STM0Drel^skips_per_fig * PRel * STM0Drel'^skips_per_fig;
%     PRoe = STM0Droe^skips_per_fig * PRoe * (STM0Droe'^skips_per_fig);
% end

% PCart
% PRel
% PRoe
% STM0Droe 

for i = 1:nSteps
    PsCart{i} = PCart;
    musCart{i} = nomcartLocs(i, :);

    PsRel{i} = PRel;
    musRel{i} = nomrelLocs(i, :);

    PsRoe{i} = PRoe;
    musRoe{i} = nomroeLocs(i, :);

    % Move it along 5 seconds
    PCart = STM0Dcart * PCart * STM0Dcart';
    PRel = STM0Drel * PRel * STM0Drel';
    PRoe = STM0Droe * PRoe * (STM0Droe');
end

% Plot x-y position covariance ellipses for each interpolation
colors = {'b', 'b', 'b', 'b'};

% figure;
% titles = {'Initial', '5s', '10s', '15s'};
% tiledlayout(2,2)
% for i = 1:4
%     j = (i-1)*skips_per_fig + 1;
%     ax = nexttile; hold on;
%     view(ax, 3); % <--- Make this axes 3D explicitly
% 
%     Cxyz = PsCart{j}(1:3,1:3);  % extract x-y-z covariance
%     mu_xyz = musCart{j}(1:3)';
%     plot_gaussian_ellipsoid_3d(mu_xyz, Cxyz, 1, colors{i});
%     plot3(cartLocs(j,:,1), cartLocs(j,:,2), cartLocs(j,:,3), 'r.');
% 
%     plot_latex([], 'X [Re]', 'Y [Re]', 'Z [Re]', titles{i}, {});
%     axis equal; grid on;
% end
% sgtitle('Cartesian ECI Linearised Covariance (3D)', 'Interpreter', 'latex');

% figure;
% titles = {'Initial', '1/3 Orbit', '2/3 Orbit', '1 Orbit'};
% tiledlayout(2,2)
% for i = 1:4
%     j = (i-1)*skips_per_fig + 1;
%     ax = nexttile; hold on;
%     view(ax, 3);
% 
%     Cxyz = PsRel{j}(1:3,1:3); % Full 3D covariance
%     mu_xyz = musRel{j}(1:3)';
%     plot_gaussian_ellipsoid_3d(mu_xyz, Cxyz, 1, colors{i});
%     plot3(relLocs(j,:,1), relLocs(j,:,2), relLocs(j,:,3), 'r.');
% 
%     plot_latex([], 'X [Re]', 'Y [Re]', 'Z [Re]', titles{i}, {});
%     axis equal; grid on;
% end
% sgtitle('Cartesian Hillframe Linearised Covariance (3D)', 'Interpreter', 'latex');
% 
% figure;
% titles = {'Initial', '1/3 Orbit', '2/3 Orbit', '1 Orbit'};
% tiledlayout(2,2)
% for i = 1:4
%     j = (i-1)*skips_per_fig + 1;
%     ax = nexttile; hold on;
%     view(ax, 3);
% 
%     Cxyz = PsRoe{j}(1:3,1:3); % Use full 3D covariance
%     mu_xyz = musRoe{j}(1:3)';
%     plot_gaussian_ellipsoid_3d(mu_xyz, Cxyz, 1, colors{i});
%     plot3(roeLocs(j,:,1), roeLocs(j,:,2), roeLocs(j,:,3), 'r.');
% 
%     plot_latex([], '\delta a [Re]', '\delta \lambda [rad]', '\delta e_x', titles{i}, {});
%     axis equal; grid on;
% end
% sgtitle('ROE Linearised Covariance (3D)', 'Interpreter', 'latex');


dim = 6;
nMC = 20;

% Preallocate result arrays
CvMstats = struct('Cart', zeros(1,nSteps), 'Rel', zeros(1,nSteps), 'Roe', zeros(1,nSteps));
pValues  = struct('Cart', zeros(1,nSteps), 'Rel', zeros(1,nSteps), 'Roe', zeros(1,nSteps));

systems = {'Cart', 'Rel', 'Roe'};
musAll = {musCart, musRel, musRoe};
PsAll  = {PsCart, PsRel, PsRoe};
locsAll = {cartLocs, relLocs, roeLocs};

for s = 1:3
    mus = musAll{s};
    Ps = PsAll{s};
    locs = locsAll{s};
    tag = systems{s};

    for i = 1:nSteps
        i
        mu = mus{i}(:);
        C = Ps{i};
        
        % Extract point cloud: [nPoints × 6]
        X = squeeze(locs(i,:,:));

        % Mahalanobis distances
        diffs = X - mu';  % each row: x_i - mu
        d2 = sum((diffs / C) .* diffs, 2);  % properly computes xᵢᵀ C⁻¹ xᵢ

        % Empirical vs theoretical CDF
        [F_emp, x_emp] = ecdf(d2);
        F_theor = chi2cdf(x_emp, dim);

        % CvM statistic
        w2 = sum((F_emp - F_theor).^2);
        CvMstats.(tag)(i) = w2;

        % Monte Carlo null distribution
        nullStats = zeros(nMC,1);
        for j = 1:nMC
            samp = chi2rnd(dim, 1000, 1);  % match sample size
            [Fnull, xnull] = ecdf(samp);
            Ftrue = chi2cdf(xnull, dim);
            nullStats(j) = sum((Fnull - Ftrue).^2);
        end

        % Empirical p-value
        pValues.(tag)(i) = mean(nullStats > w2);
    end
end
%% Plotting CVM statistics with LaTeX formatting
figure;
p = plot(stm_time*(1:nSteps), CvMstats.Cart, 'b-', ...
         stm_time*(1:nSteps), CvMstats.Rel,  'r-', ...
         stm_time*(1:nSteps), CvMstats.Roe,  'g-');
plot_latex(p, 'Time (s)', 'Cramér-von Mises Statistic', '', ...
           '', ...
           {'Cartesian', 'Hillframe', 'ROE'});

%% Plotting p-values with LaTeX formatting
figure;
p = plot(stm_time*(1:nSteps), pValues.Cart, 'b-', ...
         stm_time*(1:nSteps), pValues.Rel,  'r-', ...
         stm_time*(1:nSteps), pValues.Roe,  'g-'); hold on;
yline(0.05, '--k', 'Significance Level (0.05)', 'Interpreter', 'tex', ...
       'LabelHorizontalAlignment', 'left', 'FontSize', 10);
plot_latex(p, 'Time (s)', 'p-value', '', ...
           '', ...
           {'Cartesian', 'Hillframe', 'ROE'});
ylim([0 1]);


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
    delta_frac = 1e-6; % fractional perturbation, typical value
    J = zeros(6,6);
    for i = 1:6
        dx = zeros(6,1);
        dx(i) = delta_frac * max(abs(x(i)), 1e-3); % avoid tiny scaling near 0
        x_plus = x + dx;
        x_minus = x - dx;
        kep_plus = CoordConv.vec2orbElem(x_plus(1:3), x_plus(4:6), mu);
        kep_minus = CoordConv.vec2orbElem(x_minus(1:3), x_minus(4:6), mu);
        J(:, i) = (kep_plus - kep_minus) / (2*dx(i));
    end
end

function J = jacobianCart2Hill(xC, xD)
    % Numerically estimate Jacobian d(hill)/d(cartesian)
    delta_frac = 1e-6; % fractional perturbation, typical value
    J = zeros(6,6);
    for i = 1:6
        dx = zeros(6,1);
        dx(i) = delta_frac * max(abs(xD(i)), 1e-3); % avoid tiny scaling near 0
        xD_plus = xD + dx;
        xD_minus = xD - dx;
        [hill_plus, ~]  = Cart2Rel2([xD_plus(1:3)', xD_plus(4:6)']', xC);
        [hill_minus, ~] = Cart2Rel2([xD_minus(1:3)', xD_minus(4:6)']', xC);
        J(:, i) = (hill_plus - hill_minus) / (2*dx(i));
    end
end

function J = jacobianCart2Roe(xC, xD, mu)
    kep_chief = CoordConv.vec2orbElem(xC(1:3), xC(4:6), mu);
    delta_frac = 1e-6; % fractional perturbation, typical value
    J = zeros(6,6);
    for i = 1:6
        dx = zeros(6,1);
        dx(i) = delta_frac * max(abs(xD(i)), 1e-3); % avoid tiny scaling near 0
        xD_plus = xD + dx;
        xD_minus = xD - dx;
        kep_plus = CoordConv.vec2orbElem(xD_plus(1:3), xD_plus(4:6), mu);
        kep_minus = CoordConv.vec2orbElem(xD_minus(1:3), xD_minus(4:6), mu);
        roe_plus = kep2roe(kep_chief, kep_plus);
        roe_minus = kep2roe(kep_chief, kep_minus);
        d_roe = roe_plus - roe_minus;
        d_roe(6) = wrapToPi(roe_plus(6)) - wrapToPi(roe_minus(6));
        J(:, i) = d_roe / (2 * dx(i));
    end
end

function plot_gaussian_ellipsoid_3d(mu, C, sdwidth, color)
% PLOT_GAUSSIAN_ELLIPSOID_3D plots a 3D Gaussian ellipsoid given a mean and covariance.
%
% INPUTS:
%   mu      - 3x1 mean vector
%   C       - 3x3 covariance matrix
%   sdwidth - width of the ellipsoid in std deviations (1 = 68%, 2 = 95%)
%   color   - color string or RGB triplet for surface and edge
%
% USAGE EXAMPLE:
%   plot_gaussian_ellipsoid_3d([0;0;0], eye(3), 1, 'b');

    % Ensure inputs are column vectors/matrices
    mu = mu(:); 
    assert(length(mu) == 3, 'mu must be a 3x1 vector.');
    assert(all(size(C) == [3 3]), 'C must be a 3x3 covariance matrix.');

    % Eigen-decomposition of the covariance matrix
    [eigvec, eigval] = eig(C);
    
    % Generate a unit sphere
    [x, y, z] = sphere(30); % More points = smoother
    xyz = [x(:) y(:) z(:)]';

    % Scale the unit sphere using the eigenvalues and vectors
    scaling = sdwidth * eigvec * sqrt(eigval);
    ellipsoid_pts = scaling * xyz + mu;

    % Reshape to match sphere size
    X = reshape(ellipsoid_pts(1,:), size(x));
    Y = reshape(ellipsoid_pts(2,:), size(y));
    Z = reshape(ellipsoid_pts(3,:), size(z));

    % Plot the ellipsoid
    surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', color, 'FaceColor', color);
    hold on;
    plot3(mu(1), mu(2), mu(3), 'k.', 'MarkerSize', 15); % Mean

    % Annotate
    plot_latex([], '$x$', '$y$', '$z$', '3D Gaussian Ellipsoid', {});
    axis equal; grid on;
end