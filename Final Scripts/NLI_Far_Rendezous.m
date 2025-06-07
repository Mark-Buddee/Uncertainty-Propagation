%% Setup
clear,clc; close all;
dllFolder = 'C:\Program Files\MATLAB\R2022a\bin\win64'; % libmex.dll libMatlabDataArray.dll
daceFolder = 'C:\Program Files (x86)\DACE\lib'; % dace.dll

setenv('PATH', [getenv('PATH') ';' dllFolder ';' daceFolder]);

addpath('C:\Users\mbudd\Documents\thesis\UncertaintyProp-compiled');

clear; clc; close all;

%% Constants and Units
mu = 3.986005e14; % Earth gravitational parameter [m^3/s^2]
Re = 6378.137e3;  % Earth radius [m]
J2 =  1.08262668e-3;
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
kep_chief  = [500e3 / param.LU + param.Re, 0, deg2rad(97.4), 0, 0, 0];

% Convert to Cartesian
[rrC, vvC] = CoordConv.po2pv(kep_chief, param.mu);

x0C = [rrC; vvC] % Chief initial state

%% Generate list of trajectories into the future

%% Negate positions - now it goes back in time

%% Deputy Initial Conditions
rho = 20 / param.LU; % m
omega = sqrt(param.mu / kep_chief(1)^3);
v_rel = omega * rho;

dr = [rho; 0; 0];
dv = [0; -v_rel; v_rel]; %% Parking orbit

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
param.t0 = 0; 
param.tf = 1*Torb;

%% Non-linear index over time
% nOrbits = [1;11;21;31;41;51;61;71;81;91];
% nOrbits = [2;12;22;32;42;52];
% nOrbits = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
nOrbits = linspace(0.025, 4, 15)
% nOrbits = [1, 2, 3, 4, 5];
% nOrbits = [5, 10]
% nOrbits = [5];
for j = 1:length(nOrbits)
    param.tf = nOrbits(j)*Torb;

paramArray = [param.mu, param.J2, param.Re, param.tf];

%% Propagate Chief and Deputy
% tspan = linspace(0, param.tf, 1000); % Two orbits
% opts = odeset('RelTol',1e-9,'AbsTol',1e-12);
% [~, Xchief]  = ode45(@(t,x) twoBodyODE(t,x,param.mu), tspan, X0C, opts);
% [~, Xdeputy] = ode45(@(t,x) twoBodyODE(t,x,mu), tspan, X0D, opts);
[fx0C, STM0C] = CloudProp(paramArray, x0C);
[fx0D, STM0D] = CloudProp(paramArray, x0D);

fx0Ckep = CoordConv.vec2orbElem(fx0C(1:3), fx0C(4:6), param.mu);
fx0Dkep = CoordConv.vec2orbElem(fx0D(1:3), fx0D(4:6), param.mu);
fx0Drel = Cart2Rel2(fx0D, fx0C);
fx0Droe = kep2roe(fx0Ckep, fx0Dkep)';

%% Function to compute Jacobians and predict Keplerian STM
STM0Dcart = STM0D; % Perturb deputy only
STM0Drel = convertSTM(STM0Dcart, x0C, fx0C, x0D, fx0D, param.mu, 'hillframe');
STM0Droe = convertSTM(STM0Dcart, x0C, fx0C, x0D, fx0D, param.mu, 'roe');

%% Perturbed state propagation 

nPoints = 1000; % number of pertubed states we are looking at. 

% The covariance of the pertubation introduced in cartesian coordinates. It
% is 1 km and 1 m/s here. 
r_var = 1e3/param.LU ;
v_var = 1/param.LU*param.TU; 


for i = 1: nPoints

    r_perturb = r_var * randn(1,3);  % position noise
    v_perturb = v_var * randn(1,3);  % velocity noise
    val(i,:) = [r_perturb, v_perturb];

    % Update the cartesian initial state to the perturbed initial state
    xeC = x0C;
    xeD = x0D+val(i,:)';

    KepxeC   = CoordConv.vec2orbElem(xeC(1:3),xeC(4:6),param.mu);
    KepxeD   = CoordConv.vec2orbElem(xeD(1:3),xeD(4:6),param.mu);

    [Relxe, ~] = Cart2Rel2(xeD, xeC);
    Roexe = kep2roe(KepxeC, KepxeD);

    % In each coordinate system, whats the difference between the perturbed
    % initial state and nominal initial state? 
    valRel(i,:)   = Relxe - Relx0;
    valRoe(i,:)   = Roexe - Roex0;

%     valkep(i,:)   = KepxeD - Kepx0D;

    % Like before, we propagate the now perturbed initial state to get the
    % perturbed state transition matricies and the perturbed final states.
    [fxeC, STMeC] = CloudProp(paramArray, xeC);
    [fxeD, STMeD] = CloudProp(paramArray, xeD);

    fxeCkep = CoordConv.vec2orbElem(fxeC(1:3), fxeC(4:6), param.mu);
    fxeDkep = CoordConv.vec2orbElem(fxeD(1:3), fxeD(4:6), param.mu);
    fxeDrel = Cart2Rel2(fxeD, fxeC);
    fxeDroe = kep2roe(fxeCkep, fxeDkep)';

    % Compute Jacobians and predict Keplerian STM
    STMeDcart = STMeD; % Perturb deputy only
    STMeDrel = convertSTM(STMeDcart, xeC, fxeC, xeD, fxeD, param.mu, 'hillframe');
    STMeDroe = convertSTM(STMeDcart, xeC, fxeC, xeD, fxeD, param.mu, 'roe');

    % This is a calculaton of the nonlineartiy index based on the frobenius
    % norm from Eq. 5 in
    % https://link.springer.com/article/10.1007/BF03546420.
    nuCart(i) = norm(STMeDcart-STM0Dcart,'fro')/norm(STM0Dcart,'fro');
    nuRel(i)  = norm(STMeDrel-STM0Drel,'fro')/norm(STM0Drel,'fro');
    nuRoe(i)  = norm(STMeDroe-STM0Droe,'fro')/norm(STM0Droe,'fro');

    if(nuRoe(i) > 1e5)
        Roexe
        fxeDroe
        STMeDroe
        STM0Droe
    end

    % Now we calculate how much the perturbed state differes from the
    % nominal state after propagation in each coordinate system.
    dCartXXnlc(i,:) = fxeD - fx0D; 
    drelXXnlc(i,:) = fxeDrel - fx0Drel;
    droeXXnlc(i,:) = fxeDroe - fx0Droe;
    
end

%% Analysis 

% This plots the nonlinearity indicies - this is what you must compare for
% relative coordinate systems
nuCartl(j) = max(nuCart)
nuRell(j)  = max(nuRel)
nuRoel(j)  = max(nuRoe)

end

% nOrbits = [2, 12, 22, 32, 42, 52, 62, 72, 82, 92, 102];
% 
% nuCartl = [0.0244; 0.1554; 0.2928; 0.3249; 0.3754; 0.5790; 0.8470; 0.7796; 0.88; 0.92; 0.93];
% 
% nuMEEl = [0.0040; 0.0049; 0.0118; 0.0105; 0.0197; 0.0332; 0.0853; 0.0787];
% 
% nuGeql = 1.0e-15 * [0.1274; 0.1701; 0.1852; 0.1291; 0.1939; 0.1573; 0.1328; 0.1151];
% 
% nuCeql = 1.0e-15 * [0.1274; 0.1701; 0.1852; 0.1291; 0.1939; 0.1573; 0.1328; 0.1151];
% 
% nuKepl = [0.0029; 0.0042; 0.0123; 0.0109; 0.0203; 0.0339; 0.0862; 0.0794];
% 
% nuRell = [0.0235; 0.1472; 0.2782; 0.3083; 0.3573; 0.5520; 0.8160; 0.7564; 0.88; 0.92; 0.93];
% 
% nuRoel = [0.0164; 0.092; 0.1230; 0.1096; 0.1171; 0.1619; 0.2307; 0.1226; 0.2; 0.23; 0.4];

nuCartl(nuCartl <= 0) = 1e-12;
nuRell(nuRell <= 0) = 1e-12;
nuRoel(nuRoel <= 0) = 1e-12;

figure; hold on;
p(1) = plot(nOrbits, nuCartl, '-');
p(6) = plot(nOrbits, nuRell, '-');
p(7) = plot(nOrbits, nuRoel, '-');

set(gca, 'YScale', 'log');
% ylim([1e-12, 1e2]); % adjust upper limit as needed

plot_latex(p, 'Number of Orbits', 'Non-linearity Index', '', '', {'Cartesian', 'Hillframe', 'ROE'});

% 3D plot for Cartesian uncertainty distribution
figure;
hold on;
p = plot3(dCartXXnlc(:,1), dCartXXnlc(:,2), dCartXXnlc(:,3), 'r.');
plot3(val(:,1), val(:,2), val(:,3), 'k.');
% plot_latex(p, 'X [Re]', 'Y [Re]', 'Z [Re]', 'Cartesian 3D Distribution', {'Cloud prop', 'Initial'});
xlabel('X [Re]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Y [Re]', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('Z [Re]', 'Interpreter', 'latex', 'FontSize', 14);
title('Cartesian 3D Distribution', 'Interpreter', 'latex', 'FontSize', 16);
legend({'Cloud prop', 'Initial'}, 'Interpreter', 'latex', 'FontSize', 12);
grid on;
% axis equal;
view(30,30); % sets a nice 3D angle

% This is a distribution of uncertainity - a linear trend here is good,
% obviously the flatter the red curve is the better/less uncertainity as
% well. 
figure; 
subplot(1,5,1);hold on; 
p = plot(dCartXXnlc(:,1), dCartXXnlc(:,2), 'r.');
% p = plot(dCartXX(:,1), dCartXX(:,2), 'b.');
plot(val(:,1), val(:,2), 'k.');
plot_latex(p, 'X[Re]', 'Y[re]','', 'Cartesian' ,{'Cloud prop', 'Initial' });

% subplot(1,5,2);hold on; 
% p = plot(dMEEXXnlc(:,1), dMEEXXnlc(:,6), 'r.');
% % p = plot(dMEEXX(:,1), wrapToPi(dMEEXX(:,6)), 'b.');
% plot(valMEE(:,1), valMEE(:,6), 'k.');
% plot_latex(p, 'p [Re]', 'L[rad]','', 'Modified Equinoctials' ,{});
% 
% subplot(1,5,3);hold on; 
% p = plot(dGeqXXnlc(:,1), dGeqXXnlc(:,4), 'r.');
% % p = plot(dGEqoeXX(:,1), wrapToPi(dGEqoeXX(:,4)), 'b.');
% plot(valGeqOE(:,1), valGeqOE(:,4), 'k.');
% plot_latex(p, '\nu', 'L[rad]','', 'GEqoOE' ,{});
% 
% subplot(1,5,4);hold on; 
% p = plot(dCeqXXnlc(:,1), dCeqXXnlc(:,4), 'r.');
% % p = plot(dCEqoeXX(:,1), wrapToPi(dCEqoeXX(:,4)), 'b.');
% plot(valCeqOE(:,1), valCeqOE(:,4), 'k.');
% plot_latex(p, '\nu', 'L[rad]','', 'CEqoOE' ,{});

% New: 15 subplots for drelXXnlc, overlay valRel
figure;
counter = 1;
for i = 1:5
    for j = i+1:6
        subplot(3,5,counter); hold on;
        plot(drelXXnlc(:,i), drelXXnlc(:,j), 'r.');
        plot(valRel(:,i), valRel(:,j), 'k.');
        xlabel(['Rel ', num2str(i)], 'Interpreter', 'latex', 'FontSize', 10);
        ylabel(['Rel ', num2str(j)], 'Interpreter', 'latex', 'FontSize', 10);
        grid on;
        counter = counter + 1;
    end
end
sgtitle('Relative Coordinates Perturbations', 'Interpreter', 'latex');

% New: 15 subplots for droeXXnlc, overlay valRoe
figure;
counter = 1;
for i = 1:5
    for j = i+1:6
        subplot(3,5,counter); hold on;
        plot(droeXXnlc(:,i), droeXXnlc(:,j), 'r.');
        plot(valRoe(:,i), valRoe(:,j), 'k.');
        xlabel(['ROE ', num2str(i)], 'Interpreter', 'latex', 'FontSize', 10);
        ylabel(['ROE ', num2str(j)], 'Interpreter', 'latex', 'FontSize', 10);
        grid on;
        counter = counter + 1;
    end
end
sgtitle('Relative Orbital Elements Perturbations', 'Interpreter', 'latex');

% %% Visual Validation
% 
% figure;
% subplot(1,2,1);
% plot3(Xchief(:,1), Xchief(:,2), Xchief(:,3), 'b', 'LineWidth', 1.2); hold on;
% plot3(Xdeputy(:,1), Xdeputy(:,2), Xdeputy(:,3), 'r', 'LineWidth', 1.2);
% legend('Chief', 'Deputy');
% title('Orbits in Inertial Frame'); axis equal; grid on;
% xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
% 
% % subplot(1,2,2);
% % plot3(rel_LVLH(:,1), rel_LVLH(:,2), rel_LVLH(:,3), 'k', 'LineWidth', 1.2);
% % title('Deputy Relative to Chief in LVLH Frame');
% % xlabel('Radial (x) [m]'); ylabel('Along-track (y) [m]'); zlabel('Cross-track (z) [m]');
% % axis equal; grid on;
