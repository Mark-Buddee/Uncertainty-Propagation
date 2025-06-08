% Supervisor provided

function [rho_RTN,Cart2RelM] = Cart2Rel2(CartDep,CartChief)

rr = CartChief(1:3);
vv = CartChief(4:6);

%%  RTN matrix

rhat = rr/norm(rr);
hhat = cross(rr, vv)/norm(cross(rr, vv));

RTN = [rhat'; cross(hhat, rhat)' ; hhat'];
%% rho 
rho_cart = CartDep - CartChief;
rho = RTN*rho_cart(1:3);
rho_RTN(1:3,1) = rho;

%% rho dot

Chief_RTN = RTN*CartChief(1:3);
rc = Chief_RTN(1); 
omega_z = norm(cross(rr, vv))/rc^2;

vv_RTN =  RTN*rho_cart(4:6);

zdot = vv_RTN(end);
xdot = vv_RTN(1) + omega_z* rho(2);
ydot =  vv_RTN(2) - omega_z* rho(1);

rho_RTN(4:6,1) = [xdot;ydot;zdot];


Cart2RelM = [RTN, zeros(3,3); omega_z*([0 1 0 ; -1 0 0 ; 0 0 0])*RTN  RTN];
end