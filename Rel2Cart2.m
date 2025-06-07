function [Cart,Rel2CartM] = Rel2Cart2(rho, CartChief)

%     
% [~, ~, xc(1,:)] = getPosition(td, 'H2AF15',3.986005e14,  param.J2, 6378137);
% xc(1) = xc(1)/ 6378137;
% 
% xcMEE = CoordConv.kepler2MEOE([xc(1), xc(2), xc(5),...
%     xc(4),xc(3),xc(6)]);
% 
% CartChief = CoordConv.ep2pv(xcMEE,param.mu);



%%  RTN matrix
rr = CartChief(1:3);
vv = CartChief(4:6);

hhat = cross(rr, vv)/norm(cross(rr, vv));
rhat = rr/norm(rr);

RTN = [rhat'; cross(hhat, rhat)' ; hhat'];


%% Positon 

rhoCart = RTN\rho(1:3);

CartDep = rhoCart+rr;

Cart(1:3) = CartDep;


%% rho dot

Chief_RTN = RTN*CartChief(1:3);
rc = Chief_RTN(1); 
omega_z = norm(cross(rr, vv))/rc^2;

gammadot = rho(6);
alphadot = rho(4)-omega_z*rho(2);
betadot = rho(5)+rho(1)*omega_z;


rhod_cart =  RTN\[alphadot;betadot;gammadot];

CartDepdot = rhod_cart+vv;


Cart(4:6) = CartDepdot;


Rel2CartM = inv([RTN, zeros(3,3); omega_z*([0 1 0 ; -1 0 0 ; 0 0 0])*RTN  RTN]);


%Rel2CartM = RTN'; 

%% rho ddot

