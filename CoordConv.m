% Mostly supervisor provided

%Coordinate conversions
classdef CoordConv
    methods (Static)

        function MEEParameters = kepler2MEOE(Object)
            % input: A struct or array with keplerian orbital elements
            % output: A vector of MEE elements.
            
            if isstruct(Object)
                
                a = Object.a;
                e = Object.e;
                RAAN = Object.RAAN;
                AOP = Object.AOP;
                i = Object.i;
                
                if isfield(Object, 'MA')
                    MA = Object.MA;
                    
                    
                    E = MA;
                    E_old = 1;
                    precision  = 1e-7;
                    
                    
                    while abs(E - E_old) > precision
                        
                        % find the next eccentric anomaly value using N-R method
                        E = E - ((E - e*sin(E) - MA)./(1 - e*cos(E)));
                        E_old = E ;
                    end
                    
                    theta = 2 * atan(sqrt((1 + e)/ (1 - e)) .* tan(E/2) );
                else
                    theta = Object.theta;
                end
                
                
            else
                a = Object(1);
                e = Object(2);
                RAAN = Object(3);
                AOP = Object(4);
                i = Object(5);
                theta = Object(6);
                
            end
            
            p = a*(1 - e*e);
            f = e*cos(RAAN + AOP);
            g = e*sin(RAAN + AOP);
            h = tan(i/2)*cos(AOP);
            k = tan(i/2)*sin(AOP);
            L = wrapTo2Pi(AOP + RAAN + theta);
            MEEParameters = [p, f,g ,h,k,L];
            
            
            
        end

        function [a,e,INC,AOP,RAAN,MA,TA] = RV2OE(R,V,mu)
            %--------------------------------------------------------------------------
            %   Converts R and V to classical orbital elements
            %--------------------------------------------------------------------------
            %   Form:
            %   [a,e,INC,AOP,RAAN,MA,TA] = RV2OE( R, V, mu)
            %--------------------------------------------------------------------------
            %
            %   ------
            %   Inputs
            %   ------
            %   R               (1,3) Position vector
            %   V               (1,3) Velocity vector
            %   mu                    Gravitational parameter
            %
            %   -------
            %   Outputs
            %   -------
            %   X               (1,6) Elements vector [a,e,INC,AOP,RAAN,MA,TA]
            %
            % Magnitude of vectors (lowercase)
            r = norm(R);    v = norm(V);
            % Calculate h and e
            H = cross(R,V);     % angular momentum
            h = norm(H);
            E = -R/r - 1/mu*cross(H,V);   % eccentricity
            e = norm(E);
            % Perifocal frame
            u1 = E/e;
            u3 = H/h;
            u2 = cross(u3,u1);
            % Inertial frame
            i = [1 0 0];
            j = [0 1 0];
            k = [0 0 1];
            uln = cross(k,u3);    % node line
            % Angles from perifocal frame
            INC = acos(dot(u3,k));    % inclination
            RAAN = atan2(dot(uln,j),dot(uln,i));    % ascending node
            if RAAN < 0
                RAAN = RAAN + 2*pi;
            end
            AOP = atan2(dot(-u2,uln),dot(uln,u1));  % argument of perigee
            if AOP < 0
                AOP = AOP + 2*pi;
            end
            TA = atan2(dot(R,u2),dot(R,u1));    % true anomaly
            if TA < 0
                TA = TA + 2*pi;
            end
            u = 2*atan(sqrt((1-e)/(1+e))*tan(TA/2));    % eccentric anomaly
            MA = u - e*sin(u);  % mean anomaly (Kepler's eq.)
            % Add semimajor axis from energy equation
            a = r/(2 - v^2*r/mu);
        end
        function X = RV2CEq(R,V,mu, J2, Re)
            %--------------------------------------------------------------------------
            %   Converts R and V to Equinoctial elements
            %--------------------------------------------------------------------------
            %   Form:
            %   eq = RV2GEq( R, V, mu ,U)
            %--------------------------------------------------------------------------
            %
            %   ------
            %   Inputs
            %   ------
            %   R               (1,3) Position vector
            %   V               (1,3) Velocity vector
            %   mu                    Gravitational parameter

            %   -------
            %   Outputs
            %   -------
            %   X               (1,6) Elements vector [nu,p1,p2,varL,q1,Q2]
            %
            %--------------------------------------------------------------------------
            %   References:  Baù, G., Hernando-Ayuso, J. and Bombardelli, C., 2021.
            %   "A generalization of the equinoctial orbital elements"
            %   Celestial Mechanics and Dynamical Astronomy, 133(11), pp.1-29.
            %--------------------------------------------------------------------------


            %% Calculate GEqOE from R,V
            [~,~,INC,~,RAAN,~,~] = CoordConv.RV2OE(R,V,mu);  % obtain RAAN and INC
            % State parameters
            r = norm(R);
            v = norm(V);
            h = norm(cross(R,V));           % angular momentum module
            r_dot = dot(R,V)/r;             % radial velocity
            % Calculate energy
            epsK = 1/2*v^2 - mu/r;          % Keplerian energy


            U = 0;

            eps = epsK + U;                 % total energy
            nu = 1/mu*(-2*eps)^(3/2);       % (1)
            % Calculate q in order to obtain equinoctial frame
            q1 = tan(INC/2)*sin(RAAN);  % (5)
            q2 = tan(INC/2)*cos(RAAN);  % (6)
            % Unit vectors of equinoctial frame in inertial ref.
            Kq = 1/(1+q1^2+q2^2);
            ex = Kq*[1-q1^2+q2^2    ,2*q1*q2        ,-2*q1];
            ey = Kq*[2*q1*q2        ,1+q1^2-q2^2    ,2*q2];
            % Radial vector in inertial ref.
            er = R/r;
            % Obtain true longitude
            cL = dot(er,ex);
            sL = dot(er,ey);
            L = atan2(sL,cL);
            if L < 0
                L = L + 2*pi;
            end
            % L = AOP + RAAN + TA;
            % Calculate rho,c and a
            Ueff = h^2/2/r^2 + U;   % effective potential
            c = sqrt(2*r^2*Ueff);
            rho = c^2/mu;
            a = -mu/2/eps;
            % g vector
            p1 = (rho/r-1)*sin(L) - c*r_dot/mu*cos(L);  % (2)
            p2 = (rho/r-1)*cos(L) + c*r_dot/mu*sin(L);  % (3)
            % Obtain K for varL
            w = sqrt(mu/a);
            sK = (mu+c*w-r*r_dot^2)*sin(L) - r_dot*(c+w*r)*cos(L);
            cK = (mu+c*w-r*r_dot^2)*cos(L) + r_dot*(c+w*r)*sin(L);
            K = atan2(sK,cK);
            if K < 0
                K = K + 2*pi;
            end
            % Obtain generalized mean longitud varL
            varL = K + 1/(mu+c*w)*(cK*p1-sK*p2);  % (4)
            % GEqOE array
            X = [nu p1 p2 varL q1 q2];
        end

        function X = RV2GEq(R,V,mu,J2, Re)
            %--------------------------------------------------------------------------
            %   Converts R and V to Equinoctial elements
            %--------------------------------------------------------------------------
            %   Form:
            %   eq = RV2GEq( R, V, mu ,U)
            %--------------------------------------------------------------------------
            %
            %   ------
            %   Inputs
            %   ------
            %   R               (1,3) Position vector
            %   V               (1,3) Velocity vector
            %   mu                    Gravitational parameter
            %   U               (1,1) Perurbation potential (i.e. negative disturbing
            %                         function). For U=0 classical "alternate"
            %                         equinoctial elements are treated.
            %
            %   -------
            %   Outputs
            %   -------
            %   X               (1,6) Elements vector [nu,p1,p2,varL,q1,Q2]
            %
            %--------------------------------------------------------------------------
            %   References:  Baù, G., Hernando-Ayuso, J. and Bombardelli, C., 2021.
            %   "A generalization of the equinoctial orbital elements"
            %   Celestial Mechanics and Dynamical Astronomy, 133(11), pp.1-29.
            %--------------------------------------------------------------------------


            %% Calculate GEqOE from R,V
            [~,~,INC,~,RAAN,~,~] = CoordConv.RV2OE(R,V,mu);  % obtain RAAN and INC
            % State parameters
            r = norm(R);
            v = norm(V);
            h = norm(cross(R,V));           % angular momentum module
            r_dot = dot(R,V)/r;             % radial velocity
            % Calculate energy
            epsK = 1/2*v^2 - mu/r;          % Keplerian energy

            phi = acos(R(3)/r);
            U = J2/2*mu/r*(Re/r)^2*(3*cos(phi)^2-1);

            eps = epsK + U;                 % total energy
            nu = 1/mu*(-2*eps)^(3/2);       % (1)
            % Calculate q in order to obtain equinoctial frame
            q1 = tan(INC/2)*sin(RAAN);  % (5)
            q2 = tan(INC/2)*cos(RAAN);  % (6)
            % Unit vectors of equinoctial frame in inertial ref.
            Kq = 1/(1+q1^2+q2^2);
            ex = Kq*[1-q1^2+q2^2    ,2*q1*q2        ,-2*q1];
            ey = Kq*[2*q1*q2        ,1+q1^2-q2^2    ,2*q2];
            % Radial vector in inertial ref.
            er = R/r;
            % Obtain true longitude
            cL = dot(er,ex);
            sL = dot(er,ey);
            L = atan2(sL,cL);

            if L < 0
                L = L + 2*pi;
            end


            if R(2) <0
                if ~(L >pi && L < 2*pi)
                    "here"
                end

            end

            % L = AOP + RAAN + TA;
            % Calculate rho,c and a
            Ueff = h^2/2/r^2 + U;   % effective potential
            c = sqrt(2*r^2*Ueff);
            rho = c^2/mu;
            a = -mu/2/eps;
            % g vector
            p1 = (rho/r-1)*sin(L) - c*r_dot/mu*cos(L);  % (2)
            p2 = (rho/r-1)*cos(L) + c*r_dot/mu*sin(L);  % (3)
            % Obtain K for varL
            w = sqrt(mu/a);
            sK = (mu+c*w-r*r_dot^2)*sin(L) - r_dot*(c+w*r)*cos(L);
            cK = (mu+c*w-r*r_dot^2)*cos(L) + r_dot*(c+w*r)*sin(L);
            K = atan2(sK,cK);
            if K < 0
                K = K + 2*pi;
            end
            % Obtain generalized mean longitud varL
            varL = K + 1/(mu+c*w)*(cK*p1-sK*p2);  % (4)
            % GEqOE array
            X = [nu p1 p2 varL q1 q2];
        end



        function [rr, vv] = po2pv(PO, mu)
            
            % input: PO vector of classical orbital parameters
            %        mu gravity parameter in units consistent with a
            % output: rr position vecotor in units consistent with mu and a
            %         vv velocity vector in units consistent with mu and a
            
            a = PO(1);
            e = PO(2);
            i = PO(3);
            Om = PO(4);
            om = PO(5);
            theta = PO(6);
            
            A = [cos(om+theta) -sin(om+theta) 0;
                sin(om+theta) cos(om+theta)  0;
                0             0   1];
            
            if i<0
                i = pi+i;
            end
            
            B = [1      0       0;
                0 cos(i)  -sin(i);
                0 sin(i)   cos(i)];
            
            C =  [cos(Om) -sin(Om) 0;
                sin(Om) cos(Om)  0;
                0       0   1];
            
            p = a*(1-e^2);
            
            r = [1/(1+e*cos(theta))*p 0 0]';
            v = sqrt(mu/p)*[e*sin(theta) 1+e*cos(theta) 0]';
            
            rr = C*B*A*r;
            vv = C*B*A*v;
        end
        function OPmat = KeplerStruct(EP)
            OP = CoordConv.ep2op(EP);
            OPmat.aAU = OP(1);
            OPmat.e = OP(2);
            OPmat.IncDeg = rad2deg((OP(3)));
            OPmat.AOPDeg = rad2deg((OP(4)));
            OPmat.RAANDeg = rad2deg((OP(5)));
            OPmat.TrueAnDeg = rad2deg((OP(6)));
            OPmat.ArgofLat = (OPmat.TrueAnDeg + OPmat.RAANDeg);
            
        end
        
        function OPmat = KeplerStruct2(EP)
            OP = CoordConv.mee2coe(EP);
            OPmat.aAU = OP(1);
            OPmat.e = OP(2);
            OPmat.IncDeg = rad2deg((OP(3)));
            OPmat.AOPDeg = rad2deg((OP(4)));
            OPmat.RAANDeg = rad2deg((OP(5)));
            OPmat.TrueAnDeg = rad2deg((OP(6)));
            OPmat.ArgofLat = (OPmat.TrueAnDeg + OPmat.RAANDeg);
            
        end
        
        function coe = mee2coe(mee)
            
            % convert modified equinoctial elements to classical orbit elements
            
            % input
            
            %  mee(1) = semiparameter (kilometers)
            %  mee(2) = f equinoctial element
            %  mee(3) = g equinoctial element
            %  mee(4) = h equinoctial element
            %  mee(5) = k equinoctial element
            %  mee(6) = true longitude (radians)
            
            % output
            
            %  coe(1) = semimajor axis (kilometers)
            %  coe(2) = eccentricity
            %  coe(3) = inclination (radians)
            %  coe(4) = right ascension of ascending node (radians)
            %  coe(5) = argument of periapsis (radians)
            %  coe(6) = true anomaly (radians)
            
            % Orbital Mechanics with MATLAB
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % unload modified equinoctial orbital elements
            
            pmee = mee(1);
            fmee = mee(2);
            gmee = mee(3);
            hmee = mee(4);
            kmee = mee(5);
            lmee = mee(6);
            
            % compute classical orbital elements
            
            tani2s = sqrt(hmee * hmee + kmee * kmee);
            
            % orbital eccentricity
            
            ecc = sqrt(fmee * fmee + gmee * gmee);
            
            % semimajor axis
            
            sma = pmee / (1.0 - ecc * ecc);
            
            % orbital inclination
            
            inc = 2.0 * atan(tani2s);
            
            % right ascension of ascending node
            
            raan = atan2(kmee, hmee);
            
            % argument of periapsis
            
            atopo = atan2(gmee, fmee);
            
            argper = mod(atopo - raan, 2.0 * pi);
            
            % true anomaly
            
            tanom = mod(lmee - atopo, 2.0 * pi);
            
            % load classical orbital element array
            
            coe(1) = sma;
            coe(2) = ecc;
            coe(3) = inc;
            coe(4) = raan;
            coe(5) = argper;
            coe(6) = tanom;
        end
        
        function  OP = ep2op(EP)
            % converts equinoctial parameters in orbital parameters
            % EP = (p,f,g,h,k,L)
            % OP = (a,e,i,Om,om,theta)
            
            % Initialize:
            p = EP(1);
            f = EP(2);
            g = EP(3);
            h = EP(4);
            k = EP(5);
            L = (EP(6));
            
            % % Compute:
            % OP(1) = p/(1-f^2-g^2);
            % OP(2) = sqrt(f^2+g^2);
            % OP(3) = 2*atan(sqrt(h^2+k^2));  % Problems finding right quadrant?
            % OP(4) = atan2(g,f)-atan2(k,h);  % Problems finding right quadrant?
            % OP(5) = atan2(k,h);             % Problems finding right quadrant?
            % OP(6) = L- atan2(g,f);          % Problems finding right quadrant?, check by atan(g/f)=O+w.
            
            OP(1) = p/(1-f^2-g^2);
            OP(2) = sqrt(f^2+g^2);
            OP(3) = atan2(2*sqrt(h^2+k^2), 1-h^2-k^2);
            
            if EP(4)==0&&EP(5)==0
                OP(4) = 0;
            else
                OP(4) = atan2(k,h);
            end
            if EP(2)==0&&EP(3)==0
                OP(5) = 0;
            else
                OP(5) = atan2(g*h -f*k,f*h+g*k);
            end
            OP(6) = L - OP(4) - OP(5);
        end
        function posandvel = ep2pv(EP, mu)
            OP = CoordConv.ep2op(EP);
            
            
            [rr, vv] = CoordConv.po2pv(OP, mu);
            posandvel = [rr;vv];
        end
        
        function x = vec2kepStruct(rs,vs,mus)
            
            OP = CoordConv.vec2orbElem(rs,vs,mus);
            
            x.a = OP(1);
            x.e = OP(2);
            x.IncDeg = rad2deg((OP(3)));
            x.AOPDeg = rad2deg((OP(5)));
            x.RAANDeg = rad2deg((OP(4)));
            x.TrueAnDeg = rad2deg((OP(6)));
            x.ArgofLat = (x.TrueAnDeg + x.RAANDeg);
            
        end
        
        function mee = vec2mee(rs,vs,mus)
            
            % Get the keplerian elements
            x = CoordConv.vec2orbElem(rs,vs,mus);
            
            a = x(1);
            e = x(2);
            I = x(3); 
            RAAN = x(4);
            AOP = x(5);
            True_an = x(6);
            
            % convert kepler to mee
            mee =  CoordConv.kepler2MEOE([a,e, RAAN,AOP, I, True_an]);
        end
        
        
        function x = vec2orbElem(rs,vs,mus)
            
            mus = mus(:).';
            nplanets = numel(rs)/3;
            if mod(nplanets,1) ~= 0 || numel(vs) ~= nplanets*3 ||...
                    (length(mus) ~= nplanets && length(mus) ~= 1)
                error('vec2orbElem:inputError',['rs and vs must contain 3n ',...
                    'elements and mus must have length n or 1 for n bodies']);
            end
            if length(rs) == numel(rs)
                rs = reshape(rs,3,nplanets);
            end
            if length(vs) == numel(vs)
                vs = reshape(vs,3,nplanets);
            end
            v2s = sum(vs.^2);
            r = sqrt(sum(rs.^2)); %orbital separation
            Ws = 0.5*v2s - mus./r;
            a = -mus/2./Ws; %semi-major axis
            L = [rs(2,:).*vs(3,:) - rs(3,:).*vs(2,:);...
                rs(3,:).*vs(1,:) - rs(1,:).*vs(3,:);...
                rs(1,:).*vs(2,:) - rs(2,:).*vs(1,:)]; %angular momentum
            L2s = sum(L.^2);
            p = L2s./mus; %semi-latus rectum
            e = sqrt(1 - p./a); %eccentricity
            
            e = abs(e);
            a = abs(a);
            %ecentric anomaly
            cosE = (1 - r./a)./e;
            
            sinE = sum(rs.*vs)./(e.*sqrt(mus.*a));
            E = atan2(sinE,cosE);
            
            %inclination
            sinI = sqrt(L(1,:).^2 + L(2,:).^2)./sqrt(L2s);
            cosI = L(3,:)./sqrt(L2s);
            I = atan2(sinI,cosI);
            
            sinw = ((vs(1,:).*L(2,:) - vs(2,:).*L(1,:))./mus - ...
                rs(3,:)./r)./(e.*sinI);
            cosw = ((sqrt(L2s).*vs(3,:))./mus - (L(1,:).*rs(2,:) - ...
                L(2,:).*rs(1,:))./(sqrt(L2s).*r))./(e.*sinI);
            RAAN = atan2(sinw,cosw);
         
            
            cosO = -L(2,:)./(sqrt(L2s).*sinI);
            sinO = L(1,:)./(sqrt(L2s).*sinI);
            AOP = atan2(sinO,cosO);
            
            %orbital periods
            P = 2*pi*sqrt(a.^3./mus);
            
            %time of periapsis crossing
            tau = -(E - e.*sin(E))./sqrt(mus.*a.^-3);
            
%               True_an = 2 * atan2(sqrt(1 + e) .* sin(E/2), sqrt(1 - e) .* cos(E/2));
%             Adjusted for desired output domain -pi:pi
            True_an = 2 * atan(sqrt((1 + e)/ (1 - e)) .* tan(E/2) );
            
                
            if dot(rs,vs) < 0

                True_an = wrapTo2Pi(True_an);

                if True_an > pi && True_an < 2*pi
                else
                    True_an = 2.0*pi - True_an;
                end

               
            end

            

            
            %time of periapsis crossing
            tau = -(E - e.*sin(E))./sqrt(mus.*a.^-3);
            
            x = [a,e,I,AOP, RAAN,True_an];
            
            
            %A and B vectors
            A = [a.*(cos(AOP).*cos(RAAN) - sin(AOP).*cos(I).*sin(RAAN));...
                a.*(sin(AOP).*cos(RAAN) + cos(AOP).*cos(I).*sin(RAAN));...
                a.*sin(I).*sin(RAAN)];
            B = [-a.*sqrt(1-e.^2).*(cos(AOP).*sin(RAAN) + ...
                sin(AOP).*cos(I).*cos(RAAN));...
                a.*sqrt(1-e.^2).*(-sin(AOP).*sin(RAAN) + ...
                cos(AOP).*cos(I).*cos(RAAN));...
                a.*sqrt(1-e.^2).*sin(I).*cos(RAAN)];
            
            
        end
    end
    
    
    
    
    
end
