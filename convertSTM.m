function STM_new = convertSTM(STMcart, xC, xD, fxC, fxD, mu, targetCoordSystem)
    switch lower(targetCoordSystem)

        case 'keplerian'
            J0 = jacobianCart2Kep(      xD,  mu);
            Jf = jacobianCart2Kep(      fxD, mu);

        case 'hillframe'
            J0 = jacobianCart2Hill(xC,  xD     );
            Jf = jacobianCart2Hill(fxC, fxD    );

        case 'roe'
            J0 =  jacobianCart2Roe(xC,  xD,  mu);
            Jf =  jacobianCart2Roe(fxC, fxD, mu);

        otherwise
            error('Unsupported coordinate system: %s', targetCoordSystem);
    end

    STM_new = Jf * STMcart * pinv(J0);
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

