function doe = kep2roe(kepcheif, kepdeputy)

ud = kepdeputy(5) + true2meanAnomaly(kepdeputy(6), kepdeputy(2));
uc = kepcheif(5) + true2meanAnomaly(kepcheif(6), kepcheif(2));

doe(1) = (kepdeputy(1) - kepcheif(1)) / kepcheif(1);

% Wrap angle difference to (-pi, pi]
delta_lambda = wrapToPi((ud - uc) + (kepdeputy(4) - kepcheif(4)) * cos(kepcheif(3)));
doe(2) = delta_lambda;

doe(3) = kepdeputy(2) * cos(kepdeputy(5)) - kepcheif(2) * cos(kepcheif(5));
doe(4) = kepdeputy(2) * sin(kepdeputy(5)) - kepcheif(2) * sin(kepcheif(5));
doe(5) = kepdeputy(3) - kepcheif(3);
doe(6) = (kepdeputy(4) - kepcheif(4)) * sin(kepcheif(3));