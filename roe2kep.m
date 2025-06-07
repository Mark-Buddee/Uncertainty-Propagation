function kepdeputy = roe2kep(kepcheif, doe)

%ud = kepdeputy(5) + true2meanAnomaly(kepdeputy(6), kepdeputy(2));
uc = kepcheif(5) + true2meanAnomaly(kepcheif(6), kepcheif(2));


kepdeputy(1)= doe(1)*kepcheif(1) + kepcheif(1); 
kepdeputy(3) = doe(5) + kepcheif(3); 
kepdeputy(4) = doe(6)/sin(kepcheif(3)) + kepcheif(4); 

ud = doe(2) - (kepdeputy(4) - kepcheif(4))*cos(kepcheif(3))+uc; 
edcoswd = doe(3) + kepcheif(2)*cos(kepcheif(5));
edsinwd = doe(4) + kepcheif(2)*sin(kepcheif(5));

kepdeputy(5) = atan2(edsinwd,edcoswd);
kepdeputy(2) = sqrt(edcoswd^2 + edsinwd^2);
kepdeputy(6) = mean2trueAnomaly(ud - kepdeputy(5),kepdeputy(2));

