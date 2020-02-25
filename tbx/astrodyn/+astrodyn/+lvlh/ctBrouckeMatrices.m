function [A, B] = ctBrouckeMatrices(orbit, mc)
    e = orbit.Eccentricity;
    n = orbit.MeanMotion;
    v = @(t) orbit.TrueAnomaly(t);
    r = @(t) orbit.Radius(t);
    E = @(t) orbit.EccentricAnomaly(t);
    a = orbit.SemiMajorAxis;
    p = orbit.SemiLatusRectum;
    mu = orbit.CelestialBody.GravitationalParameter;
    
    C = sqrt( p * mu );
    thetad = @(t) C / r(t)^2;
    thetadd = @(t) -2  * mu * e * sin(v(t)) / r(t)^3;
    
    A = @(t) [ 0, 0, 0, 1, 0, 0; ...
        0, 0, 0, 0, 1, 0; ...
        0, 0, 0, 0, 0, 1; ...
        2*mu/r(t)^3 + thetad(t)^2, thetadd(t), 0, 0, 2*thetad(t), 0; ...
        -thetadd(t), -mu/r(t)^3 + thetad(t)^2, 0, -2*thetad(t), 0, 0; ...
        0, 0, -mu/r(t)^3, 0, 0, 0];
    
    B = [0, 0, 0; 0, 0, 0; 0, 0, 0; 1/mc, 0, 0; 0, 1/mc, 0; 0, 0, 1/mc];

end