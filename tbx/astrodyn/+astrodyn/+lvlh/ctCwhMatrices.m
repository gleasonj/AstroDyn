function [A, B] = ctCwhMatrices(x, y, z)
    if nargin < 3
        validateattributes(x, {'astrodyn.Orbit'}, {'nonempty'});
        validateattributes(y, {'numeric'}, {'scalar'});
    elseif nargin < 1
        error('Not enough input arguments')
    elseif nargin > 2
        error('Too many input arguments')
    else
        validateattributes(x, {'astrodyn.bodies.CelestialBody'}, {'nonempty'});
        validateattributes(y, {'numeric'}, {'scalar'});
        validateattributes(z, {'numeric'}, {'scalar'});
        if ~isprop(x, 'GravitationalParameter')
            error('Celestial body must have gravitational parameter defined');
        end
    end

    if isa(x, 'astrodyn.Orbit')
        n = sqrt(x.CelestialBody.GravitationalParameter / x.SemiMajorAxis^3);
        mc = y;
    else
        n = sqrt(x.GravitationalParameter / y^3);
        mc = z;
    end

    A = [0, 0, 0, 1, 0, 0; ...
        0, 0, 0, 0, 1, 0; ...
        0, 0, 0, 0, 0, 1; ...
        3*n^2, 0, 0, 0, 2*n, 0; ...
        0, 0, 0, -2*n, 0, 0;
        0, 0, -n^2, 0, 0, 0];

    B = [0, 0, 0; 0, 0, 0; 0, 0, 0; 1/mc, 0, 0; 0, 1/mc, 0; 0, 0, 1/mc];

end