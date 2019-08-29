classdef CelestialBody
    properties
        % Mass [kg]
        mass

        % Equitorial radius [km]
        radius
    end

    methods
        function obj = CelestialBody(varargin)
            p = inputParser();

            addParameter(p, 'Mass', [], ...
                @(x) validateattributes(x, {'numeric'}, {'scalar'}));

            addParameter(p, 'Radius', [], ...
                @(x) validateattributes(x, {'numeric'}, {'scalar'}));

            parse(p, varargin{:});

            obj.mass   = p.Results.Mass;
            obj.radius = p.Results.Radius;
        end

        function val = mu(obj)
        % Gravitational parameter of celestial body [m^3 / s^2]
            if isempty(obj.mass)
                error(['Celestial body has no mass: cannot compute ', ...
                    'gravitational parameter']);
            else
                val = astrodyn.constants.G() * obj.mass;
            end
        end
    end
end