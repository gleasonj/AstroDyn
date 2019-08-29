classdef Earth < astrodyn.bodies.CelestialBody
    properties (SetAccess = immutable)
        % Axial tilt [degrees]
        axial_tilt
    end

    methods
        function obj = Earth()
            obj = obj@astrodyn.bodies.CelestialBody('Mass', 5.97237e24, ...
                'Radius', 6378.1);

            obj.axial_tilt = 23.4392811;
        end
    end
end