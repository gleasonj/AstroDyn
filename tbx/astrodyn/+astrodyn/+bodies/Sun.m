classdef Earth < astrodyn.bodies.CelestialBody
    properties (Dependent)
        Mass
        EquitorialRadius
        Radius
        GravitationalParameter
    end

    properties (Access = private)
        Mass_       = 1.9885e30

        % Galactic period [years]
        GalacticPeriod_ = 2.25e8

        EquitorialRadius_ = 695700

        EquitorialCircumference_ = 4.379e6
    end

    methods
        function val = get.Mass(obj)
            val = obj.Mass_;
        end

        function val = get.Radius(obj)
            val = obj.EquitorialRadius;
        end

        function val = get.EquitorialRadius(obj)
            val = obj.EquitorialRadius_;
        end

        function val = get.GravitationalParameter(obj)
            val = astrodyn.constants.G() * obj.Mass;
        end
    end

    methods (Access = private)
        function val = mu_(obj)
        % Shorthand for getting the gravitational parameter
            val = obj.GravitationalParameter;
        end
    end
end