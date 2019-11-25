classdef Earth < astrodyn.bodies.CelestialBody
    properties (Dependent)
        AxialTilt
        Mass
        EquitorialRadius
        Radius
        GravitationalParameter
    end

    properties (Access = private)
        AxialTilt_  = 23.4392811
        Mass_       = 5.97237e24
        Aphelion_   = 1.521e8
        Perihelion_ = 1.47095e8
        SemiMajorAxis_ = 1.49598023e8
        Eccentricity_ = 0.0167086
        OrbitalPeriod_ = 365.256363004 * 86400

        % Average orbital speed [km/s]
        AverageOrbitalSpeed_ = 29.78 

        LongitudeAscendingNode_ = -11.2664
        ArgumentOfPeriapsis_ = 114.20783
        MeanRadius_ = 6371
        EquitorialRadius_ = 6378.1
        PolarRadius_ = 6356.8
        EquitorialCircumference_ = 40075.017
        MeridonalCircumference_  = 40007.86

        Inclination_ = 7.155
    end

    methods
        function val = get.AxialTilt(obj)
            val = obj.AxialTilt_;
        end

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

        function orb = getOrbit(obj)
            orb = astrodyn.Orbit(astrodyn.bodies.Sun(), obj.Eccentricity_, ...
                obj.SemiMajorAxis_, obj.Inclination_, ...
                obj.LongitudeAscendingNode_, obj.ArgumentOfPeriapsis_, ...
                astrodyn.constants.J2000(), 0)
        end
    end

    methods (Access = private)
        function val = mu_(obj)
        % Shorthand for getting the gravitational parameter
            val = obj.GravitationalParameter;
        end
    end
end