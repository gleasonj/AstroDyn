classdef Orbit
    properties (Dependent)
        CelestialBody

        Eccentricity

        SemiMajorAxis

        Inclination

        LongitudeOfAscendingNode

        ArgumentOfPeriapsis

        Epoch

        TrueAnomalyAtEpoch

        EccentricAnomalyAtEpoch

        OrbitPeriod

        SemiLatusRectum

        MeanMotion
    end

    properties (Dependent, Access = private)
        lan_

        mu_
    end

    properties (Access = protected)
        % Celestial body being orbited
        CelestialBody_

        % Eccentricity
        Eccentricity_

        % Semi-major axis [km]
        SemiMajorAxis_

        % Inclination [degrees]
        Inclination_

        % Longitude of the ascending node [degrees]
        LongitudeOfAscendingNode_

        % Argument of the perapsis [degrees]
        ArgumentOfPeriapsis_

        % Epoch
        Epoch_

        % True anomaly at epoch
        TrueAnomalyAtEpoch_
    end
    
    methods
        function obj = Orbit(body, e, a, inc, lon_an, w, M0, v0)
        % 
            obj.CelestialBody_       = body;
            obj.Eccentricity_        = e;
            obj.SemiMajorAxis_       = a;
            obj.Inclination_         = inc;
            obj.lan_                 = lon_an;
            obj.ArgumentOfPeriapsis_ = w;
            obj.Epoch_               = M0;
            obj.TrueAnomalyAtEpoch_  = v0;
        end

        function val = get.Eccentricity(obj)
            val = obj.Eccentricity_;
        end

        function val = get.CelestialBody(obj)
            val = obj.CelestialBody_;
        end

        function val = get.SemiMajorAxis(obj)
            val = obj.SemiMajorAxis_;
        end

        function val = get.LongitudeOfAscendingNode(obj)
            val = obj.LongitudeOfAscendingNode_;
        end

        function val = get.Inclination(obj)
            val = obj.Inclination_;
        end

        function val = get.ArgumentOfPeriapsis(obj)
            val = obj.ArgumentOfPeriapsis_;
        end

        function val = get.Epoch(obj)
            val = obj.Epoch_;
        end

        function val = get.TrueAnomalyAtEpoch(obj)
            val = obj.TrueAnomalyAtEpoch_;
        end

        function val = get.OrbitPeriod(obj)
        % oribtal period
            val = (2 * pi - obj.Eccentricity_ * sin(2 * pi)) / obj.MeanMotion;
        end

        function val = get.MeanMotion(obj)
        % Mean motion of orbit [1 / s]
            mu = obj.mu_;
            a  = obj.SemiMajorAxis;

            val = sqrt(mu / (a^3 * 1000^3));
        end

        function val = get.SemiLatusRectum(obj)
        % The parameter / semi-latus rectum [km]
            if obj.Eccentricity < 1
                val = obj.SemiMajorAxis * (1 - obj.Eccentricity^2);
            else
                error(['Semi-latus rectum parameter not defined for ', ...
                    'eccentricities >= 1']);
            end
        end

        function r = Radius(obj, t, varargin)
            [r, ~] = obj.radius_trueanomaly(obj, t, varagin{:})
        end

        function v = TrueAnomaly(obj, t, varagin)
            [~, v] = obj.radius_trueanomaly(obj, t, varagin{:})
        end

        function [r, v] = radiusAndTrueAnomaly(obj, t, varagin)
            a = obj.SemiMajorAxis;
            e = obj.Eccentricity;
            
            E = obj.EccentricAnomaly(t, varagin{:});

            % sove for v (nu)
            v = arccos((e - cos(E)) / (e * cos(E) - 1));
            r = a * (1 - e * cos(E));
        end

        function E = EccentricAnomaly(obj, t, varargin)
            p = inputParser();
            addRequired(p, 't', @(x) validateattributes(x, {'numeric'}, ...
                {'scalar'}));
            addParameter(p, 'Accuracy', 1e-9, ...
                @(x) validateattributes(x, {'numeric'}, ...
                    {'scalar', 'positive'}));

            parse(p, t, varargin{:});

            n  = obj.MeanMotion;
            E0 = obj.EccentricAnomalyAtEpoch;
            M0 = obj.Epoch;
            e  = obj.Eccentricity;

            c = E0 + e * sin(E0);

            E = 0;
            dE = 1;
            while dE > p.Results.Accuracy
                Ep = n * (t - obj.M0) + e * sin(E) + c;
                dE = abs(E - Ep);
                E = Ep;
            end
        end

        function val = get.EccentricAnomalyAtEpoch(obj)
            val = self.E0_();
        end
        
        
        function x = toGeocentricFrame(obj)
        % 
            x = [];
        end

        function plot(obj)
        % Plot the orbit
            v = linspace(0, 2 * pi, 1000);
            b = obj.b_();
            p = obj.SemiLatusRectum;

            Rinv = inv(obj.R());

            x = zeros(size(v));
            y = zeros(size(v));
            z = zeros(size(v));
            for i = 1:length(v)
                r = p / (1 + obj.e * cos(v(i)));


                temp = r * cos(v(i)) * [1, 0, 0]' + r * sin(v(i)) * [0, 1, 0]';
                temp = Rinv * temp;
                x(i) = temp(1);
                y(i) = temp(2);
                z(i) = temp(3);
            end

            plot(x, y)
            % plot3(x, y, z)
        end

        function val = R(obj)
        % Rotation matrix from geocentric frame to perifocal frame

            lan = deg2rad(obj.lon_an);
            w = deg2rad(obj.w);
            inc = deg2rad(obj.i);

            val = [
                cos(lan) * cos(w) - sin(lan) * sin(w) * cos(inc), ...
                -cos(lan) * sin(w) - sin(lan) * cos(w) * cos(inc), ...
                sin(lan) * sin(inc); ...

                sin(lan) * cos(w) + cos(lan) * sin(w) * cos(inc), ...
                -sin(lan) * sin(w) + cos(lan) * cos(w) * cos(inc), ...
                -cos(lan) * sin(inc); ...


                sin(w) * sin(inc), ...
                cos(w) * sin(inc), ...
                cos(inc)
            ];
        end
    end

    methods (Access = private)
        function val = E0_(obj)
        % Initial eccentric anomaly [degrees]
            e  = obj.Eccentricity_;
            v0 = deg2rad(obj.TrueAnomalyAtEpoch_);

            val = acos((e + cos(v0)) / (1 + e * cos(v0)));
            val = rad2deg(val);
        end

        function val = get.lan_(obj)
        % Shorthand for getting the longitude of the Ascending Node
            val = obj.LongitudeOfAscendingNode;
        end

        function set.lan_(obj, val)
            obj.LongitudeOfAscendingNode_ = val;
        end

        function val = c_(obj)
            val = sqrt(obj.e^2 * obj.a^2);
        end

        function val = get.mu_(obj)
            val = obj.CelestialBody.GravitationalParameter;
        end

        function val = b_(obj)
            if obj.e <=1
                val = sqrt(obj.a^2 * (1 - obj.e^2));
            else
                error(['Ellipse parameter b does not exist for ', ...
                    'eccentricities > 1']);
            end
        end
    end
end