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

        % Semi-major axis [m]
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

            val = sqrt(mu / a^3);
        end

        function val = get.SemiLatusRectum(obj)
        % The parameter / semi-latus rectum [m]
            if obj.Eccentricity < 1
                val = obj.SemiMajorAxis * (1 - obj.Eccentricity^2);
            else
                error(['Semi-latus rectum parameter not defined for ', ...
                    'eccentricities >= 1']);
            end
        end

        function r = Radius(obj, t, varargin)
            [r, ~] = obj.radiusAndTrueAnomaly(t, varargin{:});
        end

        function v = TrueAnomaly(obj, t, varargin)
            [~, v] = obj.radiusAndTrueAnomaly(t, varargin{:});
        end

        function [r, v] = radiusAndTrueAnomaly(obj, t, varargin)
            a = obj.SemiMajorAxis;
            e = obj.Eccentricity;
            
            E = obj.EccentricAnomaly(t, varargin{:});

            % sove for v (nu)
            v = acos((e - cos(E)) / (e * cos(E) - 1));

            if E > pi
                v = v + 2 * (pi - v);
            end

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
                Ep = n * (t - M0) + e * sin(E) + c;
                dE = abs(E - Ep);
                E = Ep;
            end

            E = mod(E, 2*pi);
        end

        function val = get.EccentricAnomalyAtEpoch(obj)
            val = obj.E0_();
        end
        
        
        function x = toGeocentricFrame(obj)
        % 
            x = [];
        end

        function varargout = plot(obj, varargin)
        % Plot the orbit
            inp = inputParser();
            inp.KeepUnmatched = true;
            addParameter(inp, 'Axis', [], ...
                @(x) isa(x, 'matlab.graphics.axis.Axes'));
            addParameter(inp, 'Plot3d', false, @(x) isa(x, 'logical'));
            addParameter(inp, 'NumberOfPoints', 1000, ...
                @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
            
            parse(inp, varargin{:});

            v = linspace(0, 2 * pi, inp.Results.NumberOfPoints);
            b = obj.b_();
            p = obj.SemiLatusRectum;

            Rinv = inv(obj.R());

            x = zeros(size(v));
            y = zeros(size(v));
            z = zeros(size(v));
            for i = 1:length(v)
                r = p / (1 + obj.Eccentricity * cos(v(i)));


                temp = r * cos(v(i)) * [1, 0, 0]' + r * sin(v(i)) * [0, 1, 0]';
                temp = Rinv * temp;
                x(i) = temp(1);
                y(i) = temp(2);
                z(i) = temp(3);
            end

            if inp.Results.Plot3d
                hl = plot3(x, y, z);
            else
                hl = plot(x, y);
            end
            
            if nargout == 1
                varargout = {hl};
            elseif nargout > 1
                error('Too many output arguments.');
            end
        end

        function varargout = plot3(obj, varargin)
            hl = plot(obj, 'Plot3d', true, varargin{:});

            if nargout == 1
                varargout = {hl};
            elseif nargout > 1
                error('Too many output arguments.');
            end
        end

        function val = R(obj)
        % Rotation matrix from geocentric frame to perifocal frame

            lan = deg2rad(obj.lan_);
            w = deg2rad(obj.ArgumentOfPeriapsis);
            inc = deg2rad(obj.Inclination);

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

        function val = get.lan_(obj)
        % Shorthand for getting the longitude of the Ascending Node
            val = obj.LongitudeOfAscendingNode_;
        end

        function obj = set.lan_(obj, val)
            obj.LongitudeOfAscendingNode_ = val;
        end

        function val = get.mu_(obj)
            val = obj.CelestialBody.GravitationalParameter;
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

        function val = c_(obj)
            val = sqrt(obj.e^2 * obj.a^2);
        end

        function val = b_(obj)
            if obj.Eccentricity <=1
                val = sqrt(obj.SemiMajorAxis^2 * (1 - obj.Eccentricity^2));
            else
                error(['Ellipse parameter b does not exist for ', ...
                    'eccentricities > 1']);
            end
        end
    end

    methods (Static)
        function orbit = fromCelestialCeneterdInertialState(body, x)
            % Get orbit from the Celestial-Body-Cenetered Inertial state vector
            % 
            % State vector must be in m and m/sec
            if size(x, 1) == 1
                x = x';
            end

            r = x(1:3);
            v = x(4:6);

            mu = body.GravitationalParameter / 1000^3;

            h = cross(r, v);
            n = cross([0, 0, 1]', h);

            p = norm(h)^2 / mu;
            e = (1 / mu) * ((norm(v)^2 - mu/norm(r)^2)*r - (r'*v)*v);
            a = p / (1 - norm(e)^2);

            i = acos(h(3) / norm(h));
            lan = acos(n(1) / norm(n));
            w = acos(n' * e / (norm(n) * norm(e)));
            v0 = acos(e' * r / (norm(e) * norm(r)));

            orbit = astrodyn.Orbit(body, norm(e), a, i, lan, w, 0, v0);
        end

        function orbit = fromEciState(x)
            % Get orbit from the Earth-Cenetered Inertial state vector
            body = astrodyn.bodies.Earth();
            orbit = astrodyn.Orbit.fromCelestialCeneterdInertialState(body, x);
        end
    end
end