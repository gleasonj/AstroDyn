classdef Orbit
    properties
        % Celestial body being orbited
        body

        % Eccentricity
        e

        % Semi-major axis [km]
        a

        % Inclination [degrees]
        i

        % Longitude of the ascending node [degrees]
        lon_an

        % Argument of the perapsis [degrees]
        w

        % Epoch
        M0

        % True anomaly at epoch
        v0

        % Orbit period [s]
        period
    end
    
    methods
        function obj = Orbit(body, e, a, inc, lon_an, w, M0, v0)
        % 
            obj.body = body;
            obj.e = e;
            obj.a = a;
            obj.i = inc;
            obj.lon_an = lon_an;
            obj.w = w;
            obj.M0 = M0;
            obj.v0 = v0;

            % oribtal period
            obj.period = (2 * pi - obj.e * sin(2 * pi)) / obj.n;
        end

        function val = n(obj)
        % Mean motion of orbit [1 / s]
            val = sqrt(obj.body.mu() / (obj.a^3 * 1000^3));
        end
        
        function val = p(obj)
        % The parameter / semi-latus rectum [km]
            if obj.e < 1
                val = obj.a * (1 - obj.e^2);
            else
                error(['Semi-latus rectum parameter not defined for ', ...
                    'eccentricities >= 1']);
            end
        end
        
        function [r, v] = getrv(obj, t)
            n = obj.n();
            E0 = obj.E0();
            a = obj.a;
            e = obj.e;
            c = E0 + e * sin(E0);

            E = 0;
            dE = 1;
            while dE > 1e-9
                Ep = n * (t - obj.M0) + e * sin(E) + c;
                dE = abs(E - Ep);
                E = Ep;
            end

            % sove for v (nu)
            v = arccos((e - cos(E)) / (e * cos(E) - 1))
            
            % solve for r
            r = a * (1 - e * cos(E));
        end

        function val = E0(obj)
        % Initial eccentric anomaly [degrees]
            val = acos((obj.e + cos(deg2rad(obj.v0))) / ...
                (1 + obj.e * cos(deg2rad(obj.v0))));
            val = rad2deg(val);
        end
        
        
        function x = toGeocentricFrame(obj)
        % 
            x = [];
        end

        function val = c(obj)
            val = sqrt(obj.e^2 * obj.a^2);
        end

        function val = b(obj)
            if obj.e <=1
                val = sqrt(obj.a^2 * (1 - obj.e^2));
            else
                error(['Ellipse parameter b does not exist for ', ...
                    'eccentricities > 1']);
            end
        end

        function plot(obj)
        % Plot the orbit
            v = linspace(0, 2 * pi, 1000);
            b = obj.b();
            p = obj.p();

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
end