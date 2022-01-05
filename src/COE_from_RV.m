function [a,e,i,RAAN,w,n,f,T] = COE_from_RV(r,v,mu)
% A function to compute classical orbital elements from position and
% velocity vectors
%
% Inputs:   r - 1x3 position vector in inertial frame, km
%           v - 1x3 velocity vector in inertial frame, km
%           mu - gravitational parameter, km^3/s^2
%
% Outputs:  a - semimajor axis, km
%           e - eccentricity
%           i - inclination, rad
%           RAAN - longitude of ascending node, rad
%           w - argument of periapsis, rad
%           n - mean motion, rad/sec
%           f - true anomaly, rad
%           T - orbital period, sec
%% 
%eps = 0;
eps = 1e-6;

rmag    = norm(r);
vmag    = norm(v);

vr   = dot(r,v)/rmag;

h    = cross(r,v);
hmag    = norm(h);

%...Equation 4.7:
i = acos(h(3)/hmag);

%...Equation 4.8:
n    = cross([0 0 1],h);
nmag    = norm(n);

%...Equation 4.9 (incorporating the case incl = 0):
if i ~= 0 %Inclined orbit
    RAAN = acos(n(1)/nmag);
    if n(2) < 0
        RAAN = 2*pi - RAAN;
    end
else %Equatorial orbit
    RAAN = 0;
end

%...Equation 4.10:
E = 1/mu*((vmag^2 - mu/rmag)*r - rmag*vr*v);
e = norm(E);

%...Equation 4.12 (incorporating the cases incl = 0 and e = 0):
if i ~= 0   %Inclined orbit
    if e > eps %Non-circular orbit
        w = acos(dot(n,E)/nmag/e);
        if E(3) < 0
            w = 2*pi - w;
        end
    else       %Circular orbit
        w = 0;
    end
else %Equatorial orbit
    if e > eps %Non-circular orbit
        w = acos(E(1)/e);
        if E(2) <0
            w = 2*pi - w;
        end
    else       %Circular orbit
        w = 0;
    end
end

%...Equation 4.13a (incorporating the cases incl = 0 and e = 0):
if i ~= 0   %Inclined orbit
    if e > eps %Non-circular orbit
        f = acos(dot(E,r)/e/rmag);
        if vr < 0
            f = 2*pi - f;
        end
    else       %Circular orbit
        f = acos(dot(n,r)/nmag/rmag);
        if r(3) < 0
            f = 2*pi - f;
        end
    end
else           %Equatorial orbit
    if e > eps %Non-circular orbit
        f = acos(dot(E,r)/e/rmag);
        if vr < 0
            f = 2*pi - f;
        end
    else       %Circular orbit
        f = acos(r(1)/rmag);
        if r(2) < 0
            f = 2*pi - f;
        end
    end
end

%...Equation 4.62 (a < 0 for a hyperbola):
a = hmag^2/mu/(1 - e^2);
n = sqrt(mu/a^3);
T = 2*pi/n;
end