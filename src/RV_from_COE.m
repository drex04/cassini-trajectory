  function [r,v] = RV_from_COE(a,e,i,RAAN,w,M,mu)
% A function to compute position and velocity from classical orbital
% elements
%
% Inputs:   a - semimajor axis, km
%           e - eccentricity
%           i - inclination, rad
%           RAAN - longitude of ascending node, rad
%           w - argument of periapsis, rad
%           M - mean anomaly, rad
%           mu - gravitational parameter, km^3/s^2
%
% Outputs:  r - 1x3 position vector in inertial frame, km
%           v - 1x3 velocity vector in inertial frame, km
%%
h = sqrt(mu*a*(1-e^2));
E = keplerE(e,M);
f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));

rp = (h^2/mu) * (1/(1 + e*cos(f))) * (cos(f)*[1;0;0] + sin(f)*[0;1;0]);
vp = (mu/h) * (-sin(f)*[1;0;0] + (e + cos(f))*[0;1;0]);

R3W = [ cos(RAAN)  sin(RAAN)  0; -sin(RAAN)  cos(RAAN)  0; 0        0     1];
R1i = [1       0          0; 0   cos(i)  sin(i); 0  -sin(i)  cos(i)];
R3w = [ cos(w)  sin(w)  0 ; -sin(w)  cos(w)  0; 0       0     1];

QpX = (R3w*R1i*R3W)';

r = QpX*rp;
v = QpX*vp;

end