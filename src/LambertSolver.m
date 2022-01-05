function [v1,v2] = LambertSolver(r1,r2,ToF,mu,string)
% Author: Drew Nollsch
% Date: 11/9/14
% ASTE 580, Professor Michael Gabor
%
% Inputs:   r1, initial position vector, km
%           r2, final position vector, km
%           ToF, time of flight from r1 to r2, sec
%           string, 'pro' if prograde, 'retro' if retrograde
%           mu, gravitational parameter of central body, km^3/s^2
%
% Outputs:  v1, initial velocity vector, km/s
%           v2, final velocity vector, km/s
%%
r1mag = norm(r1);
r2mag = norm(r2);

cross12   = cross(r1, r2);
theta = acos(dot(r1,r2)/r1mag/r2mag);

if nargin < 5 || (~strcmp(string,'retro') & (~strcmp(string,'pro')))
    string = 'pro';
    fprintf('Assume prograde orbit\n')
end

if strcmp(string,'pro')
    if cross12(3) <= 0
        theta = 2*pi - theta;
    end
elseif strcmp(string,'retro')
    if cross12(3) >= 0
        theta = 2*pi - theta;
    end
end

% Equation 5.35:
A = sin(theta)*sqrt(r1mag*r2mag/(1 - cos(theta)));

% Determine approximately where F(z,t) changes sign, and use that value of z as the starting value for Equation 5.45:
z = -100;
while F(z,ToF) < 0
    z = z + 0.1;
end

tol = 1.e-8;
imax = 1e6;

% Iterate on Equation 5.45 until tolerance is met
ratio = 1;
i = 0;
while (abs(ratio) > tol) & (i <= imax)
    i     = i + 1;
    ratio = F(z,ToF)/dFdz(z);
    z     = z - ratio;
end

% Report if the maximum number of iterations is exceeded:
if i >= imax
    fprintf('\n\nNumber of iterations exceeds %g in LambertSolver\n\n',imax)
end

%...Equation 5.46a:
f = 1 - y(z)/r1mag;

%...Equation 5.46b:
g  = A*sqrt(y(z)/mu);

%...Equation 5.46d:
gdot = 1 - y(z)/r2mag;

%...Equation 5.28:
v1 = 1/g*(r2 - f*r1);

%...Equation 5.29:
v2 = 1/g*(gdot*r2 - r1);

return


% Equation 5.38:
    function dum = y(z)
        dum = r1mag + r2mag + A*(z*S(z) - 1)/sqrt(C(z));
    end

% Equation 5.40:
    function dum = F(z,t)
        dum = (y(z)/C(z))^1.5*S(z) + A*sqrt(y(z)) - sqrt(mu)*t;
    end

% Equation 5.43:
    function dum = dFdz(z)
        if z == 0
            dum = sqrt(2)/40*y(0)^1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));
        else
            dum = (y(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) ...
                + 3*S(z)^2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z)) ...
                + A*sqrt(C(z)/y(z)));
        end
    end

% Stumpff functions:
    function dum = C(z)
        [C,S] = CS(z);
        dum = C;
    end

    function dum = S(z)
        [C,S] = CS(z);
        dum = S;
    end
end
