function F = keplerH(e, M)
% Author: Drew Nollsch
% Date: 9/28/14
% ASTE 580, Professor Michael Gabor
%
% Inputs:   e - eccentricity
%           M - hyperbolic mean anomaly, rad
% Outputs:  F - eccentric anomaly, rad
%%
tol = 1.e-8;

% Initialize F
F = M;

d = 1;
while abs(ratio) > tol
    ratio = (e*sinh(F) - F - M)/(e*cosh(F) - 1);
    F = F - ratio;
end

end