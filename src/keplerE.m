function E = keplerE(e,M)
% Author: Drew Nollsch
% Date: 9/28/14
% ASTE 580, Professor Michael Gabor
%
% Inputs:   e - eccentricity
%           M - mean anomaly, rad
% Outputs:  E - eccentric anomaly, rad
%%
tol = 1.e-8;
% Set initial E
if M < pi
    E = M + e/2;
else
    E = M - e/2;
end

diff = 1;
while abs(diff) > tol
    diff = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - diff;
end

end