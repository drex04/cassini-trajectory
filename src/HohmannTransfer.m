function [dv_total,dv1,dv2,ToF] = HohmannTransfer(r1,r2,mu)
% Author: Drew Nollsch
% Date: 11/23/14
% ASTE 580, Professor Michael Gabor
%
% A function to calculate Hohmann transfer velocity changes and time of
% flight from initial and final radii
%
% Inputs:   r1, magnitude of initial radius, km
%           r2, magnitude of final radius, km
%           mu, gravitational parameter, km^3/s^2
%
% Outputs:  dv_total, magnitude of total velocity change, km/s
%           dv1, magnitude of departure velocity change, km/s
%           dv2, magnitude of arrival velocity change, km/s
%           ToF, time of flight, days
%
%%
a1 = r1; % km
a2 = r2; % km

% Conditions before maneuver 1
v1 = sqrt(2*mu/r1 - mu/a1); % km/s
gamma1 = 0; % rad

% Determine transfer ellipse
rpt = r1; % km
rat = r2; % km
at = (rpt+rat)/2; % km
et = -(rpt/at-1);

% Determine conditions after maneuver 1
vpt = sqrt(2*(-mu/(2*at)+mu/rpt));
gammapt = 0;

% Compute dv1
dv1 = vpt-v1; % km/s

% Determine conditions prior to maneuver 2
vat = sqrt(2*(-mu/(2*at)+mu/rat)); % km/s
gammaat = 0;

% Determine conditions after maneuver 2
v2 = sqrt(mu/r2); % km/s

% Compute dv2
dv2 = v2-vat; % km/s

% Total dv
dv_total = dv1+dv2; % km/s

% Time of Flight
ToF = pi*sqrt(at^3/mu); % seconds
ToF = ToF/60/60/24; % days

end
