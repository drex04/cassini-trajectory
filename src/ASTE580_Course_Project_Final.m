% Drew Nollsch
% ASTE 580 Course Project
% Simulated Cassini Trajectory
% 11/18/14
clc
clear all

%% Planetary Orbital Elements

% Sun
mu_sun = 1.32712e11; % gravitational parameter, km^3/s^2
% Set as origin of heliocentric frame
r_sun = [0;0;0];
v_sun = [0;0;0];

% Venus
a_venus = 1.08209e8; % semimajor axis, km
e_venus = 0.00677675; % eccentricity
i_venus = 3.39477; % inclination to the ecliptic, deg
RAAN_venus = 76.7004; % right ascension of ascending node, deg
w_venus = 54.8347; % argument of periapsis, deg
M_venus = 196.094; % mean anomaly, deg
mu_venus = 3.249e5; % gravitational parameter, km^3/s^2
Rv = 6052; % mean radius, km (Note:  planetary radii are from http://nineplanets.org/data1.html)

% Earth
a_earth = 1.49598e8;
e_earth = 0.0167129;
i_earth = 0.000980125;
RAAN_earth = 349.101;
w_earth = 113.823;
M_earth = 309.536;
mu_earth = 3.986e5;
Re = 6378;

% Jupiter
a_jupiter = 7.78406e8;
e_jupiter = 0.0484018;
i_jupiter = 1.30538;
RAAN_jupiter = 100.532;
w_jupiter = 274.205;
M_jupiter = 163.185;
mu_jupiter = 1.26686e8;
Rj = 71492;

% Saturn
a_saturn = 1.42676e9;
e_saturn = 0.0541768;
i_saturn = 2.48434;
RAAN_saturn = 113.747;
w_saturn = 338.724;
M_saturn = 230.266;
mu_saturn = 3.7931e7;
Rs = 60268;

% Convert all angular values from degrees to radians
i_venus = i_venus*2*pi/360;
i_earth = i_earth*2*pi/360;
i_jupiter = i_jupiter*2*pi/360;
i_saturn = i_saturn*2*pi/360;
RAAN_venus = RAAN_venus*2*pi/360;
RAAN_earth = RAAN_earth*2*pi/360;
RAAN_jupiter = RAAN_jupiter*2*pi/360;
RAAN_saturn = RAAN_saturn*2*pi/360;
w_venus = w_venus*2*pi/360;
w_earth = w_earth*2*pi/360;
w_jupiter = w_jupiter*2*pi/360;
w_saturn = w_saturn*2*pi/360;
M_venus = M_venus*2*pi/360;
M_earth = M_earth*2*pi/360;
M_jupiter = M_jupiter*2*pi/360;
M_saturn = M_saturn*2*pi/360;

% Calculate mean motion of the planets
n_earth = sqrt(mu_sun/a_earth^3);
n_venus = sqrt(mu_sun/a_venus^3);
n_jupiter = sqrt(mu_sun/a_jupiter^3);
n_saturn = sqrt(mu_sun/a_saturn^3);

%% Mission Timing

% Convert to Julian Date
% Epoch for Planetary Orbital Elements
[POEepoch] = DateConvert(1992,11,13,0,0,0);
% Earth departure
[Earth1JD] = DateConvert(1997,10,15,12,0,0);
% First Venus flyby
[Venus1JD] = DateConvert(1998,4,26,12,0,0);
% Second Venus flyby
[Venus2JD] = DateConvert(1999,6,24,12,0,0);
% Earth flyby
[Earth2JD] = DateConvert(1999,8,18,12,0,0);
% Jupiter flyby
[JupiterJD] = DateConvert(2000,12,30,12,0,0);
% Saturn arrival
[SaturnJD] = DateConvert(2004,7,1,12,0,0);

% Calculate time of flight for each mission segment
TOF_E1toV1 = Venus1JD - Earth1JD;
TOF_V1toV2 = Venus2JD - Venus1JD;
TOF_V2toE2 = Earth2JD - Venus2JD;
TOF_E2toJ = JupiterJD - Earth2JD;
TOF_JtoS = SaturnJD - JupiterJD;
% Convert time of flight from days to seconds
TOF_E1toV1 = TOF_E1toV1*24*60*60;
TOF_V1toV2 = TOF_V1toV2*24*60*60;
TOF_V2toE2 = TOF_V2toE2*24*60*60;
TOF_E2toJ = TOF_E2toJ*24*60*60;
TOF_JtoS = TOF_JtoS*24*60*60;

% Initialize time vector from initial epoch to Saturn arrival
t = POEepoch:0.5:SaturnJD;
% Find the index of the time element for Earth departure
timetestE1 = t==Earth1JD;
indexE1 = find(timetestE1);
% Concatenate t to start at Earth departure
t = t(indexE1:end);
% Reset the index of Earth departure to match concatenated t array
timetestE1 = t==Earth1JD;
indexE1 = find(timetestE1);
% Find the indexes of the other rendezvous
timetestV1 = t==Venus1JD;
timetestV2 = t==Venus2JD;
timetestE2 = t==Earth2JD;
timetestJ = t==JupiterJD;
timetestS = t==SaturnJD;
indexV1 = find(timetestV1);
indexV2 = find(timetestV2);
indexE2 = find(timetestE2);
indexJ = find(timetestJ);
indexS = find(timetestS);

%% Calculate mean anomalies, positions, and velocities of the bodies for the entire mission
% Note:  Dates of rendezvous points as well as the epoch are in days, so
% have to convert to seconds
for i = 1:length(t)
    % Earth
    M_E(i) = M_earth + n_earth*(t(i)-POEepoch)*24*60*60;
    [r_E(i,:),v_E(i,:)] = RV_from_COE(a_earth,e_earth,i_earth,RAAN_earth,w_earth,M_E(i),mu_sun);
    % Venus
    M_V(i) = M_venus + n_venus*(t(i)-POEepoch)*24*60*60;
    [r_V(i,:),v_V(i,:)] = RV_from_COE(a_venus,e_venus,i_venus,RAAN_venus,w_venus,M_V(i),mu_venus);
    % Jupiter
    M_J(i) = M_jupiter + n_jupiter*(t(i)-POEepoch)*24*60*60;
    [r_J(i,:),v_J(i,:)] = RV_from_COE(a_jupiter,e_jupiter,i_jupiter,RAAN_jupiter,w_jupiter,M_J(i),mu_jupiter);
    % Saturn
    M_S(i) = M_saturn + n_saturn*(t(i)-POEepoch)*24*60*60;
    [r_S(i,:),v_S(i,:)] = RV_from_COE(a_saturn,e_saturn,i_saturn,RAAN_saturn,w_saturn,M_S(i),mu_saturn);
end

%% Calculate transfer angle for each mission segment
phi_E1toV1 = 2*pi - acos(dot(r_E(indexE1,:),r_V(indexV1,:))/(norm(r_E(indexE1,:))*norm(r_V(indexV1,:))));
phi_V1toV2 = acos(dot(r_V(indexV1,:),r_V(indexV2,:))/(norm(r_V(indexV1,:))*norm(r_V(indexV2,:))));
phi_V2toE2 = acos(dot(r_V(indexV2,:),r_E(indexE2,:))/(norm(r_V(indexV2,:))*norm(r_E(indexE2,:))));
phi_E2toJ = acos(dot(r_E(indexE2,:),r_J(indexJ,:))/(norm(r_E(indexE2,:))*norm(r_J(indexJ,:))));
phi_JtoS = acos(dot(r_J(indexJ,:),r_S(indexJ,:))/(norm(r_J(indexJ,:))*norm(r_S(indexS,:))));

%% Earth Departure to Venus Flyby (E1toV1)
% Use Lambert function to find transfer orbit and departure and arrival
% velocities
[v1_E1toV1,v2_E1toV1] = lambert(r_E(indexE1,:),r_V(indexV1,:),TOF_E1toV1,'pro',mu_sun);
% Find COEs of transfer orbit
[a_E1toV1,e_E1toV1,i_E1toV1,RAAN_E1toV1,w_E1toV1,n_E1toV1,fdep_E1toV1,T_E1toV1] = COE_from_RV(r_E(indexE1,:),v1_E1toV1,mu_sun);
[a_E1toV1,e_E1toV1,i_E1toV1,RAAN_E1toV1,w_E1toV1,n_E1toV1,farr_E1toV1,T_E1toV1] = COE_from_RV(r_V(indexV1,:),v2_E1toV1,mu_sun);

for j = indexE1:indexV1
    % Calculate departure eccentric anomaly from departure true anomaly
    E_SC_E1 = 2*atan(sqrt((1-e_E1toV1)/(1+e_E1toV1))*tan(fdep_E1toV1/2)); % rad
    % Calculate departure mean anomaly from departure eccentric anomaly
    M_SC_E1 = E_SC_E1 - e_E1toV1*sin(E_SC_E1);
    % Calculate mean anomaly at each time point
    M_SC_E1toV1(j) = M_SC_E1 + n_E1toV1*(t(j)-Earth1JD)*24*60*60;
    % Calculate spacecraft position at each time point
    [r_SC(j,:),v_SC(j,:)] = RV_from_COE(a_E1toV1,e_E1toV1,i_E1toV1,RAAN_E1toV1,w_E1toV1,M_SC_E1toV1(j),mu_sun);
end
% for q = indexE1:10:indexV1
    % Plot XY
%     figure(1)
%     hold on;
%     plot(r_SC(q,1),r_SC(q,2),'k','markersize',10)
%     plot(r_E(q,1),r_E(q,2),'b','markersize',30)
%     plot(r_V(q,1),r_V(q,2),'r','markersize',20)
%     plot(r_sun(1),r_sun(2),'y.', 'markersize', 50)
%     title('Earth Departure to First Venus Flyby')
%     xlabel('X (km)')
%     ylabel('Y (km)')
%     legend('Transfer','Earth','Venus','Location','southeast')
% end

%% First Venus Flyby to Second Venus Flyby (V1toV2)
% Use Lambert function to find transfer orbit and departure and arrival
% velocities
[v1_V1toV2,v2_V1toV2] = lambert(r_V(indexV1,:),r_V(indexV2,:),TOF_V1toV2,'pro',mu_sun);
% Find COEs of transfer orbit
[a_V1toV2,e_V1toV2,i_V1toV2,RAAN_V1toV2,w_V1toV2,n_V1toV2,fdep_V1toV2,T_V1toV2] = COE_from_RV(r_V(indexV1,:),v1_V1toV2,mu_sun);
[a_V1toV2,e_V1toV2,i_V1toV2,RAAN_V1toV2,w_V1toV2,n_V1toV2,farr_V1toV2,T_V1toV2] = COE_from_RV(r_V(indexV2,:),v2_V1toV2,mu_sun);

for k = indexV1:indexV2
    % Calculate departure eccentric anomaly from departure true anomaly
    E_SC_V1 = 2*atan(sqrt((1-e_V1toV2)/(1+e_V1toV2))*tan(fdep_V1toV2/2)); % rad
    % Calculate departure mean anomaly from departure eccentric anomaly
    M_SC_V1 = E_SC_V1 - e_V1toV2*sin(E_SC_V1);
    % Calculate mean anomaly at each time point
    M_SC_V1toV2(k) = M_SC_V1 + n_V1toV2*(t(k)-Venus1JD)*24*60*60;
    % Calculate spacecraft position at each time point
    [r_SC(k,:),v_SC(k,:)] = RV_from_COE(a_V1toV2,e_V1toV2,i_V1toV2,RAAN_V1toV2,w_V1toV2,M_SC_V1toV2(k),mu_sun);
    % Plot XY
%     figure(2)
%     hold on;
%     plot(r_SC(k,1),r_SC(k,2),'k','markersize',10)
%     plot(r_E(k,1),r_E(k,2),'b','markersize',30)
%     plot(r_V(k,1),r_V(k,2),'r','markersize',20)
%     plot(r_sun(1),r_sun(2),'y.', 'markersize', 50)
%     title('First Venus Flyby to Second Venus Flyby')
%     xlabel('X (km)')
%     ylabel('Y (km)')
%     legend('Transfer','Earth','Venus','Location','southeast')
end

%% Second Venus Flyby to Earth Flyby (V2toE2)
% Use Lambert function to find transfer orbit and departure and arrival
% velocities
[v1_V2toE2,v2_V2toE2] = lambert(r_V(indexV2,:),r_E(indexE2,:),TOF_V2toE2,'pro',mu_sun);
% Find COEs of transfer orbit
[a_V2toE2,e_V2toE2,i_V2toE2,RAAN_V2toE2,w_V2toE2,n_V2toE2,fdep_V2toE2,T_V2toE2] = COE_from_RV(r_V(indexV2,:),v1_V2toE2,mu_sun);
[a_V2toE2,e_V2toE2,i_V2toE2,RAAN_V2toE2,w_V2toE2,n_V2toE2,farr_V2toE2,T_V2toE2] = COE_from_RV(r_E(indexE2,:),v2_V2toE2,mu_sun);

for l = indexV2:indexE2
    % Calculate departure eccentric anomaly from departure true anomaly
    E_SC_V2 = 2*atan(sqrt((1-e_V2toE2)/(1+e_V2toE2))*tan(fdep_V2toE2/2)); % rad
    % Calculate departure mean anomaly from departure eccentric anomaly
    M_SC_V2 = E_SC_V2 - e_V2toE2*sin(E_SC_V2);
    % Calculate mean anomaly at each time point
    M_SC_V2toE2(l) = M_SC_V2 + n_V2toE2*(t(l)-Venus2JD)*24*60*60;
    % Calculate spacecraft position at each time point
    [r_SC(l,:),v_SC(l,:)] = RV_from_COE(a_V2toE2,e_V2toE2,i_V2toE2,RAAN_V2toE2,w_V2toE2,M_SC_V2toE2(l),mu_sun);
    % Plot XY
%     figure(3)
%     hold on;
%     plot(r_SC(l,1),r_SC(l,2),'k','markersize',10)
%     plot(r_E(l,1),r_E(l,2),'b','markersize',30)
%     plot(r_V(l,1),r_V(l,2),'r','markersize',20)
%     plot(r_sun(1),r_sun(2),'y.', 'markersize', 50)
%     title('Second Venus Flyby to Earth Flyby')
%     xlabel('X (km)')
%     ylabel('Y (km)')
%     legend('Transfer','Earth','Venus','Location','southeast')
end

%% Earth Flyby to Jupiter Flyby (E2toJ)
% Use Lambert function to find transfer orbit and departure and arrival
% velocities
[v1_E2toJ,v2_E2toJ] = lambert(r_E(indexE2,:),r_J(indexJ,:),TOF_E2toJ,'pro',mu_sun);
% Find COEs of transfer orbit
[a_E2toJ,e_E2toJ,i_E2toJ,RAAN_E2toJ,w_E2toJ,n_E2toJ,fdep_E2toJ,T_E2toJ] = COE_from_RV(r_E(indexE2,:),v1_E2toJ,mu_sun);
[a_E2toJ,e_E2toJ,i_E2toJ,RAAN_E2toJ,w_E2toJ,n_E2toJ,farr_E2toJ,T_E2toJ] = COE_from_RV(r_J(indexJ,:),v2_E2toJ,mu_sun);

for m = indexE2:indexJ
    % Calculate departure eccentric anomaly from departure true anomaly
    E_SC_E2 = 2*atan(sqrt((1-e_E2toJ)/(1+e_E2toJ))*tan(fdep_E2toJ/2)); % rad
    % Calculate departure mean anomaly from departure eccentric anomaly
    M_SC_E2 = E_SC_E2 - e_E2toJ*sin(E_SC_E2);
    % Calculate mean anomaly at each time point
    M_SC_E2toJ(m) = M_SC_E2 + n_E2toJ*(t(m)-Earth2JD)*24*60*60;
    % Calculate spacecraft position at each time point
    [r_SC(m,:),v_SC(m,:)] = RV_from_COE(a_E2toJ,e_E2toJ,i_E2toJ,RAAN_E2toJ,w_E2toJ,M_SC_E2toJ(m),mu_sun);
    % Plot XY
%     figure(4)
%     hold on;
%     plot(r_SC(m,1),r_SC(m,2),'k','markersize',10)
%     plot(r_E(m,1),r_E(m,2),'b','markersize',30)
%     plot(r_J(m,1),r_J(m,2),'g','markersize',40)
%     plot(r_sun(1),r_sun(2),'y.', 'markersize', 50)
%     title('Earth Flyby to Jupiter Flyby')
%     xlabel('X (km)')
%     ylabel('Y (km)')
%     legend('Transfer','Earth','Jupiter','Location','southeast')
end

%% Jupiter Flyby to Saturn (JtoS)
% Use Lambert function to find transfer orbit and departure and arrival
% velocities
[v1_JtoS,v2_JtoS] = lambert(r_J(indexJ,:),r_S(indexS,:),TOF_JtoS,'pro',mu_sun);
% Find COEs of transfer orbit
[a_JtoS,e_JtoS,i_JtoS,RAAN_JtoS,w_JtoS,n_JtoS,fdep_JtoS,T_JtoS] = COE_from_RV(r_J(indexJ,:),v1_JtoS,mu_sun);
[a_JtoS,e_JtoS,i_JtoS,RAAN_JtoS,w_JtoS,n_JtoS,farr_JtoS,T_JtoS] = COE_from_RV(r_S(indexS,:),v2_JtoS,mu_sun);

for n = indexJ:indexS
    % Calculate departure eccentric anomaly from departure true anomaly
    E_SC_J = 2*atan(sqrt((1-e_JtoS)/(1+e_JtoS))*tan(fdep_JtoS/2)); % rad
    % Calculate departure mean anomaly from departure eccentric anomaly
    M_SC_J = E_SC_J - e_JtoS*sin(E_SC_J);
    % Calculate mean anomaly at each time point
    M_SC_JtoS(n) = M_SC_J + n_JtoS*(t(n)-JupiterJD)*24*60*60;
    % Calculate spacecraft position at each time point
    [r_SC(n,:),v_SC(n,:)] = RV_from_COE(a_JtoS,e_JtoS,i_JtoS,RAAN_JtoS,w_JtoS,M_SC_JtoS(n),mu_sun);
end
% for p = indexJ:10:indexS
%     Plot XY
%     figure(5)
%     hold on;
%     plot(r_SC(p,1),r_SC(p,2),'k','markersize',10)
%     plot(r_J(p,1),r_J(p,2),'g','markersize',40)
%     plot(r_S(p,1),r_S(p,2),'m','markersize',35)
%     plot(r_sun(1),r_sun(2),'y.', 'markersize', 50)
%     title('Jupiter Flyby to Saturn')
%     xlabel('X (km)')
%     ylabel('Y (km)')
%     legend('Transfer','Jupiter','Saturn','Location','southeast')
% end

%% Plot Entire Orbit in XY Plane
% for o = indexE1:50:indexS
%     figure(6)
%     hold on;
%     plot(r_SC(o,1),r_SC(o,2),'k.','markersize',10)
%     plot(r_E(o,1),r_E(o,2),'b.','markersize',10)
%     plot(r_V(o,1),r_V(o,2),'r.','markersize',7)
%     plot(r_J(o,1),r_J(o,2),'g.','markersize',15)
%     plot(r_S(o,1),r_S(o,2),'m.','markersize',12)
%     plot(r_sun(1),r_sun(2),'y.', 'markersize',50)
%     title('Cassini Mission Trajectory')
%     xlabel('X (km)')
%     ylabel('Y (km)')
%     legend('Cassini','Earth','Venus','Jupiter','Saturn','Location','southeast')
% end

%% Plot Entire Orbit in XZ Plane
% for o = indexE1:50:indexS
%     figure(7)
%     hold on;
%     plot(r_SC(o,1),r_SC(o,3),'k.','markersize',10)
%     plot(r_E(o,1),r_E(o,3),'b.','markersize',10)
%     plot(r_V(o,1),r_V(o,3),'r.','markersize',7)
%     plot(r_J(o,1),r_J(o,3),'g.','markersize',15)
%     plot(r_S(o,1),r_S(o,3),'m.','markersize',12)
%     plot(r_sun(1),r_sun(3),'y.', 'markersize',20)
%     title('Cassini Mission Trajectory')
%     xlabel('X (km)')
%     ylabel('Z (km)')
%     legend('Cassini','Earth','Venus','Jupiter','Saturn','Location','southeast')
% end

%% Plot Entire Orbit in YZ Plane
% for o = indexE1:50:indexS
%     figure(8)
%     hold on;
%     plot(r_SC(o,2),r_SC(o,3),'k.','markersize',10)
%     plot(r_E(o,2),r_E(o,3),'b.','markersize',10)
%     plot(r_V(o,2),r_V(o,3),'r.','markersize',7)
%     plot(r_J(o,2),r_J(o,3),'g.','markersize',15)
%     plot(r_S(o,2),r_S(o,3),'m.','markersize',12)
%     plot(r_sun(2),r_sun(3),'y.','markersize',20)
%     title('Cassini Mission Trajectory')
%     xlabel('Y (km)')
%     ylabel('Z (km)')
%     legend('Cassini','Earth','Venus','Jupiter','Saturn','Location','southeast')
% end

%% Earth Departure
disp('Earth Departure')
r_park = 500 + Re; % km
v_park = sqrt(mu_earth/r_park); % parking velocity at Earth
v1_E1toV1mag = norm(v1_E1toV1); % magnitude of heliocentric Earth departure velocity
v_E_departure = norm(v_E(indexE1,:)); % heliocentric velocity of Earth at departure
vinf_E1toV1 = v1_E1toV1mag - v_E_departure;
dv_departure = sqrt(vinf_E1toV1^2 + 2*mu_earth/r_park) - v_park;
fprintf('Departure velocity change: %.3f km/s\nV_infinity in the planetary frame: %.3f km/s\n\n',dv_departure,vinf_E1toV1)
%% First Venus Gravity Assist
dv_V1 = v1_V1toV2 - v2_E1toV1;
dv_V1mag = norm(dv_V1)
gamma_V1 = acos(dot(v1_V1toV2,v2_E1toV1)/norm(v1_V1toV2)/norm(v2_E1toV1))
v_V1mag = norm(v_V(indexV1,:));
vplus_V1 = norm(v1_V1toV2)
vinf_V1mag = sqrt(v_V1mag^2 + vplus_V1^2 - 2*v_V1mag*vplus_V1*cos(gamma_V1))
delta_V1 = 2*asin(dv_V1mag/(2*vinf_V1mag))

e_swing_V1 = 1/sin(delta_V1/2);
a_swing_V1 = -mu_sun/(vinf_V1mag^2/2)/2;
rf_swing_V1 = a_swing_V1*(1-e_swing_V1);

if rf_swing_V1 > Rv
    rf_test_V1 = 'OK';
else
    rf_test_V1 = 'Too Low';
end

%% Second Venus Gravity Assist
dv_V2 = v1_V2toE2 - v2_V1toV2;
dv_V2mag = norm(dv_V2)
gamma_V2 = acos(dot(v1_V2toE2,v2_V1toV2)/norm(v1_V2toE2)/norm(v2_V1toV2))
v_V2mag = norm(v_V(indexV2,:))
vplus_V2 = norm(v1_V2toE2);
vinf_V2mag = sqrt(v_V2mag^2 + vplus_V2^2 - 2*v_V2mag*vplus_V2*cos(gamma_V2))
delta_V2 = 2*asin(dv_V2mag/(2*vinf_V2mag))

e_swing_V2 = 1/sin(delta_V2/2);
a_swing_V2 = -mu_sun/(vinf_V2mag^2/2)/2;
rf_swing_V2 = a_swing_V2*(1-e_swing_V2);

if rf_swing_V2 > Rv
    rf_test_V2 = 'OK';
else
    rf_test_V2 = 'Too Low';
end

%% Earth Gravity Assist
dv_E2 = v1_E2toJ - v2_V2toE2;
dv_E2mag = norm(dv_E2)
gamma_E2 = acos(dot(v1_E2toJ,v2_V2toE2)/norm(v1_E2toJ)/norm(v2_V2toE2))
v_E2mag = norm(v_E(indexE2,:))
vplus_E2 = norm(v1_E2toJ);
vinf_E2mag = sqrt(v_E2mag^2 + vplus_E2^2 - 2*v_E2mag*vplus_E2*cos(gamma_E2))
delta_E2 = 2*asin(dv_E2mag/(2*vinf_E2mag))

e_swing_E2 = 1/sin(delta_E2/2);
a_swing_E2 = -mu_sun/(vinf_E2mag^2/2)/2;
rf_swing_E2 = a_swing_E2*(1-e_swing_E2);
if rf_swing_E2 > Re
    rf_test_E2 = 'OK';
else
    rf_test_E2 = 'Too Low';
end

%% Jupiter Gravity Assist
dv_J = v1_JtoS - v2_E2toJ;
dv_Jmag = norm(dv_J)
gamma_J = acos(dot(v1_JtoS,v2_E2toJ)/norm(v1_JtoS)/norm(v2_E2toJ))
v_Jmag = norm(v_J(indexJ,:))
vplus_J = norm(v1_JtoS)
vinf_Jmag = sqrt(v_Jmag^2 + vplus_J^2 - 2*v_Jmag*vplus_J*cos(gamma_J))
delta_J = 2*asin(dv_Jmag/(2*vinf_Jmag))

e_swing_J = 1/sin(delta_J/2);
a_swing_J = -mu_sun/(vinf_Jmag^2/2)/2;
rf_swing_J = a_swing_J*(1-e_swing_J);
if rf_swing_J > Rj
    rf_test_J = 'OK';
else
    rf_test_J = 'Too Low';
end

%% Saturn Gravity Assist to Attempt Solar System Escape
disp('Saturn Flyby to Attempt Solar System Escape')

v_helio_escape_mag = sqrt(2*mu_sun/norm(r_S(indexS,:))); % Magnitude of heliocentric escape velocity
dv_escape_Smag = v_helio_escape_mag - norm(v2_JtoS); % required velocity change to achieve heliocentric escape
e_escape_S = 1 + 1e-6; % e > 1 for hyperbolic escape trajectory
delta_escape_S = 2*asin(1/e_escape_S);
vinf_escape_Smag =dv_escape_Smag/2/sin(delta_escape_S/2);
a_escape_S = -mu_sun/(vinf_escape_Smag^2/2)/2;
rf_escape_S = a_escape_S*(1-e_escape_S);
fprintf('Required flyby radius for heliocentric escape: %.5e km\n',rf_escape_S);

if rf_escape_S > Rs
    disp('Flyby transfer radius OK, heliocentric escape achieved')
else
    disp('Error: flyby transfer radius is too low!')
    disp('Calculate maximum possible heliocentric orbit energy using')
    disp('Saturn''s radius as the flyby radius')
    % Maximum heliocentric energy
    rf_swing_S = Rs;
    vinf_swing_S = v2_JtoS - v_S(indexS,:);
    a_swing_S = -mu_sun/norm(vinf_swing_S);
    e_swing_S = 1-rf_swing_S/a_swing_S;
    delta_swing_S = 2*asin(1/e_swing_S);
    vplus_swing_S = sqrt(norm(v_S(indexS,:))^2 + norm(vinf_swing_S)^2 - 2*norm(v_S(indexS,:))*norm(vinf_swing_S)*cos(delta_swing_S));
    Emax_swing_S = vplus_swing_S^2/2 - mu_sun/norm(r_S(indexS,:));
    fprintf('Maximum heliocentric orbit energy: %.5e km\n\n',Emax_swing_S);
    % Compare with heliocentric energy before Saturn assist
    E_before_S = norm(v2_JtoS)^2/2 - mu_sun/norm(r_S(indexS,:));
end

%% Hohmann Transfer from Earth to Saturn
% Hohmann transfer assumes circular orbits, use semimajor axes
[dv_total,dv_hohmann_departure,dv_hohmann_arrival,ToF_hohmann] = HohmannTransfer(a_earth,a_saturn,mu_sun);
disp('Hohmann transfer orbit from Earth to Saturn');
disp('Note: for Hohmann transfers, orbits are assumed circular and coplanar');
fprintf('Departure velocity change: %.3f km/s\nTime of Flight: %.3f days\n\n',dv_hohmann_departure,ToF_hohmann);
disp('Cassini trajectory')
Cassini_TOF = indexS-indexE1; % days
fprintf('Departure velocity change: %.3f km/s\nTime of Flight: %.3f days\n\n',dv_departure,Cassini_TOF);

%% Construct tables for workspace variables
% E1toV1
data1 = {
    'a', a_E1toV1, 'km', 'E1toV1 - Semi-major axis';
    'e', e_E1toV1, '-', 'E1toV1 - Eccentricity';
    'i', i_E1toV1, 'rad', 'E1toV1 - Inclination';
    'RAAN',RAAN_E1toV1, 'rad', 'E1toV1 - Argument of Node';
    'w',w_E1toV1, 'rad', 'E1toV1 - Argument of Periapsis';
    'ToF', TOF_E1toV1, 'sec', 'E1toV1 - Time of Flight';
    'fd', fdep_E1toV1, 'rad', 'E1toV1 - Departure True Anomaly';
    'fa', farr_E1toV1, 'rad', 'E1toV1 - Arrival True Anomaly';
    'T', T_E1toV1, 'sec', 'E1toV1 - Orbit Period';
    'Rd_I', r_E(indexE1,1), 'km', 'E1toV1 - Departure Radius I-component';
    'Rd_J', r_E(indexE1,2), 'km', 'E1toV1 - Departure Radius J-component';
    'Rd_K', r_E(indexE1,3), 'km', 'E1toV1 - Departure Radius K-component';
    'Vd_I', v1_E1toV1(1), 'km/s', 'E1toV1 - Departure Velocity I-component';
    'Vd_J', v1_E1toV1(2), 'km/s', 'E1toV1 - Departure Velocity J-component';
    'Vd_K', v1_E1toV1(3), 'km/s', 'E1toV1 - Departure Velocity K-component';
    'Ra_I', r_V(indexV1,1), 'km', 'E1toV1 - Arrival Radius I-component';
    'Ra_J', r_V(indexV1,1), 'km', 'E1toV1 - Arrival Radius J-component';
    'Ra_K', r_V(indexV1,3), 'km', 'E1toV1 - Arrival Radius K-component';
    'Va_I', v2_E1toV1(1), 'km/s', 'E1toV1 - Arrival Velocity I-component';
    'Va_J', v2_E1toV1(2), 'km/s', 'E1toV1 - Arrival Velocity J-component';
    'Va_K', v2_E1toV1(3), 'km/s', 'E1toV1 - Arrival Velocity K-component';
    };
dfig1 = figure(9);
T1 = uitable('data',data1, 'ColumnName',{'variable', 'value', 'unit',...
    'description'},'ColumnWidth',{64,64,64,180});
reformatTable(T1,dfig1);

% V1toV2
data2 = {
    'a', a_V1toV2, 'km', 'V1toV2 - Semi-major axis';
    'e', e_V1toV2, '-', 'V1toV2 - Eccentricity';
    'i', i_V1toV2, 'rad', 'V1toV2 - Inclination';
    'RAAN',RAAN_V1toV2, 'rad', 'V1toV2 - Argument of Node';
    'w',w_V1toV2, 'rad', 'V1toV2 - Argument of Periapsis';
    'ToF', TOF_V1toV2, 'sec', 'V1toV2 - Time of Flight';
    'fd', fdep_V1toV2, 'rad', 'V1toV2 - Departure True Anomaly';
    'fa', farr_V1toV2, 'rad', 'V1toV2 - Arrival True Anomaly';
    'T', T_V1toV2, 'sec', 'V1toV2 - Orbit Period';
    'Rd_I', r_V(indexV1,1), 'km', 'V1toV2 - Departure Radius I-component';
    'Rd_J', r_V(indexV1,2), 'km', 'V1toV2 - Departure Radius J-component';
    'Rd_K', r_V(indexV1,3), 'km', 'V1toV2 - Departure Radius K-component';
    'Vd_I', v1_V1toV2(1), 'km/s', 'V1toV2 - Departure Velocity I-component';
    'Vd_J', v1_V1toV2(2), 'km/s', 'V1toV2 - Departure Velocity J-component';
    'Vd_K', v1_V1toV2(3), 'km/s', 'V1toV2 - Departure Velocity K-component';
    'Ra_I', r_V(indexV2,1), 'km', 'V1toV2 - Arrival Radius I-component';
    'Ra_J', r_V(indexV2,1), 'km', 'V1toV2 - Arrival Radius J-component';
    'Ra_K', r_V(indexV2,3), 'km', 'V1toV2 - Arrival Radius K-component';
    'Va_I', v2_V1toV2(1), 'km/s', 'V1toV2 - Arrival Velocity I-component';
    'Va_J', v2_V1toV2(2), 'km/s', 'V1toV2 - Arrival Velocity J-component';
    'Va_K', v2_V1toV2(3), 'km/s', 'V1toV2 - Arrival Velocity K-component';
    };
dfig2 = figure(10);
T2 = uitable('data',data2, 'ColumnName',{'variable', 'value', 'unit',...
    'description'},'ColumnWidth',{64,64,64,180});
reformatTable(T2,dfig2);

% V2toE2
data3 = {
    'a', a_V2toE2, 'km', 'V2toE2 - Semi-major axis';
    'e', e_V2toE2, '-', 'V2toE2 - Eccentricity';
    'i', i_V2toE2, 'rad', 'V2toE2 - Inclination';
    'RAAN',RAAN_V2toE2, 'rad', 'V2toE2 - Argument of Node';
    'w',w_V2toE2, 'rad', 'V2toE2 - Argument of Periapsis';
    'ToF', TOF_V2toE2, 'sec', 'V2toE2 - Time of Flight';
    'fd', fdep_V2toE2, 'rad', 'V2toE2 - Departure True Anomaly';
    'fa', farr_V2toE2, 'rad', 'V2toE2 - Arrival True Anomaly';
    'T', T_V2toE2, 'sec', 'V2toE2 - Orbit Period';
    'Rd_I', r_V(indexV2,1), 'km', 'V2toE2 - Departure Radius I-component';
    'Rd_J', r_V(indexV2,2), 'km', 'V2toE2 - Departure Radius J-component';
    'Rd_K', r_V(indexV2,3), 'km', 'V2toE2 - Departure Radius K-component';
    'Vd_I', v1_V2toE2(1), 'km/s', 'V2toE2 - Departure Velocity I-component';
    'Vd_J', v1_V2toE2(2), 'km/s', 'V2toE2 - Departure Velocity J-component';
    'Vd_K', v1_V2toE2(3), 'km/s', 'V2toE2 - Departure Velocity K-component';
    'Ra_I', r_E(indexE2,1), 'km', 'V2toE2 - Arrival Radius I-component';
    'Ra_J', r_E(indexE2,1), 'km', 'V2toE2 - Arrival Radius J-component';
    'Ra_K', r_E(indexE2,3), 'km', 'V2toE2 - Arrival Radius K-component';
    'Va_I', v2_V2toE2(1), 'km/s', 'V2toE2 - Arrival Velocity I-component';
    'Va_J', v2_V2toE2(2), 'km/s', 'V2toE2 - Arrival Velocity J-component';
    'Va_K', v2_V2toE2(3), 'km/s', 'V2toE2 - Arrival Velocity K-component';
    };
dfig3 = figure(11);
T3 = uitable('data',data3, 'ColumnName',{'variable', 'value', 'unit',...
    'description'},'ColumnWidth',{64,64,64,180});
reformatTable(T3,dfig3);

% E2toJ
data4 = {
    'a', a_E2toJ, 'km', 'E2toJ - Semi-major axis';
    'e', e_E2toJ, '-', 'E2toJ - Eccentricity';
    'i', i_E2toJ, 'rad', 'E2toJ - Inclination';
    'RAAN',RAAN_E2toJ, 'rad', 'E2toJ - Argument of Node';
    'w',w_E2toJ, 'rad', 'E2toJ - Argument of Periapsis';
    'ToF', TOF_E2toJ, 'sec', 'E2toJ - Time of Flight';
    'fd', fdep_E2toJ, 'rad', 'E2toJ - Departure True Anomaly';
    'fa', farr_E2toJ, 'rad', 'E2toJ - Arrival True Anomaly';
    'T', T_E2toJ, 'sec', 'E2toJ - Orbit Period';
    'Rd_I', r_E(indexE2,1), 'km', 'E2toJ - Departure Radius I-component';
    'Rd_J', r_E(indexE2,2), 'km', 'E2toJ - Departure Radius J-component';
    'Rd_K', r_E(indexE2,3), 'km', 'E2toJ - Departure Radius K-component';
    'Vd_I', v1_E2toJ(1), 'km/s', 'E2toJ - Departure Velocity I-component';
    'Vd_J', v1_E2toJ(2), 'km/s', 'E2toJ - Departure Velocity J-component';
    'Vd_K', v1_E2toJ(3), 'km/s', 'E2toJ - Departure Velocity K-component';
    'Ra_I', r_J(indexJ,1), 'km', 'E2toJ - Arrival Radius I-component';
    'Ra_J', r_J(indexJ,1), 'km', 'E2toJ - Arrival Radius J-component';
    'Ra_K', r_J(indexJ,3), 'km', 'E2toJ - Arrival Radius K-component';
    'Va_I', v2_E2toJ(1), 'km/s', 'E2toJ - Arrival Velocity I-component';
    'Va_J', v2_E2toJ(2), 'km/s', 'E2toJ - Arrival Velocity J-component';
    'Va_K', v2_E2toJ(3), 'km/s', 'E2toJ - Arrival Velocity K-component';
    };
dfig4 = figure(12);
T4 = uitable('data',data4, 'ColumnName',{'variable', 'value', 'unit',...
    'description'},'ColumnWidth',{64,64,64,180});
reformatTable(T4,dfig4);

% JtoS
data5 = {
    'a', a_JtoS, 'km', 'JtoS - Semi-major axis';
    'e', e_JtoS, '-', 'JtoS - Eccentricity';
    'i', i_JtoS, 'rad', 'JtoS - Inclination';
    'RAAN',RAAN_JtoS, 'rad', 'JtoS - Argument of Node';
    'w',w_JtoS, 'rad', 'JtoS - Argument of Periapsis';
    'ToF', TOF_JtoS, 'sec', 'JtoS - Time of Flight';
    'fd', fdep_JtoS, 'rad', 'JtoS - Departure True Anomaly';
    'fa', farr_JtoS, 'rad', 'JtoS - Arrival True Anomaly';
    'T', T_JtoS, 'sec', 'JtoS - Orbit Period';
    'Rd_I', r_J(indexJ,1), 'km', 'JtoS - Departure Radius I-component';
    'Rd_J', r_J(indexJ,2), 'km', 'JtoS - Departure Radius J-component';
    'Rd_K', r_J(indexJ,3), 'km', 'JtoS - Departure Radius K-component';
    'Vd_I', v1_JtoS(1), 'km/s', 'JtoS - Departure Velocity I-component';
    'Vd_J', v1_JtoS(2), 'km/s', 'JtoS - Departure Velocity J-component';
    'Vd_K', v1_JtoS(3), 'km/s', 'JtoS - Departure Velocity K-component';
    'Ra_I', r_S(indexS,1), 'km', 'JtoS - Arrival Radius I-component';
    'Ra_J', r_S(indexS,1), 'km', 'JtoS - Arrival Radius J-component';
    'Ra_K', r_S(indexS,3), 'km', 'JtoS - Arrival Radius K-component';
    'Va_I', v2_JtoS(1), 'km/s', 'JtoS - Arrival Velocity I-component';
    'Va_J', v2_JtoS(2), 'km/s', 'JtoS - Arrival Velocity J-component';
    'Va_K', v2_JtoS(3), 'km/s', 'JtoS - Arrival Velocity K-component';
    };
dfig5 = figure(13);
T5 = uitable('data',data5, 'ColumnName',{'variable', 'value', 'unit',...
    'description'},'ColumnWidth',{64,64,64,180});
reformatTable(T5,dfig5);

% First Venus Gravity Assist
data6 = {
    'dV_I', dv_V1(1), 'km/s', 'First Venus Flyby - Departure Velocity I-component';
    'dV_J', dv_V1(2), 'km/s', 'First Venus Flyby - Departure Velocity J-component';
    'dV_K', dv_V1(3), 'km/s', 'First Venus Flyby - Departure Velocity K-component';
    'Vinf', vinf_V1mag, 'km/s', 'First Venus Flyby - V_infinity in planetary frame';
    'Rf', rf_swing_V1, 'km', 'First Venus Flyby - Required Flyby Altitude';
    '-', rf_test_V1, '-', 'Reasonable flyby altitude?';
    };
dfig6 = figure(14);
T6 = uitable('data',data6, 'ColumnName',{'variable', 'value', 'unit',...
    'description'},'ColumnWidth',{64,64,64,180});
reformatTable(T6,dfig6);

% Second Venus Gravity Assist
data7 = {
    'dV_I', dv_V2(1), 'km/s', 'Second Venus Flyby - Departure Velocity I-component';
    'dV_J', dv_V2(2), 'km/s', 'Second Venus Flyby - Departure Velocity J-component';
    'dV_K', dv_V2(3), 'km/s', 'Second Venus Flyby - Departure Velocity K-component';
    'Vinf', vinf_V2mag, 'km/s', 'Second Venus Flyby - V_infinity in planetary frame';
    'Rf', rf_swing_V2, 'km', 'Second Venus Flyby - Required Flyby Altitude';
    '-', rf_test_V2, '-', 'Reasonable flyby altitude?';
    };
dfig7 = figure(15);
T7 = uitable('data',data7, 'ColumnName',{'variable', 'value', 'unit',...
    'description'},'ColumnWidth',{64,64,64,180});
reformatTable(T7,dfig7);

% Earth Gravity Assist
data8 = {
    'dV_I', dv_E2(1), 'km/s', 'Earth Flyby - Departure Velocity I-component';
    'dV_J', dv_E2(2), 'km/s', 'Earth Flyby - Departure Velocity J-component';
    'dV_K', dv_E2(3), 'km/s', 'Earth Flyby - Departure Velocity K-component';
    'Vinf', vinf_E2mag, 'km/s', 'Earth Flyby - V_infinity in planetary frame';
    'Rf', rf_swing_E2, 'km', 'Earth Flyby - Required Flyby Altitude';
    '-', rf_test_E2, '-', 'Reasonable flyby altitude?';
    };
dfig8 = figure(16);
T8 = uitable('data',data8, 'ColumnName',{'variable', 'value', 'unit',...
    'description'},'ColumnWidth',{64,64,64,180});
reformatTable(T8,dfig8);

% Jupiter Gravity Assist
data9 = {
    'dV_I', dv_J(1), 'km/s', 'Jupiter Flyby - Departure Velocity I-component';
    'dV_J', dv_J(2), 'km/s', 'Jupiter Flyby - Departure Velocity J-component';
    'dV_K', dv_J(3), 'km/s', 'Jupiter Flyby - Departure Velocity K-component';
    'Vinf', vinf_Jmag, 'km/s', 'Jupiter Flyby - V_infinity in planetary frame';
    'Rf', rf_swing_J, 'km', 'Jupiter Flyby - Required Flyby Altitude';
    '-', rf_test_J, '-', 'Reasonable flyby altitude?';
    };
dfig9 = figure(17);
T9 = uitable('data',data9, 'ColumnName',{'variable', 'value', 'unit',...
    'description'},'ColumnWidth',{64,64,64,180});
reformatTable(T9,dfig9);
