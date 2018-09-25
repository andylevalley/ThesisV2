function r_RSO2moon_unit = RSO2moon(JD,r_RSO,v_RSO)

% Calculate moon vector from Earth to Moon in Earth radii (ER) units
[r_E2moon, ~, ~] = moon(JD); % in ERs (1 ER = 6378.137 km)

% Define AU unit
ER = 6378.137; % km

% FROM p. 604 of SCHAUB
% Given current position and velocity of RSO, calculate DCM:
% RSO_nu0 = 0; % deg
% r_RSO = [sma;0;0];
% v_RSO = [0;sqrt(mu/r_RSO);0];
% Momentum vector:
ohat_r = r_RSO/norm(r_RSO);
h_RSO = cross(r_RSO,v_RSO);
ohat_h = h_RSO/norm(h_RSO);
ohat_theta = cross(ohat_h,ohat_r);
ON = [ohat_r';ohat_theta';ohat_h'];

% Find moon vector from RSO to moon in relative frame
%r_RSO2moon_ECI_0 = r_EtoS_0-r_RSO/ER;
r_RSO2moon = ON*(r_E2moon'-r_RSO/ER); % this is in ERs
r_RSO2moon_unit = r_RSO2moon/norm(r_RSO2moon);

% % Flip direction to find vector from moon to RSO in relative frame
% r_moon2RSO = -r_RSO2moon;
% r_moon2RSO_unit = r_moon2RSO/norm(r_moon2RSO); % This is moon shine vector for current time