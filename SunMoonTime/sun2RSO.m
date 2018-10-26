function r_sun2RSO_unit = sun2RSO(year,mo,d,h,mins,s,r_RSO,v_RSO)

% Calculate Julian date
JulianDate = JD(year,mo,d,h,mins,s);
JD_UT1 = JulianDate;

% Calculate sun vector from Earth to Sun in AU units
r_E2sun = sunvector(JD_UT1); % in AUs (1 AU = 149,597,870 km)

% Define AU unit
AU = 149597870; % km


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

% Find sun vector from RSO to sun in relative frame
%r_RSO2sun_ECI_0 = r_EtoS_0-r_RSO/AU;
RSO2Sun = ON*(r_E2sun-r_RSO/AU);

% Flip direction to find vector from sun to RSO in relative frame
r_sun2RSO = -RSO2Sun;
r_sun2RSO_unit = r_sun2RSO/norm(r_sun2RSO); % This is sun shine angle for current time