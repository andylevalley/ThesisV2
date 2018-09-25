function jd = JD(yr,mo,d,h,mins,s)

% Taken from VALLADO, algorithm 14

sixty = 60;
% UNLESS day contains a leap second, then sixty = 61;

jd = 367*yr-floor(7*(yr+floor((mo+9)/12))/4)+floor(275*mo/9)+d+1721013.5+((s/sixty+mins)/60+h)/24;
