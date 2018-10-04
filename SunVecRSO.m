function [RSO2Sun] = SunVecRSO(time)

% current working epoch 
epoch = [2011 08 16 17 00 00];

%Compute the Julian Date (JD) from a Gregorian Date (UT1) input
yr=epoch(1);  %(year)
mo=epoch(2);  %(month)
d=epoch(3);  %(day)
h=epoch(4);  %(hr)
min=epoch(5);  %(min)
s=epoch(6);  %(sec)
JD=367*yr-floor(7/4*(yr+floor((mo+9)/12)))+floor(275*mo/9)+d+1721013.5+(((s/60+min)/60)+h)/24;  %(JD)

Re = 6378.1363e3;   %Earth's radius (m)

%%%Placeholder for propagation method
w = 2*pi/(24*60*60); % for GEO rad/s
a = 42164*10^3; % semi-major axis (m) 
X_ECI = a*cos(w*time); % (m)
Y_ECI = a*sin(w*time); % (m)
Z_ECI = 0;     % (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CWH_origin_position = [X_ECI Y_ECI Z_ECI];
%%%

%%%%% ECI to LVLH (RIC) conversions
%x along position vector
%y along velocity vector
%z along orbit normal
r = CWH_origin_position.';
v = (cross(r,[0 0 -1])).' ;%V_CWH_Origin;
r_norm = norm(r);
x = r/r_norm;
h_norm = norm(skew(r)*v);
z = skew(r)*v/h_norm;
y = skew(z)*x;
R_LVLHtoECI = [x y z];
DCM_ECItoLVLH = R_LVLHtoECI.';
%%%%%

%Compute the Universal Time (UT) from JD
T_UT1 = (JD-2451545.0)/36525;

%Compute the Greenwich Mean Sidereal Time (GMST) from UT
GMST = 67310.54841 + (876600.0*3600.0 + 8640184.812866)*T_UT1 + 0.093104*T_UT1^2 - 6.2E-6*T_UT1^3;
G0 = rem(GMST,86400)/240;
if G0 < 0
    G0 = -(-G0-360);
end
%truncate GMST between 0 and 2pi
GMST = G0*pi/180;  %(rad)

%Compute the Direction Cosine Matrix (DCM) for converting vectors in 
%Earth-Centered Inertial (ECI) frame to vectors in Earth-Centered 
%Earth-Fixed (ECEF) frame
R_ECItoECEF = ROT3(GMST);

%Compute the Sun vector based on the UT
RSO2Sun = DCM_ECItoLVLH*sunvector(T_UT1);


end

function A = ROT3(theta)

A = [cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
     0          0          1];
 
end

function r = r_sun(T_UT1)

lambdaM = 280.4606184+36000.77005361*T_UT1;  %(°)
M = (357.5277233+35999.05034*T_UT1)*pi/180.0;  %(rad)
lambdaecliptic = (lambdaM+1.914666471*sin(M)+0.019994643*sin(2*M))*pi/180.0;  %(rad)
E = (23.439291-0.0130042*T_UT1)*pi/180.0;  %(rad)
Rs = (1.000140612-0.016708617*cos(M)-0.000139589*cos(2*M))*149598000*10^3;  %(m)
r = Rs*[cos(lambdaecliptic);
        cos(E)*sin(lambdaecliptic);
        sin(E)*sin(lambdaecliptic)];  %(m)
  
end

function A = skew(x)

A = [0    -x(3)  x(2);
     x(3)  0    -x(1);
    -x(2)  x(1)  0];

end
