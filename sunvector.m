function r_EtoS = sunvector(JD_UT1)

% Taken from VALLADO, algorithm 29

T_UT1 = (JD_UT1-2451545)/36525;
lambda_M_S = 280.460+36000.771*T_UT1;
T_TDB = T_UT1;
M_S = 357.5277233+35999.05034*T_TDB;
lambda_ecliptic = lambda_M_S+1.914666471*sind(M_S)+.019994643*sind(2*M_S);
r_S = 1.000140612-.016708617*cosd(M_S)-.000139589*cosd(2*M_S);
epsangle = 23.439291-.0130042*T_TDB;
r_EtoS = [r_S*cosd(lambda_ecliptic);
          r_S*cosd(epsangle)*sind(lambda_ecliptic);
          r_S*sind(epsangle)*sind(lambda_ecliptic)];
      