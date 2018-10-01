function X = LROE2X(a_e,x_d,y_d0,z_max,gamma,beta,w)

x = -a_e/2*cos(beta)+x_d;
y_d = y_d0-3/2*w*x_d*beta/w;
y = a_e*sin(beta)+y_d;
z = z_max*sin(gamma+beta);
xd = a_e/2*w*sin(beta);
yd = a_e*w*cos(beta)-3/2*w*x_d;
zd = z_max*w*cos(gamma+beta);

X = [x;y;z;xd;yd;zd];