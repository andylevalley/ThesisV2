function TargetInfo = Target(Omega,Marks,beta,i)
% Converts relative orbital elements (as defined by Lovell) to
% Hill-Clohessy-Wiltshire (HCW)
% n = mean-motion of chief (1/sec)

Type = Marks{i,1};

switch Type
    case 'NMC'
        n = Omega;
        ae = Marks{i,2};
        xd = Marks{i,3};
        yd = Marks{i,4};
        zmax = Marks{i,5};
        upsilon = Marks{i,6};
        beta = beta(i);

        x = (-ae/2)*cos(beta) + xd;
        y = ae*sin(beta) + yd;
        z = zmax*sin(upsilon);
        x_dot = (ae/2)*n*sin(beta);
        y_dot = ae*n*cos(beta) - 3*n*xd/2;
        z_dot = zmax*n*cos(upsilon);
        TargetInfo = [x y z x_dot y_dot z_dot];
end
end