function TargetInfo = Target(Omega,Marks,beta,i,NumberMarks,Order,TransferTimes)
% Converts relative orbital elements (as defined by Lovell) to
% Hill-Clohessy-Wiltshire (HCW)
% n = mean-motion of chief (1/sec)

% propagate mark from initial initial state to prox ops state

index = find(Order == i);
propagateTime = sum(TransferTimes(1:index+1));
[x_drift,v_drift] = CWHPropagator(Marks{i,8}(1:3)',Marks{i,8}(4:6)',Omega,propagateTime);
Marks{i,8} = [x_drift',v_drift'];

    

Type = Marks{i,1};

switch Type
    case 'NMC'
        n = Omega;
        ae = Marks{i,2};
        xd = Marks{i,8}(1);
        yd = Marks{i,8}(2);
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