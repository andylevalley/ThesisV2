% BVP ODES, HCW rendezvous



function dxdt = BVP_ode(t,x,CurrentState,TargetInfo,omega)


dxdt = [x(4) % x
        x(5) % y
        x(6) % z
        3*omega*omega*x(1)+2*omega*x(5) + 0.5*x(10) % xdot
        -2*omega*x(4) + 0.5*x(11) % ydot
        -omega^2*x(3) + 0.5*x(12) % zdot
        -3*omega^2*x(10) % lambda_x
        0 % lambda_y
        omega^2*x(12) % lambda_z
        -x(7)+2*omega*x(11)
        -x(8)-2*omega*x(10)
        -x(9)];
end
