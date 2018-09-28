% ODE For HCW States and Co-States Equations: Coplanar Case, Control Input applied along x and y axis 
% x(1)=x1, x(2)=x2, x(3)=x3, x(4)=x4, x(5)=L1, x(6)=L2, x(7)=L3, x(8)=L4

function dxdt = BVP_ode(t,x,CurrentState,TargetInfo,omega)

% a0 = .0686; % N/kg - initial acceleration
% c = 3000; % m/s
% a = a0/(1-t*(a0/c));
% 
% dxdt = [x(2)
%     3*omega*omega*x(1)+2*omega*x(4) - x(6).*a
%     x(4)
%     -2*omega*x(2) - x(8).*a
%     -3*omega*omega*x(6)
%     -x(5)+2*omega*x(8)
%     0
%     -2*omega*x(6)-x(7)];

dxdt = [x(4)
        x(5)
        x(6)
        3*omega*omega*x(1)+2*omega*x(5) + 0.5*x(10)
        -2*omega*x(4) + 0.5*x(11)
        -omega^2*x(3) + 0.5*x(12)
        -3*omega^2*x(10)
        0
        omega^2*x(12)
        -x(7)+2*omega*x(11)
        -x(8)-2*omega*x(10)
        -x(9)];
end
