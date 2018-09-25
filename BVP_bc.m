% Boundary conditions
% x(1)=x position, x(2)=x velocity, x(3)=y position, x(4)=y velocity
% [x y z xdot ydot zdot]

function res = BVP_bc(xi,xf,CurrentState,TargetInfo,omega)

x_init = CurrentState;
x_target = TargetInfo;

res = [xi(1)-x_init(1)
    xi(2)-x_init(4)
    xi(3)-x_init(2)
    xi(4)-x_init(5)
    xf(1)-x_target(1)
    xf(2)-x_target(4)
    xf(3)-x_target(2)
    xf(4)-x_target(5)];
