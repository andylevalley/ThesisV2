% Boundary conditions


function res = BVP_bc(xi,xf,CurrentState,TargetInfo,omega)

x_init = CurrentState;
x_target = TargetInfo;

res = [xi(1)-x_init(1)
    xi(2)-x_init(2)
    xi(3)-x_init(3)
    xi(4)-x_init(4)
    xi(5)-x_init(5)
    xi(6)-x_init(6)
    xf(1)-x_target(1)
    xf(2)-x_target(2)
    xf(3)-x_target(3)
    xf(4)-x_target(4)
    xf(5)-x_target(5)
    xf(6)-x_target(6)];
