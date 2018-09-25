function [u_tot] = ParameterizedFitness_BVP_rendezvous(dvar,obj_traj,omega,num_objects)
%MYFUNINT1 integration function for piecewise constant control


order = dvar(1:num_objects);
t_times = dvar(num_objects+1:end);
sat_state = [0 0 0 0 0 0];
u_tot = 0;
n = 1;
clock = 0;

for i = 1:length(order)
    
    
    % drift sat for first time index at current state
    [x_drift,v_drift] = CWHPropagator(sat_state(1:3)',sat_state(4:6)',omega,t_times(n));
    sat_state = [x_drift',v_drift'];
    
    % target obj we are travling to already drifted
    x_target = obj_traj{order(i)}(1:6,floor(clock+t_times(n+1)+t_times(n)))';
    % current location of the sat, not drifted
    x_init = sat_state(1,:);
    
    solinit = bvpinit(linspace(0,t_times(n+1),150),[0 0 0 0 0 0 0 0]);
    options = bvpset('RelTol',1e-6);

    sol = bvp4c(@BVP_ode,@BVP_bc,solinit,options,x_target,x_init,omega);
    t=[sol.x];
    x=[sol.y]; 

    sat_state(1,:) = [x(1,end) x(3,end) 0 x(2,end) x(4,end) 0];
    u_x = [-x(6,:)];
    u_y = [-x(8,:)];
    u_tot = trapz(0.5*(u_x.^2+u_y.^2)) + u_tot;
    
    clock = clock + t_times(n) + t_times(n+1);
    n = n + 2;
    
end

% drift sat for first time index at current state
[x_drift,v_drift] = CWHPropagator(sat_state(1:3)',sat_state(4:6)',omega,t_times(n));
sat_state = [x_drift',v_drift'];

% target obj we are travling to already drifted
x_target = [0 0 0 0 0 0];
% current location of the sat, not drifted
x_init = sat_state(1,:);

solinit = bvpinit(linspace(0,t_times(n+1),150),[0 0 0 0 0 0 0 0]);
options = bvpset('RelTol',1e-6);

sol = bvp4c(@BVP_ode,@BVP_bc,solinit,options,x_target,x_init,omega);
t=[sol.x];
x=[sol.y]; 

sat_state(1,:) = [x(1,end) x(3,end) 0 x(2,end) x(4,end) 0];
u_x = [-x(6,:)];
u_y = [-x(8,:)];
u_tot = trapz(0.5*(u_x.^2+u_y.^2))+u_tot;

 