function phaseout = GPOPS_Continuous(input)

a0 = input.auxdata.a0;
c = input.auxdata.c;
w = input.auxdata.w;

t     = input.phase(1).time;

x     = input.phase(1).state(:,1);
y     = input.phase(1).state(:,2);
z     = input.phase(1).state(:,3);
xdot  = input.phase(1).state(:,4);
ydot  = input.phase(1).state(:,5);
zdot  = input.phase(1).state(:,6);
alpha = input.phase(1).state(:,7);
beta  = input.phase(1).state(:,8);

alphadot = input.phase(1).control(:,1);
betadot  = input.phase(1).control(:,2);
th = input.phase(1).control(:,3);

% Now perform calculations with inputs ***********************************
% a       = T./(m0 - dm.*t); % acceleration
th_t = cumtrapz(t,th); % to account for when throttle is turned off
a = a0./(1-th_t.*(a0./c));
% a = a0./(1-t1.*(a0./c));

% HERE'S THE DYNAMICS %
dx = xdot;
dy = ydot;
dz = zdot;
dxdot = 2*w*ydot+3*w^2*x+th.*a.*cos(beta).*cos(alpha);
dydot = -2*w*xdot + th.*a.*cos(beta).*sin(alpha);
dzdot = -w^2*z+th.*a.*sin(beta);
dalpha = alphadot;
dbeta = betadot;


phaseout(1).dynamics  = [dx,dy,dz,dxdot,dydot,dzdot,dalpha,dbeta];
phaseout(1).integrand = th.^2;

end


