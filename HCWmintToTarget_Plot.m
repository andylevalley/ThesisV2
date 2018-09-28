%------------------------------%
% Extract Solution from Output %
%------------------------------%
close all;
solution = output.result.solution;
time = solution.phase(1).time/60; % in MINUTES

w = auxdata.w;

x     = solution.phase(1).state(:,1)./1000; % scaled to km
y     = solution.phase(1).state(:,2)./1000; % scaled to km
z     = solution.phase(1).state(:,3)./1000; % scaled to km
xdot  = solution.phase(1).state(:,4); % m/s
ydot  = solution.phase(1).state(:,5); % m/s
zdot  = solution.phase(1).state(:,6); % m/s

u1  = solution.phase(1).control(:,1); % xhat direction of thrust
u2  = solution.phase(1).control(:,2); % yhat direction of thrust
u3  = solution.phase(1).control(:,3); % zhat direction of thrust

% r         = solution.phase(1).state(:,1);
% theta     = solution.phase(1).state(:,2);
% vr        = solution.phase(1).state(:,3);
% vtheta    = solution.phase(1).state(:,4);
% u1        = solution.phase(1).control(:,1);
% u2        = solution.phase(1).control(:,2);

% alpha     = unwrap(atan2(u1,u2))*180/pi; % Do I care about something
% similar with my problem?

% state = [r, theta, vr, vtheta];
% control = [u1, u2];

state = [x, y, z, xdot, ydot, zdot]; 
control = [u1, u2, u3];
u_tot = trapz(sqrt(u1.^2+u2.^2));

% Natural motion resulting from initial conditions and zero control
X = xdot(1)/1000/w*sin(w*time)-(3*x(1)+2*ydot(1)/1000/w)*cos(w*time)+(4*x(1)+2*ydot(1)/1000/w);
Y = (6*x(1)+4*ydot(1)/1000/w)*sin(w*time)+2*xdot(1)/1000/w*cos(w*time)-(6*w*x(1)+3*ydot(1)/1000)*time+(y(1)-2*xdot(1)/1000/w);
Z = z(1)*cos(w*time)+zdot(1)/1000/w*sin(w*time);

%---------------%
% Plot Solution %
%---------------%
figure(1)
pp = plot(time,state,'-o');
xl = xlabel('$t$ (min)','Interpreter','LaTeX');
% yl = ylabel('$(r(t),\theta(t),v_r(t),v_\theta(t))$','Interpreter','LaTeX');
yl = ylabel('$x,y,z,\dot{x},\dot{y},\dot{z}$ in $km$ and $\frac{m}{s}$','Interpreter','LaTeX');
% ll = legend('$r(t)$','$\theta(t)$','$v_r(t)$','$v_\theta(t)$','Location','Northwest');
ll = legend('$x$','$y$','$z$','$\dot{x}$','$\dot{y}$','$\dot{z}$','Location','Northwest');
tt = title('States vs. Time');
set(tt,'FontSize',18);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',18,'Interpreter','LaTeX');
set(gca,'FontSize',16,'FontName','Times');
set(pp,'LineWidth',1.25);
grid on
print -depsc2 HCWmintToTarget_State.eps;
print -dpng HCWmintToTarget_State.png;

vidObj1 = VideoWriter('3D_Trajectory');
vidObj1.FrameRate = 10;
% vidObj1.Quality = 100;
open(vidObj1);

t_int = output.result.interpsolution.phase.time;
x_int = output.result.interpsolution.phase.state(:,1)/1000; % km
y_int = output.result.interpsolution.phase.state(:,2)/1000; % km
z_int = output.result.interpsolution.phase.state(:,3)/1000; % km

for i=2:length(t_int)
    t_diffs(i-1) = t_int(i)-t_int(i-1); % in seconds
end

tf = t_int(end); % in seconds

% Resolution of plot - plot every max of t_diffs
t_diff_max = max(t_diffs);
t_plot = 0:t_diff_max:tf;

% Find equally spaced (in time) x, y, and z
xq = interp1(t_int,x_int,t_plot);
yq = interp1(t_int,y_int,t_plot);
zq = interp1(t_int,z_int,t_plot);

azvec = linspace(70,6,length(t_plot));
% elvec = linspace(-20,-10,length(t_plot));
elvec = linspace(-20,-20,length(t_plot));


for i = 1:length(t_plot)
    figure(2);
    if i == 1;
        pp = plot3(xq(1:i),yq(1:i),zq(1:i),'d');
    else
        pp = plot3(xq(1:i),yq(1:i),zq(1:i));
    end
    hold on;
    target = plot3(0,0,0,'x','MarkerSize',15);
    hold off;
    grid on;
    xl = xlabel('x (km)');
    yl = ylabel('y (km)');
    zl = zlabel('z (km)');
    tt = title('3-D Trajectory of Inspector');
%     ll = legend('Trajectory','Target');
    set(xl,'FontSize',18);
    set(yl,'FontSize',18);
    set(zl,'FontSize',18);
    set(tt,'FontSize',18);
    set(gca,'FontSize',16,'FontName','Times');
    set(pp,'LineWidth',1.25);
    set(target,'LineWidth',1.25);
    axis equal
    axis([-40 0 -100 0 0 10]);
    view([azvec(i),elvec(i)]);
    
    drawnow
    Frame = getframe(gcf);
    writeVideo(vidObj1,Frame);
end
    
close(vidObj1);
    

figure(3);
% p = plot3(X,Y,Z,'-d'); hold on;
pp = plot3(x,y,z,'-o'); hold on;
pptarget = plot3(0,0,0,'x','MarkerSize',15);
xl = xlabel('x (km)');
yl = ylabel('y (km)');
zl = zlabel('z (km)');
tt = title('3-D Trajectory of Inspector');
% ll = legend('No Control','Controlled Trajectory','Target');
ll = legend('Trajectory','Target');

% pp = plot(time,r,'-o');
% xl = xlabel('$t$','Interpreter','LaTeX');
% yl = ylabel('$r(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(zl,'FontSize',18);
set(tt,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
% set(p,'LineWidth',1.25);
set(pp,'LineWidth',1.25);
set(pptarget,'LineWidth',1.25);
grid on

%print -depsc2 orbitRaisingRadius.eps
print -dpng HCWmintToTarget_3D.png

% figure(3)
% pp = plot(time,theta,'-o');
% xl = xlabel('$t$','Interpreter','LaTeX');
% yl = ylabel('$\theta(t)$','Interpreter','LaTeX');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontSize',16,'FontName','Times');
% set(pp,'LineWidth',1.25);
% grid on
% %print -depsc2 orbitRaisingTheta.eps
% print -dpng orbitRaisingTheta.png
% 
% figure(4)
% pp = plot(time,vr,'-o');
% xl = xlabel('$t$','Interpreter','LaTeX');
% yl = ylabel('$v_r(t)$','Interpreter','LaTeX');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontSize',16,'FontName','Times');
% set(pp,'LineWidth',1.25);
% grid on
% %print -depsc2 orbitRaisingVr.eps
% print -dpng orbitRaisingVr.png
% 
% figure(5)
% pp = plot(time,vtheta,'-o');
% xl = xlabel('$t$','Interpreter','LaTeX');
% yl = ylabel('$v_\theta(t)$','Interpreter','LaTeX');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontSize',16,'FontName','Times');
% set(pp,'LineWidth',1.25,'MarkerSize',8);
% grid on
% %print -depsc2 orbitRaisingVtheta.eps
% print -dpng orbitRaisingVtheta.png

figure(6)
pp = plot(time,control,'-o');
hold on;
for i = 1:size(control,1)
    control_norm(i) = norm(control(i,:));
end
pp2 = plot(time,control_norm);
xl = xlabel('$t$ (min)','Interpreter','LaTeX');
yl = ylabel('Control Directions $u_1,u_2,u_3$','Interpreter','LaTeX');
tt = title('Control in $\hat{x}, \hat{y}, \hat{z}$ Directions','Interpreter','LaTeX');
ll = legend('$u_1$','$u_2$','$u_3$','Norm');
set(ll,'FontSize',18,'Interpreter','LaTeX');
set(tt,'FontSize',18);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
set(pp2,'LineWidth',1.25);
grid on
print -depsc2 HCWmintToTarget_Control.eps;
print -dpng HCWmintToTarget_Control.png;

alpha = atan2(u1,u2)*180/pi;
beta = atan2(u3,sqrt(u1.^2+u2.^2))*180/pi;

figure(7)
pp = plot(time,alpha,'-o'); hold on;
pp2 = plot(time,beta,'-o');
xl = xlabel('$t$ (min)','Interpreter','LaTeX');
% yl = ylabel('$\alpha(t)=\displaystyle\tan^{-1}\left(\frac{u_1(t)}{u_2(t)}\right)$','Interpreter','LaTeX');
yl = ylabel('$\alpha, \beta$ in degrees','Interpreter','LaTeX');
tt = title('Control Angles vs. Time');
ll = legend('$\alpha$','$\beta$');
set(ll,'FontSize',18,'Interpreter','LaTeX');
set(tt,'FontSize',18);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
set(gca,'FontSize',16);
set(pp,'LineWidth',1.25);
set(pp2,'LineWidth',1.25);
grid on
print -depsc2 HCWmintToTarget_ControlAngles.eps;
print -dpng HCWmintToTarget_ControlAngles.png;
