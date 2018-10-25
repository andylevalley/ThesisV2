clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%% Shell Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scenario Parameters
Problem.TimeTotal = 60*60*24*3; % seconds
Problem.NumberMarks = 2;
Problem.Omega = 7.291e-5; % mean motion (rad/sec)
Problem.InitState = [0 0 0 0 0 0]; % [x,y,z,xdot,ydot,zdot];
Problem.TSP = 'yes';

% UTC and JD time at beginning of maneuver, i.e. at t_0
yr0 = 2017;
mo0 = 8;
day0 = 15;
hr0 = 18;
min0 = 7;
sec0 = 57.84;
Problem.Time.UTCO = [yr0,mo0,day0,hr0,min0,sec0]; 
Problem.Time.JD = JD(yr0,mo0,day0,hr0,min0,sec0);

%% Bounds
Problem.lb_init = 1;
Problem.ub_init = 60*60*24;
Problem.lb_loiter = 60*60*24;
Problem.lb_transfer = 60*30;
Problem.ub_loiter = 60*60*30;
Problem.ub_transfer = 60*60*24;
Problem.ub_beta = 2*pi;
Problem.lb_beta = 0;

%% Optimization Parameters
Problem.GA.Popsize = 200;
Problem.GA.EliteCount = ceil(0.1*Problem.GA.Popsize);
Problem.GA.MaxGenerations = 200;
Problem.GA.UseParallel = false;
Problem.GA.CrossoverFraction = 0.9;

%% Define Marks
% ['NMC', ae, xd, yd, zmax, upsilon, constraint, initState]
% ['TD', size, yd, xd, constraint]
ae = 5;
halfconeangle = deg2rad(20);
zmax1 = ae*tan(halfconeangle);
Problem.Mark.Info = [{'NMC',ae,0,0,0,0,'sun',[5 5 0 0 0 0]};
                     {'NMC',ae,0,0,zmax1,pi/2,'sun',[0 0 0 0 0 0]}];
                                  
Problem.Mark.Constraints.SunAngle = deg2rad(45);

                 
%% Define RSO or Virtual Chief
Problem.mu = 398600.5;
Problem.RSO.Parms.sma = 6378.137 + 35786; % semi-major axis of *circular* RSO or virtual RSO (km)
Problem.RSO.Parms.nu0 = deg2rad(0); % initial true anomaly of RSO or virtual RSO (deg)
Problem.RSO.Parms.ecc = 0; % eccentricity - this should be or stay very close to zero - we are using the HCW equations!
Problem.RSO.Parms.incl = deg2rad(0); % inclination (deg) - this combined with RAAN should be appropriate for use (see dissertation)
Problem.RSO.Parms.RAAN = deg2rad(0); % right ascension of the ascending node (deg) - see note above
Problem.RSO.Parms.argp = deg2rad(0); % argument of perigee (deg)
Problem.RSO.Parms.arglat = deg2rad(0); % for ci orbit
Problem.RSO.Parms.truelon = deg2rad(0); % for ce orbit
Problem.RSO.Parms.lonper = deg2rad(0); % for ee orbit

Problem.RSO.Parms.w = sqrt(Problem.mu/(Problem.RSO.Parms.sma)^3);
Problem.RSO.Parms.p = 2*pi/Problem.RSO.Parms.w;
Problem.RSO.Parms.p = Problem.RSO.Parms.sma*(1-Problem.RSO.Parms.ecc^2);

%% Fitness Function
% 1 = HCWimpulsive, 2 = BVP
Problem.GA.FitnessFunction = 1;

%% Unpack Shell 
Problem = UnpackShell(Problem);

%% Run GA 
Solution = RunGA(Problem);
close(gcf);

%% Propagate Solution 
Solution = PropagateSolution(Solution,Problem);

%% Plot solution

dvar = Solution.dvar;
NumberMarks = Problem.NumberMarks;
TransferTimes = dvar(NumberMarks+1:end-NumberMarks);
UTCO = Problem.Time.UTCO;

Time1 = ceil(sum(TransferTimes(1:2)));
Time2 = ceil(sum(TransferTimes(1:3)));
Time3 = ceil(sum(TransferTimes(1:4)));
Time4 = ceil(sum(TransferTimes(1:5)));
Time5 = ceil(sum(TransferTimes(1:6)));

figure(1)
plot3(Solution.GPOPS.TrajState(2,1:Time1),Solution.GPOPS.TrajState(1,1:Time1),Solution.GPOPS.TrajState(3,1:Time1),'r')
hold on
plot3(Solution.GPOPS.TrajState(2,Time1:Time2),Solution.GPOPS.TrajState(1,Time1:Time2),Solution.GPOPS.TrajState(3,Time1:Time2),'b')
hold on
plot3(Solution.GPOPS.TrajState(2,Time2:Time3),Solution.GPOPS.TrajState(1,Time2:Time3),Solution.GPOPS.TrajState(3,Time2:Time3),'r')
hold on
plot3(Solution.GPOPS.TrajState(2,Time3:Time4),Solution.GPOPS.TrajState(1,Time3:Time4),Solution.GPOPS.TrajState(3,Time3:Time4),'b')
hold on
plot3(Solution.GPOPS.TrajState(2,Time4:end),Solution.GPOPS.TrajState(1,Time4:end),Solution.GPOPS.TrajState(3,Time4:end),'r')
xlabel('y');
ylabel('x');
zlabel('z');
hold on

scatter3(Solution.GPOPS.TrajState(2,1),Solution.GPOPS.TrajState(1,1),Solution.GPOPS.TrajState(3,1),'MarkerFaceColor','g','MarkerEdgeColor','k')
hold on
scatter3(Solution.GPOPS.TrajState(2,Time1),Solution.GPOPS.TrajState(1,Time1),Solution.GPOPS.TrajState(3,Time1),'MarkerFaceColor','r','MarkerEdgeColor','k')
hold on
scatter3(Solution.GPOPS.TrajState(2,Time2),Solution.GPOPS.TrajState(1,Time2),Solution.GPOPS.TrajState(3,Time2),'MarkerFaceColor','g','MarkerEdgeColor','k')
hold on
scatter3(Solution.GPOPS.TrajState(2,Time3),Solution.GPOPS.TrajState(1,Time3),Solution.GPOPS.TrajState(3,Time3),'MarkerFaceColor','r','MarkerEdgeColor','k')
hold on
scatter3(Solution.GPOPS.TrajState(2,Time4),Solution.GPOPS.TrajState(1,Time4),Solution.GPOPS.TrajState(3,Time4),'MarkerFaceColor','g','MarkerEdgeColor','k')
hold on

% SunVec = Solution.GPOPS.RSO2Sun;
% q = quiver3(-SunVec(2,Time1)*5,-SunVec(1,Time1)*5,-SunVec(3,Time1)*5,SunVec(2,Time1)*5,SunVec(1,Time1)*5,SunVec(3,Time1)*5);
% q.Color = 'r';
% q = quiver3(0,0,0,SunVec(1,Time2),SunVec(2,Time1),SunVec(3,Time2));
% q.Color = 'r';
% 
% q = quiver3(0,0,0,SunVec(1,Time3),SunVec(2,Time3),SunVec(3,Time3));
% q.Color = 'g';
% q = quiver3(0,0,0,SunVec(1,Time4),SunVec(2,Time4),SunVec(3,Time4));
% q.Color = 'g';

figure(2)
SunAngle = Solution.GPOPS.SunAngle;
SA = plot(rad2deg(SunAngle),'-r');
SA.LineWidth = 2;
axis tight
title('Sun Angle Throughout Trajectory')
ylabel('Angle (degs)')
xlabel('Time (hrs)')

ticks = 0:3600*5:length(Solution.GPOPS.SunAngle);
for i = 1:length(ticks)
    labels(i) = ticks(i)/3600;
end
xticklabels(labels)
xticks(ticks)

x_points = [Time1, Time1, Time2, Time2];  
y_points = [0, 180, 180, 0];
color = [0, 0, 1];
hold on;
a = fill(x_points, y_points, color);
a.FaceAlpha = 0.1;

hold on
x_points = [Time3, Time3, Time4, Time4];  
y_points = [0, 180, 180, 0];
color = [0, 0, 1];
a = fill(x_points, y_points, color);
a.FaceAlpha = 0.3;

[~,h_legend] = legend('Sun Angle','NMC 1','NMC 2');
PatchInLegend = findobj(h_legend, 'type', 'patch');
set(PatchInLegend(1), 'FaceAlpha', 0.1);
set(PatchInLegend(2), 'FaceAlpha', 0.3);

ticks = 0:3600*5:length(Solution.GPOPS.SunAngle);
for i = 1:length(ticks)
    labels(i) = ticks(i)/3600;
end
xticklabels(labels)
xticks(ticks)


figure(3)
ylim auto
Traj = Solution.GPOPS.TrajState;
time = 1:1:length(Traj);
xplot = plot(time,Traj(1,:),'-.',time,Traj(2,:),'-.',time,Traj(3,:),'-.');
xplot(1).LineWidth = 2;
xplot(2).LineWidth = 2;
xplot(3).LineWidth = 2;


ticks = 0:3600*5:length(Traj);

for i = 1:length(ticks)
    labels(i) = ticks(i)/3600;
end

xticklabels(labels)
xticks(ticks)
xlabel('Time (hours)')
ylabel('km')
axis tight

x_points = [Time1, Time1, Time2, Time2];  
y_points = [-6, 6, 6, -6];
color = [0, 0, 1];
hold on;
a = fill(x_points, y_points, color);
a.FaceAlpha = 0.1;

hold on
x_points = [Time3, Time3, Time4, Time4];  
y_points = [-6, 6, 6, -6];
color = [0, 0, 1];
a = fill(x_points, y_points, color);
a.FaceAlpha = 0.3;

[~,h_legend] = legend('x','y','z','NMC 1','NMC 2');
PatchInLegend = findobj(h_legend, 'type', 'patch');
set(PatchInLegend(1), 'FaceAlpha', 0.1);
set(PatchInLegend(2), 'FaceAlpha', 0.3);


