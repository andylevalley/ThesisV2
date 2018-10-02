clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%% Shell Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scenario Parameters
Problem.TimeTotal = 60*60*24*4; % seconds
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
Problem.GA.MaxGenerations = 500;
Problem.GA.UseParallel = false;
Problem.GA.CrossoverFraction = 0.9;

%% Define Marks
% ['NMC', ae, yd, xd, zmax, upsilon, constraint]
% ['TD', size, yd, xd, constraint]
ae = 5;
halfconeangle = deg2rad(20);
zmax1 = ae*tan(halfconeangle);
Problem.Mark.Info = [{'NMC',ae,0,0,0,0,'sun'};
                     {'NMC',ae,0,0,zmax1,-pi/2,'sun'}];
                 
% Problem.Mark.Info = [{'NMC',ae,0,0,0,0,'sun'};
%                     {'NMC',ae+5,0,0,0,0,'sun'}];

                 
Problem.Mark.Constraints.SunAngle = deg2rad(30);

                 
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

% r_RSO = [42164;0;0];
% v_RSO = [0;3.0747;0];
% Sun2RSO1 = sun2RSO(UTCO(1),UTCO(2),UTCO(3),UTCO(4),UTCO(5),UTCO(6),r_RSO,v_RSO)';
% Sun2RSO2 = sun2RSO(UTCO(1),UTCO(2),UTCO(3),UTCO(4),UTCO(5),UTCO(6),r_RSO,v_RSO)';

% pull out relevent values from Problem structure
NumberMarks = Problem.NumberMarks;
TimeJD = Problem.Time.JD;
UTCO = Problem.Time.UTCO;
TimeTotal = Problem.TimeTotal;
Omega = Problem.Omega;

sma = Problem.RSO.Parms.sma; % semi-major axis of *circular* RSO or virtual RSO (km)
nu0 = Problem.RSO.Parms.nu0; % initial true anomaly of RSO or virtual RSO (deg)
ecc = Problem.RSO.Parms.ecc; % eccentricity - this should be or stay very close to zero - we are using the HCW equations!
incl = Problem.RSO.Parms.incl; % inclination (deg) - this combined with RAAN should be appropriate for use (see dissertation)
RAAN = Problem.RSO.Parms.RAAN; % right ascension of the ascending node (deg) - see note above
argp = Problem.RSO.Parms.argp; % argument of perigee (deg)
arglat = Problem.RSO.Parms.arglat; % for ci orbit
truelon = Problem.RSO.Parms.truelon; % for ce orbit
lonper = Problem.RSO.Parms.lonper; % for ee orbit
w = Problem.RSO.Parms.w;
p = Problem.RSO.Parms.p;
mu = Problem.mu;

% Seperate chromosome
Order = dvar(1:NumberMarks);
TransferTimes = dvar(NumberMarks+1:end-2);


% Sun constraint
InitState = Problem.InitState;
Marks = Problem.Mark.Info;
beta = dvar(end-NumberMarks+1:end);

count = size(Marks,1);

% Targets states defined by user
TargetInfo = zeros(count,6);
for i = 1:count
    TargetInfo(i,1:6) = Target(Omega,Marks,beta,i);
end

% Starting after initial wait and transfer
clock = sum(TransferTimes(1:2));
n = 3;


TimeJD = TimeJD + clock/3600/24; % clock is in seconds
UTCO = JDtoGregorianDate(TimeJD);
nu_current = nu0 + w*clock;

for i = 1:NumberMarks
    
    tgt = dvar(i);

    [r_RSO,v_RSO] = coe2rvh(p,ecc,incl,Omega,argp,nu_current,arglat,truelon,lonper,mu);
    Sun2RSO(i,:) = sun2RSO(UTCO(1),UTCO(2),UTCO(3),UTCO(4),UTCO(5),UTCO(6),r_RSO,v_RSO)';
    StartState = TargetInfo(tgt,1:6);
    
    SunAngleStart(i,:) = acos(dot(StartState(1:3),Sun2RSO(i,:))/(norm(StartState(1:3))*norm(Sun2RSO(i,:))));
    
    clock = clock + TransferTimes(n) + TransferTimes(n+1);
    TimeJD = TimeJD + clock/3600/24; % clock is in seconds
    UTCO = JDtoGregorianDate(TimeJD); 
    nu_current = nu_current + w*clock;
    n = n + 2;
    
end

q = quiver3(0,0,0,Sun2RSO(1,2),Sun2RSO(1,1),Sun2RSO(1,3));
q.Color = 'r';
hold on
q = quiver3(0,0,0,Sun2RSO(2,2),Sun2RSO(2,1),Sun2RSO(2,3));
q.Color = 'r';
hold off
