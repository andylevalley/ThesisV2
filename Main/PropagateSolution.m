function Solution = PropagateSolution(Solution,Problem)

%% BVP Propagation (HCW) !!Needs Work!! 
% pull out relevent values from Problem structure
NumberMarks = Problem.NumberMarks;
TimeJD = Problem.Time.JD;
UTCO = Problem.Time.UTCO;
Omega = Problem.Omega;
dvar = Solution.dvar;

% Seperate chromosome
Order = dvar(1:NumberMarks);
TransferTimes = dvar(NumberMarks+1:end-NumberMarks);

% Targeting
InitState = Problem.InitState;
Marks = Problem.Mark.Info;
beta = dvar(end-NumberMarks+1:end);

count = NumberMarks;

TargetInfo = zeros(count,6);
for i = 1:count
    TargetInfo(i,1:6) = Target(Omega,Marks,beta,i,NumberMarks,Order,TransferTimes);
end

% Calculate Fitness
% initialize variables
CurrentState = InitState;
n = 1;
TrajState = [];
TrajControl = [];
TrajTime = [];
clock = 0;

for i = 1:length(Order)
    
    [x_drift,v_drift] = CWHPropagator(CurrentState(1:3)',CurrentState(4:6)',Omega,0:TransferTimes(n));
    CurrentState = [x_drift(1:3,end)',v_drift(1:3,end)'];
    
    TrajState = horzcat(TrajState, [x_drift(1:3,:);v_drift(1:3,:)]);
    TrajControl = horzcat(TrajControl,zeros(3,ceil(TransferTimes(n))));
    TrajTime = horzcat(TrajTime,clock:TransferTimes(n));
    
    solinit = bvpinit(linspace(0,ceil(TransferTimes(n+1)),150),[0 0 0 0 0 0 0 0 0 0 0 0]);
    options = bvpset('RelTol',1e-6);
    sol = bvp4c(@BVP_ode,@BVP_bc,solinit,options,CurrentState,TargetInfo(Order(i),:),Omega);
    t = [sol.x];
    x = [sol.y]; 
    u_x = x(10,:);
    u_y = x(11,:);
    u_z = x(12,:);
    CurrentState = [x(:,end)'];
    
    
    
    
    state = interp1(t,x',0:1:ceil(TransferTimes(n+1)));
    TrajState = horzcat(TrajState, state(:,1:6)');
    TrajControl = horzcat(TrajControl,state(:,10:12)');
    TrajTime = horzcat(TrajTime,clock+TransferTimes(n):1:clock+TransferTimes(n+1));
    
    n = n + 2;
    clock = ceil(TransferTimes(n))+ceil(TransferTimes(n+1))+clock;
    
end

[x_drift,v_drift] = CWHPropagator(CurrentState(1:3)',CurrentState(4:6)',Omega,0:TransferTimes(n));
CurrentState = [x_drift(1:3,end)',v_drift(1:3,end)'];

TrajState = horzcat(TrajState, [x_drift(1:3,:);v_drift(1:3,:)]);
TrajControl = horzcat(TrajControl,zeros(3,ceil(TransferTimes(n))));
TrajTime = horzcat(TrajTime,clock:TransferTimes(n));

ReturnState = InitState;

solinit = bvpinit(linspace(0,ceil(TransferTimes(n+1)),150),[0 0 0 0 0 0 0 0 0 0 0 0]);
options = bvpset('RelTol',1e-6);
sol = bvp4c(@BVP_ode,@BVP_bc,solinit,options,CurrentState,ReturnState,Omega);
t = [sol.x];
x = [sol.y];
u_x = x(10,:);
u_y = x(11,:);
u_z = x(12,:);
CurrentState = [x(:,end)'];

state = interp1(t,x',0:1:ceil(TransferTimes(n+1)));
TrajState = horzcat(TrajState, state(:,1:6)');
TrajControl = horzcat(TrajControl,state(:,10:12)');
TrajTime = horzcat(TrajTime,clock+TransferTimes(n):1:clock+TransferTimes(n+1));

f = trapz(0.5*(TrajControl(1,:).^2+TrajControl(2,:).^2));

alpha = atan2(TrajControl(1,:),TrajControl(2,:)).*180/pi;
beta = atan2(TrajControl(3,:),sqrt(TrajControl(1,:).^2+TrajControl(2,:).^2))*180/pi;

Solution.BVP.TrajState = TrajState;
Solution.BVP.TrajControl = [TrajControl; alpha;beta];
Solution.BVP.TrajTime = TrajTime;

%% Propagate GPOPS
% Calculate Fitness
% initialize variables

% auxdata.T  = 0.1405; % thrust
auxdata.a0 = .05; %.0686; % N/kg - initial acceleration .05 originally
% auxdata.m0 = 1; % initial mass
auxdata.c = 3.333e3; %3000; % m/s
% auxdata.dm = 0.0749; % mass flow rate
R_e = 6378.137; % km - Earth radius
alt = 35786; % km - GEO altitude
w = sqrt(398600.5/(R_e+alt)^3);
auxdata.w = w;

CurrentState = InitState;
n = 1;
TrajState = [];
TrajControl = [];
TrajTime = [];
clock = 0;

for i = 1:length(Order)
    
    [x_drift,v_drift] = CWHPropagator(CurrentState(1:3)',CurrentState(4:6)',Omega,0:TransferTimes(n));
    CurrentState = [x_drift(1:3,end)',v_drift(1:3,end)'];
    
    TrajState = horzcat(TrajState, [x_drift(1:3,:);v_drift(1:3,:)]);
    TrajControl = horzcat(TrajControl,zeros(5,ceil(TransferTimes(n))));
    TrajTime = horzcat(TrajTime,clock:TransferTimes(n));
    
    output = PropagateGPOPS(CurrentState(1:6),TargetInfo(Order(i),:),0,ceil(TransferTimes(n+1)),auxdata);
    solution = output.result.solution;
    
    x1    = solution.phase(1).state(:,1); % km
    y1    = solution.phase(1).state(:,2); % km
    z1    = solution.phase(1).state(:,3); % km
    xdot  = solution.phase(1).state(:,4); % km/s
    ydot  = solution.phase(1).state(:,5); % km/s
    zdot  = solution.phase(1).state(:,6); % km/s
    alpha = solution.phase(1).state(:,7); % rad
    beta  = solution.phase(1).state(:,8); % rad
    alphadot  = solution.phase(1).control(:,1); 
    betadot   = solution.phase(1).control(:,2); 
    thr       = solution.phase(1).control(:,3); 
    
    x = [x1';y1';z1';xdot';ydot';zdot';alpha';beta';alphadot';betadot';thr'];
    t = solution.phase(1).time; % NOW IN MINUTES
    
    CurrentState = [x(:,end)'];
    
    state = interp1(t',x',0:1:ceil(TransferTimes(n+1)));
    TrajState = horzcat(TrajState, state(:,1:6)');
    TrajControl = horzcat(TrajControl,state(:,7:11)');
    
    n = n + 2;
    clock = ceil(TransferTimes(n))+ceil(TransferTimes(n+1))+clock;
    
end

[x_drift,v_drift] = CWHPropagator(CurrentState(1:3)',CurrentState(4:6)',Omega,0:TransferTimes(n));
CurrentState = [x_drift(1:3,end)',v_drift(1:3,end)'];

TrajState = horzcat(TrajState, [x_drift(1:3,:);v_drift(1:3,:)]);
TrajControl = horzcat(TrajControl,zeros(5,ceil(TransferTimes(n))));
TrajTime = horzcat(TrajTime,clock:TransferTimes(n));


ReturnState = InitState;
output = PropagateGPOPS(CurrentState(1:6),ReturnState,0,ceil(TransferTimes(n+1)),auxdata);
solution = output.result.solution;

x1    = solution.phase(1).state(:,1); % km
y1    = solution.phase(1).state(:,2); % km
z1    = solution.phase(1).state(:,3); % km
xdot  = solution.phase(1).state(:,4); % km/s
ydot  = solution.phase(1).state(:,5); % km/s
zdot  = solution.phase(1).state(:,6); % km/s
alpha = solution.phase(1).state(:,7); % rad
beta  = solution.phase(1).state(:,8); % rad
alphadot  = solution.phase(1).control(:,1); 
betadot   = solution.phase(1).control(:,2); 
thr       = solution.phase(1).control(:,3); 

x = [x1';y1';z1';xdot';ydot';zdot';alpha';beta';alphadot';betadot';thr'];
t = solution.phase(1).time; % NOW IN MINUTES

CurrentState = [x(:,end)'];

state = interp1(t',x',0:1:ceil(TransferTimes(n+1)));
TrajState = horzcat(TrajState, state(:,1:6)');
TrajControl = horzcat(TrajControl,state(:,7:11)');

Solution.GPOPS.TrajState = TrajState;
Solution.GPOPS.TrajControl = TrajControl;
Solution.GPOPS.TrajTime = TrajTime;

%% Find Sun Angle throughout trajectory
Traj = Solution.GPOPS.TrajState;

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

TimeJD = Problem.Time.JD;
clock = 0;
n = 1;
i = 1;
Time(n) = clock;

nu = nu0 + w*clock;
[r,v] = coe2rvh (p,ecc,incl,Omega,argp,nu,nu,nu,lonper,mu);
RSO2Sun = sun2RSO(UTCO(1),UTCO(2),UTCO(3),UTCO(4),UTCO(5),UTCO(6),r,v);
theta(n,1) = acos(dot(Traj(1:3,i)',RSO2Sun(1:3,n)')/(norm(Traj(1:3,i)')*norm(RSO2Sun(1:3,n)')));


for i = 100:100:length(Traj)
    
    n = n + 1;
    clock = clock + 100;
    Time(n) = clock;

    nu = nu0 + w*clock;
    [r,v] = coe2rvh (p,ecc,incl,Omega,argp,nu,nu,nu,lonper,mu);
    RSO2Sun(1:3,n) = sun2RSO(UTCO(1),UTCO(2),UTCO(3),UTCO(4),UTCO(5),UTCO(6),r,v);
    
    theta(n,1) = acos(dot(-Traj(1:3,i)',RSO2Sun(1:3,n)')/(norm(-Traj(1:3,i)')*norm(RSO2Sun(1:3,n)')));
    
    if isnan(theta(n,1)) == 1
        theta(n,1) = 0;
    end

    
end

SunVec = interp1(Time,RSO2Sun',0:1:length(Traj))';
thetaInterp = interp1(Time,theta,0:1:length(Traj));

Solution.GPOPS.SunAngle = thetaInterp;
Solution.GPOPS.RSO2Sun = SunVec;


end