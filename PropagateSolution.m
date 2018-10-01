function Solution = PropagateSolution(Solution,Problem)

%% BVP Propagation (HCW) !!Needs Work!! 
% pull out relevent values from Problem structure
NumberMarks = Problem.NumberMarks;
TimeJD = Problem.Time.JD;
TimeUTCO = Problem.Time.UTCO;
Omega = Problem.Omega;
dvar = Solution.dvar;

% Seperate chromosome
Order = dvar(1:NumberMarks);
TransferTimes = dvar(NumberMarks+1:end-NumberMarks);

% Targeting
InitState = Problem.InitState;
Marks = Problem.Mark.Info;
beta = dvar(end-NumberMarks+1:end);

count = size(Marks,1);

TargetInfo = zeros(count,6);
for i = 1:count
    TargetInfo(i,1:6) = Target(Omega,Marks,beta,i);
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

ReturnState = [0 0 0 0 0 0];

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


end