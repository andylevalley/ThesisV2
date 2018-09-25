function [f] = Fitness_BVP(dvar,Problem)

%% pull out relevent values from Problem structure
NumberMarks = Problem.NumberMarks;
TimeJD = Problem.Time.JD;
TimeUTCO = Problem.Time.UTCO;
Omega = Problem.Omega;

%% Seperate chromosome
Order = dvar(1:NumberMarks);
TransferTimes = dvar(NumberMarks+1:end);

%% Targeting
InitState = Problem.InitState;
Marks = Problem.Mark.Info;
beta = pi;

count = size(Marks,1);

TargetInfo = zeros(count,6);
for i = 1:count
    TargetInfo(i,1:6) = Target(Omega,Marks,beta,i);
end

%% Calculate Fitness

% initialize variables
CurrentState = InitState;
n = 1;
TrajState = [];
TrajControl = [];
TrajTime = [];

for i = 1:length(Order)
    
    [x_drift,v_drift] = CWHPropagator(CurrentState(1:3)',CurrentState(4:6)',Omega,TransferTimes(n));
    CurrentState = [x_drift',v_drift'];
    
    solinit = bvpinit(linspace(0,TransferTimes(n+1),150),[0 0 0 0 0 0 0 0]);
    options = bvpset('RelTol',1e-6);
    sol = bvp4c(@BVP_ode,@BVP_bc,solinit,options,CurrentState,TargetInfo(i,:),Omega);
    t = [sol.x];
    x = [sol.y]; 
    u_x = -x(6,:);
    u_y = -x(8,:);
    CurrentState = [x(1,end) x(3,end) 0 x(2,end) x(4,end) 0];
    
    TrajState = horzcat(TrajState,x);
    TrajControl = horzcat(TrajControl,[u_x;u_y]);
    TrajTime = horzcat(TrajTime,t);
    
    n = n + 2;
    
end

[x_drift,v_drift] = CWHPropagator(CurrentState(1:3)',CurrentState(4:6)',Omega,TransferTimes(n));
CurrentState = [x_drift',v_drift'];

ReturnState = [0 0 0 0 0 0];

solinit = bvpinit(linspace(0,TransferTimes(n+1),150),[0 0 0 0 0 0 0 0]);
options = bvpset('RelTol',1e-6);
sol = bvp4c(@BVP_ode,@BVP_bc,solinit,options,CurrentState,ReturnState,Omega);
t = [sol.x];
x = [sol.y];
u_x = -x(6,:);
u_y = -x(8,:);
CurrentState = [x(1,end) x(3,end) 0 x(2,end) x(4,end) 0];

TrajState = horzcat(TrajState,x);
TrajControl = horzcat(TrajControl,[u_x;u_y]);
TrajTime = horzcat(TrajTime,t);

f = trapz(0.5*(TrajControl(1,:).^2+TrajControl(2,:).^2));

end
    
    


