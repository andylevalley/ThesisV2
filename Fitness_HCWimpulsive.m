function f = Fitness_HCWimpulsive(dvar,Problem)

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
f = 0;

for i = 1:length(Order)
    
    tgt = dvar(i);
    
    [x_drift,v_drift] = CWHPropagator(CurrentState(1:3)',CurrentState(4:6)',Omega,TransferTimes(n));
    CurrentState = [x_drift',v_drift'];
    
    DeltaV = HCW_DeltaV(CurrentState(1:3)',TargetInfo(tgt,1:3)',...
        CurrentState(4:6)',TargetInfo(tgt,4:6)',TransferTimes(n+1),Omega);
    
    CurrentState = TargetInfo(tgt,:);
    
    f = f + norm(DeltaV);
    
    n = n + 2;
    
end

[x_drift,v_drift] = CWHPropagator(CurrentState(1:3)',CurrentState(4:6)',Omega,TransferTimes(n));
CurrentState = [x_drift',v_drift'];

ReturnState = [0 0 0 0 0 0];

DeltaV = HCW_DeltaV(CurrentState(1:3)',ReturnState(1:3)',...
    CurrentState(4:6)',ReturnState(4:6)',TransferTimes(n+1),Omega);

CurrentState = ReturnState;

f = f + norm(DeltaV);

end