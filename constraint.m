function [c ceq] = constraint(dvar,Problem)

%% pull out relevent values from Problem structure
NumberMarks = Problem.NumberMarks;
TimeJD = Problem.Time.JD;
TimeUTCO = Problem.Time.UTCO;
TimeTotal = Problem.TimeTotal;
Omega = Problem.Omega;

r_RSO = [42,164;0;0];
v_RSO = [0;sqrt(mu/r_RSO);0];

%% Seperate chromosome
Order = dvar(1:NumberMarks);
TransferTimes = dvar(NumberMarks+1:end);

test = 1:1:NumberMarks;
A = sum(Order == test,1);
B = ones(1,NumberMarks);


% c = [sum(dvar(num_objects+1:end))- (t_total); % change t_total and num_objects 
%     -(sum(order(:)==1)*sum(order(:)==2)*sum(order(:)==3)*sum(order(:)==4)*sum(order(:)==5))+1]; % this needs to change depending on number of objects 
% ceq = [];

%% Time and visit constraints
c = [sum(TransferTimes)-(TimeTotal); % change t_total and num_objects 
    -(isequal(A,B))+1]; % this needs to change depending on number of objects 

%% Sun constraint
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
clock = TransferTimes(1);
n = 3;

for i = 1:NumberMarks
    
    tgt = dvar(i);
    secs = UTCO(6) + clock;
    
    Sun2RSO_Start = sun2RSO(UTCO(1),UTCO(2),UTCO(3),UTCO(4),UTCO(5),secs,r_RSO,v_RSO);
    
    StartState = TargetInfo(tgt,1:6);
    
    [x_drift,v_drift] = CWHPropagator(StartState(1:3)',StartState(4:6)',Omega,TransferTimes(n));
    EndState = [x_drift',v_drift'];
    
    secs = secs + TransferTimes(n);
    
    Sun2RSO_End = sun2RSO(UTCO(1),UTCO(2),UTCO(3),UTCO(4),UTCO(5),secs,r_RSO,v_RSO);
    
    
    clock = clock + TrasnferTimes(n-1) + TranferTimes(n);
    n = n + 2;
    
    
    
        


ceq = [];
end