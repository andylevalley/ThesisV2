function [c ceq] = constraint(dvar,Problem)

%% pull out relevent values from Problem structure
NumberMarks = Problem.NumberMarks;
TimeJD = Problem.Time.JD;
UTCO = Problem.Time.UTCO;
TimeTotal = Problem.TimeTotal;
Omega = Problem.Omega;
SunAngleCon = Problem.Mark.Constraints.SunAngle;

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

%% Seperate chromosome
Order = dvar(1:NumberMarks);
TransferTimes = dvar(NumberMarks+1:end-2);

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
clock = sum(TransferTimes(1:2));
n = 3;


for i = 1:NumberMarks
    
    tgt = dvar(i);

    [RSO2Sun] = SunVecRSO(clock);
    StartState = TargetInfo(tgt,1:6);
    
    thetaStart = acos(dot(-StartState(1:3),RSO2Sun')/(norm(-StartState(1:3))*norm(RSO2Sun')));
    
    if strcmp(Marks{i,end},'sun') == 1
        c(end+1:end+2,1) = [-SunAngleCon + thetaStart];
        
    end
    
    clock = clock + TransferTimes(n) + TransferTimes(n+1);
    n = n + 2;
    
end
    
    
    
        


ceq = [];
end