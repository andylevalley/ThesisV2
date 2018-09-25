function [c ceq] = constraint(dvar,Problem)

%% pull out relevent values from Problem structure
NumberMarks = Problem.NumberMarks;
TimeJD = Problem.Time.JD;
TimeUTCO = Problem.Time.UTCO;
TimeTotal = Problem.TimeTotal;
Omega = Problem.Omega;


%% Seperate chromosome
Order = dvar(1:NumberMarks);
TransferTimes = dvar(NumberMarks+1:end);

test = 1:1:NumberMarks;
A = sum(Order == test,1);
B = ones(1,NumberMarks);


% c = [sum(dvar(num_objects+1:end))- (t_total); % change t_total and num_objects 
%     -(sum(order(:)==1)*sum(order(:)==2)*sum(order(:)==3)*sum(order(:)==4)*sum(order(:)==5))+1]; % this needs to change depending on number of objects 
% ceq = [];

c = [sum(dvar(NumberMarks+1:end))-(TimeTotal); % change t_total and num_objects 
    -(isequal(A,B))+1]; % this needs to change depending on number of objects 
ceq = [];
end