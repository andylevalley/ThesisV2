function [c ceq] = constraint(dvar,t_total,num_objects)

order = dvar(1:num_objects);
test = 1:1:num_objects;
A = sum(order(:) == test,1);
B = ones(1,num_objects);


% c = [sum(dvar(num_objects+1:end))- (t_total); % change t_total and num_objects 
%     -(sum(order(:)==1)*sum(order(:)==2)*sum(order(:)==3)*sum(order(:)==4)*sum(order(:)==5))+1]; % this needs to change depending on number of objects 
% ceq = [];

c = [sum(dvar(num_objects+1:end))- (t_total); % change t_total and num_objects 
    -(isequal(A,B))+1]; % this needs to change depending on number of objects 
ceq = [];
end