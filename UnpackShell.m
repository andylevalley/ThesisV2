function Problem = UnpackShell(Problem)


%% Derived GA Parameters (DO NOT EDIT)
NumObjects  = Problem.NumObjects;
lb_init = Problem.lb_init;
ub_init = Problem.ub_init;
lb_loiter = Problem.lb_loiter;
lb_transfer = Problem.lb_transfer;
ub_loiter = Problem.ub_loiter;
ub_transfer = Problem.ub_transfer;

numvars = NumObjects + NumObjects*2+2;

for i = 1:Popsize
    IntPop(i,:) = randperm(NumObjects);
end

% IntPop = randi([1,4],Popsize,NumObjects);
IntPop = [IntPop, 1+(60*60-1).*rand(Popsize,NumObjects*2+2)]; 

IntCon = [1:1:NumObjects];

lb = [repmat(1,1,NumObjects), repmat([lb_loiter,lb_transfer],1,(numvars-NumObjects)/2)];
ub = [repmat(NumObjects,1,NumObjects), repmat([ub_loiter,ub_transfer],1,(numvars-NumObjects)/2)];
lb(NumObjects+1) = lb_init;
ub(NumObjects+1) = ub_init;

Problem.GA.lb = lb;
Problem.GA.ub = ub;
Problem.GA.IntCon = IntCon;
Problem.GA.IntPop = IntPop;

%% Marker Properties






end
