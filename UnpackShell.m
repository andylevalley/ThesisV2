function Problem = UnpackShell(Problem)


%% Derived GA Parameters (DO NOT EDIT)
NumberMarks  = Problem.NumberMarks;
lb_init = Problem.lb_init;
ub_init = Problem.ub_init;
lb_loiter = Problem.lb_loiter;
lb_transfer = Problem.lb_transfer;
ub_loiter = Problem.ub_loiter;
ub_transfer = Problem.ub_transfer;
lb_beta = Problem.lb_beta;
ub_beta = Problem.ub_beta;

numvars = NumberMarks + NumberMarks*2+2+2;
Popsize = Problem.GA.Popsize;

for i = 1:Popsize
    IntPop(i,:) = randperm(NumberMarks);
end

% IntPop = randi([1,4],Popsize,NumObjects);
IntPop = [IntPop, 1+(60*60-1).*rand(Popsize,NumberMarks*2+2)]; 

IntCon = [1:1:NumberMarks];

lb = [repmat(1,1,NumberMarks), repmat([lb_loiter,lb_transfer],1,(numvars-NumberMarks-2)/2),lb_beta,lb_beta];
ub = [repmat(NumberMarks,1,NumberMarks), repmat([ub_loiter,ub_transfer],1,(numvars-NumberMarks-2)/2),ub_beta,ub_beta];
lb(NumberMarks+1) = lb_init;
ub(NumberMarks+1) = ub_init;

Problem.GA.lb = lb;
Problem.GA.ub = ub;
Problem.GA.IntCon = IntCon;
Problem.GA.IntPop = IntPop;
Problem.GA.NumberVars = numvars;




end
