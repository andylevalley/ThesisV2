function Solution = RunGA(Problem)

lb = Problem.GA.lb;
ub = Problem.GA.ub;
IntCon = Problem.GA.IntCon;
IntPop = Problem.GA.IntPop;

Popsize = Problem.GA.Popsize;
EliteCount = Problem.GA.EliteCount;
MaxGenerations = Problem.GA.MaxGenerations;
UseParallel = Problem.GA.UseParallel;
CrossoverFraction = Problem.GA.CrossoverFraction;
numvars = Problem.GA.NumberVars;

opts = optimoptions('ga','CrossoverFraction',CrossoverFraction,...
                    'PopulationSize',Popsize,'EliteCount',EliteCount,...
                    'MaxGenerations',MaxGenerations,'PlotFcn',@gaplotbestf,...
                    'UseParallel',UseParallel);
                
% Fitness and Penalty Function definition    
FitnessFunction = @(dvar) Fitness_BVP(dvar,Problem);
ConstraintFunction = @(dvar) constraint(dvar,Problem);

[dvar,fval,exitflag,output] = ga(FitnessFunction,numvars,[],[],[],[],...
    lb,ub,ConstraintFunction,IntCon,opts); 

Solution.dvar = dvar;
Solution.fval = fval;
Solution.exitflag = exitflag;
Solution.output = output;

end