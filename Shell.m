clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%% Shell Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scenario Parameters
Problem.TimeTotal = 60*60*24*2.5; % seconds
Problem.NumberMarks = 2;
Problem.Omega = 7.291e-5; % mean motion (rad/sec)
Problem.InitState = [30 30 0 0 0 0]; % [x,y,z,xdot,ydot,zdot];

% UTC and JD time at beginning of maneuver, i.e. at t_0
yr0 = 2017;
mo0 = 8;
day0 = 15;
hr0 = 18;
min0 = 7;
sec0 = 57.84;
Problem.Time.UTCO = [yr0,mo0,day0,hr0,min0,sec0]; 
Problem.Time.JD = JD(yr0,mo0,day0,hr0,min0,sec0);

%% Bounds
Problem.lb_init = 1;
Problem.ub_init = 60*60;
Problem.lb_loiter = 60*60*24;
Problem.lb_transfer = 60*30;
Problem.ub_loiter = 60*60*25;
Problem.ub_transfer = 60*60*24;

%% Optimization Parameters
Problem.GA.Popsize = 50;
Problem.GA.EliteCount = ceil(0.1*Problem.GA.Popsize);
Problem.GA.MaxGenerations = 500;
Problem.GA.UseParallel = false;
Problem.GA.CrossoverFraction = 0.9;

%% Define Marks
% ['NMC', ae, yd, xd, constraint]
% ['TD', size, yd, xd, constraint]
Problem.Mark.Info = [{'NMC',5,0,0,'sun'};
                     {'NMC',10,0,0,'sun'}];

%% Fitness Function
% 1 = HCWimpulsive, 2 = BVP
Problem.GA.FitnessFunction = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Unpack Shell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Problem = UnpackShell(Problem);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Run GA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Solution = RunGA(Problem);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Solution = PropagateSolution(Solution,Problem);

