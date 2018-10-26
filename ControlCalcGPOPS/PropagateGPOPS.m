function output = PropagateGPOPS(state0, statef, t0, tf, auxdata)

%--------------------------------------------------------------------------%
%--------------- Set Up Bounds on State, Control, and Time ----------------%
%--------------------------------------------------------------------------%

x0 = state0(1);
y0 = state0(2);
z0 = state0(3);
xdot0 = state0(4);
ydot0 = state0(5);
zdot0 = state0(6);

xf = statef(1);
yf = statef(2);
zf = statef(3);
xdotf = statef(4);
ydotf = statef(5);
zdotf = statef(6);


tfmin = tf; % 1 min (in seconds)
tfmax = tf; % 2 hr (in seconds)

slewrate_max = deg2rad(.5);

% x0 = -30*1000; % m % these distances need to be in bounds to make HCW applicable, -30, -60, 15 is baseline for initial conditions
% y0 = -60*1000; % m
% z0 = 15*1000; % m
% xdot0 = 0; % m/s
% ydot0 = 0; % m/s
% zdot0 = 0; % m/s
alpha0 = atan2(-y0,-x0); 
beta0 = atan2(-z0,sqrt(x0^2+y0^2));

% Hard and/or safety bounds
xmin = -200*1000; % these are all in (m)
xmax = 200*1000; 
ymin = -200*1000;
ymax = 200*1000;
zmin = -200*1000;
zmax = 200*1000;
xdotmin = -1*1000; % 1 km/s
xdotmax = 1*1000; % 1 km/s
ydotmin = -1*1000; % 1 km/s
ydotmax = 1*1000; % 1 km/s
zdotmin = -1*1000; % 1 km/s
zdotmax = 1*1000; % 1 km/s
alphamin = -6*pi;
alphamax = 6*pi;
betamin = -6*pi;
betamax = 6*pi;

% Control bounds - don't change from phase to phase
alphadotmin = -(slewrate_max);
alphadotmax = (slewrate_max);
betadotmin = -(slewrate_max);
betadotmax = (slewrate_max);
thmin = 0;
thmax = 1;

% Put values into GPOPS bounds structure
iphase = 1;

bounds.phase(iphase).initialtime.lower = t0; % fixed
bounds.phase(iphase).initialtime.upper = t0; % fixed
bounds.phase(iphase).finaltime.lower = tfmin; % lower bound
bounds.phase(iphase).finaltime.upper = tfmax; % upper bound

bounds.phase(iphase).initialstate.lower = [x0,y0,z0,xdot0,ydot0,zdot0,alpha0,beta0]; % fixed
bounds.phase(iphase).initialstate.upper = [x0,y0,z0,xdot0,ydot0,zdot0,alpha0,beta0]; % fixed
bounds.phase(iphase).state.lower = [xmin,ymin,zmin,xdotmin,ydotmin,zdotmin,alphamin,betamin]; % lower bound
bounds.phase(iphase).state.upper = [xmax,ymax,zmax,xdotmax,ydotmax,zdotmax,alphamax,betamax]; % upper bound
bounds.phase(iphase).finalstate.lower = [xf,yf,zf,xdotf,ydotf,zdotf, alphamax,betamax]; % lower bound
bounds.phase(iphase).finalstate.upper = [xf,yf,zf,xdotf,ydotf,zdotf, alphamax,betamax]; % upper bound

bounds.phase(iphase).control.lower = [alphadotmin,betadotmin,thmin]; % lower bound
bounds.phase(iphase).control.upper = [alphadotmax,betadotmax,thmax]; % upper bound

% bounds.phase(iphase).path.lower  = 1; % unit norm
% bounds.phase(iphase).path.upper  = 1; % unit norm
bounds.phase(iphase).integral.lower = [0]; 
bounds.phase(iphase).integral.upper = [100];

%--------------------------------------------------------------------------%
%------------------------- Set Up Initial Guess ---------------------------%
%--------------------------------------------------------------------------%
tGuess = [t0; tfmax]; % fixed initial time, guessing half of upper bound

xGuess = [x0;0];
yGuess = [y0;0];
zGuess = [z0;0];
xdotGuess = [xdot0;0]; 
ydotGuess = [ydot0;0]; 
zdotGuess = [zdot0;0]; % make this better? 
alphaGuess = [alpha0;alpha0+pi]; % EHHH, this could cause problems...
betaGuess = [beta0; beta0+pi]; % likewise

alphadotGuess = [0; 0]; % These aren't good guesses
betadotGuess = [0; 0];
thGuess = [1;1];

% Put values into GPOPS structure
guess.phase(iphase).time = tGuess;
guess.phase(iphase).state = [xGuess,yGuess,zGuess,xdotGuess,ydotGuess,zdotGuess,alphaGuess,betaGuess];
guess.phase(iphase).control = [alphadotGuess,betadotGuess,thGuess];
guess.phase(iphase).integral = 5;


%--------------------------------------------------------------------------%
%------------------------- Set Up Initial Mesh ----------------------------%
%--------------------------------------------------------------------------%
N = 10;
meshphase.colpoints = 4*ones(1,N);
meshphase.fraction   = ones(1,N)/N;

%--------------------------------------------------------------------------%
%-------------------------- Set Up for Solver -----------------------------%
%--------------------------------------------------------------------------%
setup.name = 'GPOPS_Mark';
setup.functions.continuous = @GPOPS_Continuous; % the dynamics (and path constraints?)
setup.functions.endpoint = @GPOPS_Endpoint;

setup.displaylevel = 2; % ?
setup.nlp.solver = 'ipopt';
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;

setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';

setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-4;
setup.mesh.phase = meshphase;
setup.scales.method = 'automatic-hybrid';

%--------------------------------------------------------------------------%
%-------------------- Solve Problem and Extract Solution ------------------%
%--------------------------------------------------------------------------%
output = gpops2(setup);


end
