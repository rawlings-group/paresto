%%
%% jbr,  2/20/21
%%
%  simple dae to debug NaNs in sensitivities
%  model:       dx/dt = -k z;  z = x; x(0)=x_0;
%  measurement: x(t).
%

daemodel = struct;
%% daemodel.transcription = 'simultaneous';
%% daemodel.ord = 1;
daemodel.transcription = 'shooting';
daemodel.nlp_solver_options.ipopt.linear_solver = 'ma27';
%daemodel.nlp_solver_options.ipopt.mumps_scaling = 0;
% set eps to zero for daeebraic model
%%daemodel.nlp_solver_options.sens_linsol_options.eps = 0;
daemodel.print_level = 0;
%
daemodel.x = {'x'};
daemodel.z = {'z'};
%daemodel.p = {'k', 'b'};
daemodel.p = {'k'};
daemodel.d = {'xmeas'};
tplot = linspace(0,1,3);

%% get information on the multipliers
%%daemodel.nlp_solver_options.print_in = true;
%%daemodel.nlp_solver_options.print_out = true;
daemodel.nlp_solver_options.sens_linsol = 'csparse';

daemodel.tout = tplot;

daemodel.ode = @(t, y, p) {-p.k*y.x};
%daemodel.alg = @(t, y, p) {y.z - p.b*y.x};
daemodel.alg = @(t, y, p) {y.z - y.x};

%% NaNs in est if use the following without regularization term;
daemodel.lsq = @(t, y, p) {y.xmeas - y.x};
%%No NaNs when using regularization term
%%daemodel.lsq = @(t, y, p) {y.xmeas - y.x, 1e-10*y.z};
%% The following also works without generating NaNs
%%daemodel.lsq = @(t, y, p) {y.xmeas - y.z};

%% create measurements
x0 = 1;
k = 1;
b = 1;
z0 = b*x0;
randn('seed', 0);
meas = x0*exp(-k*tplot) + 0.1*randn(size(tplot));

p.k = k;
%p.b = b;

%% set up parameter initial guesses and bounds
theta0 = struct;
theta0.k = 2*k;
%theta0.b = b;
theta0.x = 1.1*x0;
theta0.z = z0;

lb = struct;
lb.k = sqrt(1E-3);
%lb.b = 0;
lb.x = theta0.x;

ub = struct;
ub.k = sqrt(5);
%ub.b = 200;
ub.x = theta0.x;

pe = paresto(daemodel);


%% check model
ysim = pe.simulate(0, x0, p, z0);

%% estimate the parameters

est = pe.optimize(meas, theta0, lb, ub);

conf =  pe.confidence(est, 0.95);

disp('Estimated parameters')
disp(est.theta)
disp('Bounding box intervals')
disp(conf.bbox)
