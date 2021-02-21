%%
%% jbr,  2/20/21
%%
%  simple dae to debug NaNs in sensitivities
%


daemodel = struct;
daemodel.transcription = 'shooting';
%daemodel.nlp_solver_options.ipopt.linear_solver = 'ma27';
daemodel.nlp_solver_options.ipopt.mumps_scaling = 0;
% set eps to zero for daeebraic model
daemodel.nlp_solver_options.sens_linsol_options.eps = 0;
daemodel.print_level = 1;
% use a dummy differential state to measure time
daemodel.x = {'x'};
daemodel.z = {'z'};
%daemodel.p = {'k', 'b'};
daemodel.p = {'k'};
daemodel.d = {'xmeas'};
tplot = linspace(0,1,3);


daemodel.tout = tplot;

daemodel.ode = @(t, y, p) {-p.k*y.z};
daemodel.alg = @(t, y, p) {y.z - y.x};
daemodel.lsq = @(t, y, p) {y.xmeas - y.x};

%% create measurements
x0 = 1;
k = 1;
b = 1;
z0 = b*x0;
randn('seed', 0);
meas = x0*exp(-k*tplot) + 0.1*randn(size(tplot));

p.k = 1;
%p.b = 2;

%% set up parameter initial guesses and bounds
theta0 = struct;
theta0.k = 2*k;
%theta0.b = b;
theta0.x = 0.5*x0;
theta0.z = z0;

lb = struct;
lb.k = sqrt(1E-3);
%lb.b = 10;
lb.x = 0;
lb.z = 0;

ub = struct;
ub.k = sqrt(5);
%ub.b = 200;
ub.x = 100;
ub.z = 100;

pe = paresto(daemodel);


%% check model
%% need a vector of parameters for simulate function
%pp = [p.k; p.b];
pp = [p.k];
ysim = pe.simulate(0, x0, pp, z0);

%% estimate the parameters

est = pe.optimize([meas], theta0, lb, ub);

conf =  pe.confidence(est, 0.95);

disp('Estimated parameters')
disp(est.theta)
disp('Bounding box intervals')
disp(conf.bbox)

