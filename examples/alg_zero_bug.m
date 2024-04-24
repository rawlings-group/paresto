%% This file has a parameter, p.k1, that has been removed from the model
%% generates a casadi failure in the call to the sensitivity solver
%% should just generate zero derivatives with respect to this parameter
%%

%%
%% Simple algebraic model to investigate zero derivatives
%%
%% xA = k1; xB = k2

%%Model
model = struct;
model.print_level = 1;
model.transcription = "shooting";
%model.nlp_solver_options.ipopt.linear_solver = 'ma27';
model.nlp_solver_options.ipopt.mumps_scaling = 0;
%% set eps to zero for algebraic model
algmodel.nlp_solver_options.sens_linsol_options.eps = 0;
model.z = {'a', 'b'};
model.p = {'k1', 'k2'};
%% measurement list
model.d = {'m_a', 'm_b'};

model.alg = @(t, y, p) { y.a - p.k1, y.b - p.k2 };

model.lsq = @(t, y, p) { y.a - y.m_a, p.k2 - y.m_b };
%model.lsq = @(t, y, p) {p.k1^2 - y.m_a, p.k2^2 - y.m_b };

%% create two time point measurements
%%tout = [0, 1];
tout = [0];
model.tout = tout;

k1 = 1;
k2 = 2;

%ymeas = [ k1^2 + 0.5, k1^2 - 0.25; sqrt(k2) + 1, sqrt(k2) - 0.5 ];
%%ymeas = [ k1 + 1, k1 - 0.5;  k2 + 2, k2 - 1 ];
ymeas = [ k1 + 1;  k2 + 2 ];

% Create a paresto instance
pe = paresto(model);

%% parameters; initial guesses and bounds;
thetaic.k1 = 0;
thetaic.k2 = 0;
thetaic.a = 0;
thetaic.b = 0;

lb = thetaic;
ub = thetaic;
%% loosen algebraic IC and estimate parameters
lb.k1 = -Inf;
lb.k2 = -Inf;
lb.a = -Inf;
lb.b = -Inf;

ub.k1 = Inf;
ub.k2 = Inf;
ub.a = Inf;
ub.b = Inf;

%% estimate the parameters
est = pe.optimize(ymeas, thetaic, lb, ub);

% Also calculate confidence intervals with 95% confidence

conf = pe.confidence(est, 0.95);

disp('Estimated parameters')
disp(est.theta)
disp('Bounding box intervals')
disp(conf.bbox)


%%plot the model fit to the noisy measurements

