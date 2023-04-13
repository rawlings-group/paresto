%% Note: This file doesn't  make a data file and plot.
%% It just estimates parameter values that appear in the text of ch9.
%%
%% jbr, elh,  11/28/01
%%
%% A = cccDNA, B = rcDNA, C = env
%%
%% Reactions
%% 1: A -> A + B
%% 2: B -> A
%% 3: A -> A + C
%% 4: A -> 0 (degradation)
%% 5: C -> 0 (degradation)
%% 6: B + C -> 0 (secreted virus)
%%
%% reduced hbv model
%%
%% converted to use parest, jbr, 4/24/2007
%% converted to use paresto (casadi), jbr, 4/13/2018
%% updated to paresto structs for parameters, jbr, 6/1/2019
%% updated to revised paresto, jbr, 6/7/2020
%% pass constants with anonymous functions, not as parameters, jbr, 4/3/2021
%%

model=struct;
model.transcription = 'shooting';
model.print_level = 1;
%model.nlp_solver_options.ipopt.mumps_scaling = 0;
model.nlp_solver_options.ipopt.linear_solver = 'ma27';
%model.nlp_solver_options.ipopt.print_level = 5;

model.x = {'ca', 'cb', 'cc'};
model.p = {'k1', 'k2', 'k3', 'k4', 'k5', 'k6'};
model.d = {'ma', 'mb', 'mc'};

% Non-scalar dimensions
%model.dim.k = 6;

function rhs = hbv_rxs(t, y, p)
  kr = [ 10.^(p.k1) + p.k4, ...
	 10.^(p.k2) * p.k4, ...
	 p.k3, ...
	 p.k4, ...
	 p.k3 / 10.^p.k5, ...
	 10.^p.k6 ];
  rhs = {kr(2)*y.cb - kr(4)*y.ca, ...
         kr(1)*y.ca - kr(2)*y.cb - kr(6)*y.cb*y.cc, ...
         kr(3)*y.ca - kr(5)*y.cc - kr(6)*y.cb*y.cc};
end%function

model.ode = @hbv_rxs;

%% measurement weights
mweight = sqrt([1; 1e-2; 1e-4]);
w.a = mweight(1); w.b = mweight(2); w.c = mweight(3);

function retval = lsq_weighted(t, y, p, w)
  retval = {w.a*(y.ca-y.ma), w.b*(y.cb-y.mb), w.c*(y.cc-y.mc)};
endfunction

model.lsq = @(t, y, p) lsq_weighted(t, y, p, w);

%% Set the reaction rate constant vector kr; use log transformation
krac = [2; 0.025; 1000; 0.25; 1.9985; 7.5E-6];
thetaac = [log10(krac(1)-krac(4)); ...
           log10(krac(2)/krac(4)); ...
           krac(3); ...
	   krac(4); ...
	   log10(krac(3)/krac(5));
           log10(krac(6))];

for i = 1:numel(model.p)
  fn = model.p{i};
  p.(fn) = thetaac(i);
endfor

%% output times
tfinal = 100;
ndata = 51;
tdata = linspace(0, tfinal, ndata)';
model.tout = tdata;

pe = paresto(model);

%% Set the initial condition
small = 0;
x0_ac  = [1; small; small];

y_ac = pe.simulate(zeros(3, 1), x0_ac, p);

randn('seed',1);
R = diag([0.1 0.1 0.1])^2;
noise = sqrt(R)*randn(size(y_ac));
%% proportional error
y_noisy = y_ac .* (1 + noise);
y_noisy = max(y_noisy, 0);

%% initialize parameters and bounds
theta0 = p;
theta0.ca = x0_ac(1);
theta0.cb = x0_ac(2);
theta0.cc = x0_ac(3);
lb = theta0;
ub = theta0;

%% list of parameters to estimate
ind = [1, 2, 5, 6];
del = 1;

%% perturb initial guess and loosen bounds bounds on estimated parameters 
for i = 1:numel(ind)
  name = ['k' num2str(ind(i))];
  theta0.(name) = theta0.(name) + 0.75*del;
  lb.(name) = lb.(name) - del;
  ub.(name) = ub.(name) + del;
endfor

[est, y, p] = pe.optimize(y_noisy, theta0, lb, ub);

conf = pe.confidence(est, 0.95);

disp('Initial guess')
disp(theta0)

disp('Optimal parameters')
disp(est.theta)

disp('Confidence intervals')
disp(conf.bbox)

disp('True parameters')
disp(thetaac(ind))

%% optimal fit
subplot(3,1,1)
plot(tdata, y.ca, tdata, y_noisy(1,:), 'o')
axis ([0, 100, 0, 60]);

subplot(3,1,2)
plot(tdata, y.cb, tdata, y_noisy(2,:), 'o')
axis ([0, 100, 0, 600]);

subplot(3,1,3)
plot(tdata, y.cc, tdata, y_noisy(3,:), 'o')
axis ([0, 100, 0, 30000]);

