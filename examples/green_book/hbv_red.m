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
%%

more off

model=struct;
model.nlp_solver_options.ipopt.linear_solver = 'ma27';
model.nlp_solver_options.ipopt.print_level = 5;

model.x = {'ca', 'cb', 'cc'};
model.p = {'k', 'wa', 'wb', 'wc'};
model.d = {'ma', 'mb', 'mc'};

% Non-scalar dimensions
model.dim.k = 6;

function rhs = hbv_rxs(t, y, p)
  kr = [ 10.^(p.k(1)) + p.k(4), ...
	 10.^(p.k(2)) * p.k(4), ...
	 p.k(3), ...
	 p.k(4), ...
	 p.k(3) / 10.^p.k(5), ...
	 10.^p.k(6) ];
  rhs = {kr(2)*y.cb - kr(4)*y.ca, ...
         kr(1)*y.ca - kr(2)*y.cb - kr(6)*y.cb*y.cc, ...
         kr(3)*y.ca - kr(5)*y.cc - kr(6)*y.cb*y.cc};
end%function

model.ode = @hbv_rxs;
model.lsq = @(t, y, p) {p.wa*(y.ca-y.ma), p.wb*(y.cb-y.mb), p.wc*(y.cc-y.mc)};

%% Set the reaction rate constant vector kr; use log transformation
krac = [2; 0.025; 1000; 0.25; 1.9985; 7.5E-6];
thetaac = [log10(krac(1)-krac(4)); ...
           log10(krac(2)/krac(4)); ...
           krac(3); ...
	   krac(4); ...
	   log10(krac(3)/krac(5));
           log10(krac(6))];
p.k = thetaac;

%% output times
tfinal = 100;
ndata = 51;
tdata = linspace(0, tfinal, ndata)';
model.tout = tdata;

%% measurement weights
mweight = sqrt([1; 1e-2; 1e-4]);
p.wa = mweight(1); p.wb = mweight(2); p.wc = mweight(3);
p_ac = [thetaac; mweight];

pe = paresto(model);

%% Set the initial condition
small = 0;
x0_ac  = [1; small; small];

y_ac = pe.simulate(zeros(3, 1), x0_ac, p_ac);

randn('seed',1);
R = diag([0.1 0.1 0.1])^2;
noise = sqrt(R)*randn(size(y_ac));
%% proportional error
y_noisy = y_ac .* (1 + noise);
y_noisy = max(y_noisy, 0);

%% initialize all parameters
%% index of estimated parameters
ind = [1, 2, 5, 6];
del = 1;

theta0 = p;
theta0.k(ind) = theta0.k(ind) + 0.75*del;
theta0.ca = x0_ac(1);
theta0.cb = x0_ac(2);
theta0.cc = x0_ac(3);

%% loosen bounds bounds on estimated parameters
lb = struct();
lb.k = theta0.k;
lb.k(ind) = lb.k(ind) - del;
lb.wa = theta0.wa; lb.wb = theta0.wb; lb.wc = theta0.wc;

ub = struct();
ub.k = theta0.k;
ub.k(ind) = ub.k(ind) + del;
ub.wa = theta0.wa; ub.wb = theta0.wb; ub.wc = theta0.wc;

[est, y, p] = pe.optimize(y_noisy, theta0, lb, ub);

theta_conf = pe.confidence(est, ind, 0.95);

disp('Initial guess, optimal parameters, confidence intervals, true parameters')
[theta0.k(ind), est.theta(ind), theta_conf, thetaac(ind)]

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

