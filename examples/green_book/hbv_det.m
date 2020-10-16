% [makes] hbv_det.dat

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

%% converted to use parest, jbr, 4/24/2007
%% converted to use paresto (casadi), jbr, 4/13/2018
%% updated to paresto structs for parameters, jbr, 5/24/2019
% updated to revised paresto, jbr, 6/7/2020


model=struct();
%model.print_level = 1;
%model.nlp_solver_options.ipopt.linear_solver = 'ma27';
model.nlp_solver_options.ipopt.mumps_scaling = 0;
%model.nlp_solver_options.ipopt.print_level = 0;
%model.nlp_solver_options.print_time = false;
%model.nlp_solver_options.ipopt.max_iter = 4000;
model.transcription = 'shooting';

model.x = {'ca', 'cb', 'cc'};
model.p = {'k1', 'k2', 'k3', 'k4', 'k5', 'k6', 'wa', 'wb', 'wc'};
model.d = {'ma', 'mb', 'mc'};

function rhs = hbv_rxs(t, y, p)
  kr = [ 10.^p.k1, ...
	 10.^p.k2, ...
	 10.^p.k3, ...
	 10.^p.k4, ...
	 10.^p.k5, ...
	 10.^p.k6 ];
  rhs = {kr(2)*y.cb - kr(4)*y.ca, ...
         kr(1)*y.ca - kr(2)*y.cb - kr(6)*y.cb*y.cc, ...
         kr(3)*y.ca - kr(5)*y.cc - kr(6)*y.cb*y.cc};
end%function

model.ode = @hbv_rxs;
model.lsq = @(t, y, p) {p.wa*(y.ca-y.ma), p.wb*(y.cb-y.mb), p.wc*(y.cc-y.mc)};

%% Set the reaction rate constant vector kr; use log transformation
krac = [2; 0.025; 1000; 0.25; 1.9985; 7.5E-6];
thetaac = log10(krac);

for i = 1:numel(thetaac)
  p.(['k' num2str(i)]) = thetaac(i);
endfor

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


%% initial guess for rate constants, page 544, edition 2.2.
logkrinit = [0.80, -1.13, 3.15, -0.77, -0.16, -5.46]';
p_init = [logkrinit; mweight];
%% outputs for the initial guess
y_init = pe.simulate(zeros(3, 1), x0_ac, p_init);


%% initialize all parameters
%% index of estimated rate constants
ind = [1, 2, 3, 4, 5, 6];

%%del = 1;
del = 0.1;

theta0 = p;
for i = 1:numel(ind)
  theta0.(['k' num2str(ind(i))]) = logkrinit(i);
endfor
theta0.ca = x0_ac(1);
theta0.cb = x0_ac(2);
theta0.cc = x0_ac(3);

lb = theta0;
ub = theta0;
%% loosen bounds bounds on estimated parameters

del = 1;

for i = 1: numel(ind)
  name = ['k' num2str(ind(i))];
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
disp(thetaac)

%% check eigenvalues/eigenvectors for badly determined parameter
%% directions

[u, lam] = eig(conf.H);

%%
%%  watch out, eigenvalues may be unordered; 
%%  run sort so i(1) will be the index of the smallest eigenvalue;
%%  still have a +/- sign ambiguity on u(:,i(1)), but thetap with either
%%  sign on u(:,i(1)) should give the same fit to the data

[s,i] = sort(diag(lam));
logkrp = cell2mat(struct2cell(est.theta)) + 0.5*u(:,i(1));

%% simulate with parameters perturbed in weak direction
pp = [logkrp; mweight];
yp = pe.simulate(zeros(3,1), x0_ac, pp);

store1 = [model.tout,  [y.ca;y.cb;y.cc]', y_init', yp'];
store2 = [model.tout, y_noisy'];
save hbv_det.dat store2 store1;

if (~ strcmp (getenv ('OMIT_PLOTS'), 'true')) %% PLOTTING
subplot (3, 2, 1);
plot (model.tout, y_noisy(1,:), '+', model.tout, [y.ca; y_init(1,:)]);
axis ([0, 100, 0, 60]);
%% TITLE hbv_det_cccdna

subplot (3, 2, 3);
plot (model.tout, y_noisy(2,:), '+', model.tout, [y.cb; y_init(2,:)]);
axis ([0, 100, 0, 600]);
%% TITLE hbv_det_rcdna

subplot (3, 2, 5);
plot (model.tout, y_noisy(3,:), '+', model.tout, [y.cc; y_init(3,:)]);
axis ([0, 100, 0, 30000]);
%% TITLE hbv_det_env

subplot (3, 2, 2);
plot (model.tout, y_noisy(1,:), '+', model.tout, [y.ca; yp(1,:)]);
axis ([0, 100, 0, 60]);
%% TITLE hbv_det_cccdnap

subplot (3, 2, 4);
plot (model.tout, y_noisy(2,:), '+', model.tout, [y.cb; yp(2,:)]);
axis ([0, 100, 0, 600]);
%% TITLE hbv_det_rcndap

subplot (3, 2, 6);
plot (model.tout, y_noisy(3,:), '+', model.tout, [y.cc; yp(3,:)]);
axis ([0, 100, 0, 30000]);
%% TITLE hbv_det_envp
endif %% PLOTTING
