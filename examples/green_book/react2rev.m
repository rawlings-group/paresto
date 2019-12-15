%%
%%
%% react2rev.m, solves Exercise 9.11 of the 2nd edtion
%%
%% converted to paresto.m
%% Joel Andersson, jbr, 4/24/2018
%%
%%

% Model
model = struct;
model.transcription = 'shooting';
model.nlp_solver_options.ipopt.linear_solver = 'ma27';
%model.nlp_solver_options.ipopt.print_level = 0;
%model.nlp_solver_options.print_time = false;

model.x = {'ca', 'cb'};
model.p = {'k', 'n', 'm'};

% Measurements of the state enter the problem as data
model.d = {'m_ca', 'm_cb'};

% ODE right-hand-side
model.ode = @(t, y, p) {-p.k * y.ca^p.n * y.cb^p.m, ...
                        -p.k * y.ca^p.n * y.cb^p.m};

% Least squares objective function
model.lsq = @(t, y, p) {y.ca-y.m_ca, y.cb-y.m_cb};

% Options
load "react2rev_data.dat"
tout = react2rev_data(:,1)';
nts  = numel(tout);
ymeas = react2rev_data(:,2:5);

model.tout = tout;

% True parameter values (not known to optimizer)
kac   = 0.1;
nac   = 2;
mac   = 1;
p_ac = [kac; nac; mac];

parts = {"a", "b", "c"};

for npart = 1:numel(parts)

if (parts{npart} == "a")
  %% part a, using only first experiment
  use_first = true; use_second = false;
elseif(parts{npart} == "b")
  %% part b, using only second experiment
  use_first = false; use_second = true;
elseif (parts{npart} == "c")
  %% part c, using both experiments
  use_first = true; use_second = true;
endif

% Number of experiments
nsets = use_first + use_second;

model.nsets = nsets;

% Create a paresto instance
pe = paresto(model);

% Actual values for ca0 and cb0; load data sets into y_noisy
x_ac = NaN(numel(model.x), nsets);
y_noisy = NaN(numel(model.d),nts, nsets);
if use_first
  x_ac(:,1) = [2; 4];
  y_noisy(:,:,1) = ymeas(:,1:2)';
end
if use_second
  x_ac(:,nsets) = [5; 3];
  y_noisy(:,:,nsets) = ymeas(:,3:4)';
end

% Initial guess, upper and lower bounds for the estimated parameters
theta0 = struct;
theta0.k = 1.1*kac;
theta0.n = 1.1*nac;
theta0.m = 1.1*mac;
theta0.ca = NaN(1,1,nsets);
theta0.ca(:,:,1:nsets) = x_ac(1,:);
theta0.cb = NaN(1,1,nsets);
theta0.cb(:,:,1:nsets) = x_ac(2,:);

lbtheta = struct;
lbtheta.k = 0.5*kac;
lbtheta.n = 0.5*nac;
lbtheta.m = 0.5*mac;
lbtheta.ca = 0.5*theta0.ca;
lbtheta.cb = 0.5*theta0.cb;
ubtheta = struct;
ubtheta.k = 1.5*kac;
ubtheta.n = 1.5*nac;
ubtheta.m = 1.5*mac;
ubtheta.ca = 1.5*theta0.ca;
ubtheta.cb = 1.5*theta0.cb;

est_ind = 1:(3+ nsets*2);

% Estimate parameters
[est,y,p] = pe.optimize(y_noisy, theta0, lbtheta, ubtheta);

% Also calculate confidence intervals with 95 % confidence
theta_conf = pe.confidence(est, est_ind, 0.95);

disp('Estimated parameters and confidence intervals')
[est.theta(est_ind), theta_conf]

% Plot optimized trajectories
figure(npart)
clf
for i=1:nsets
  subplot(nsets,1,i)
  plot(tout, y.ca(:,:,i));
  hold on;
  plot(tout, y.cb(:,:,i));
  plot(tout, y_noisy(:,:,i), 'o');
  legend({'c_A (estimated)', 'c_B (estimated)', 'c_A (measured)', 'c_B (measured)'})
  title('Optimal fit')
  xlabel('t (min)')
  ylabel('c (mol/L)')
end

endfor
