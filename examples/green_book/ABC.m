% [depends] ABC_data.dat
%%
%% Estimate parameters for the A->B->C model using paresto.m
%% jbr,  4/2018
%%

%%Model
model = struct;
%%model.print_level = 1;
model.transcription = "shooting";
model.nlp_solver_options.ipopt.linear_solver = 'ma27';
model.x = {'ca', 'cb', 'cc'};
model.p = {'k1', 'k2'};
%% measurement list
model.d = {'m_ca', 'm_cb', 'm_cc'};

function xdot = massbal_ode(t, x, p)
  r1 = p.k1*x.ca;
  r2 = p.k2*x.cb;
  xdot = {-r1, r1-r2, r2};
end%function

model.ode = @massbal_ode;

ncase = 1;
switch ncase
  case 1
    %% fit A,B,C
    model.lsq = @(t, y, p) {y.ca-y.m_ca, y.cb-y.m_cb, y.cc-y.m_cc};
    nmeas = [1;2;3];
  case 2
    %% fit A only
    model.lsq = @(t, y, p) {y.ca-y.m_ca};
    nmeas = [1];
  case 3
    %% fit B only
    model.lsq = @(t, y, p) {y.cb-y.m_cb};
    nmeas = [2];
  case 4
    %% fit C only
    model.lsq = @(t, y, p) {y.cc-y.m_cc};
    nmeas = [3];
  otherwise
    error ('ncase out of range: %f\n', ncase)
end%switch

tfinal = 6;
nplot = 75;
tplot = linspace(0, tfinal, nplot)';

%% For reference; actual parameters generating data (see ABC_data.m)
k1 = 2; k2 = 1;
ca0 = 1; cb0 = 0; cc0 = 0;

%% load measurements from a file
tabledat = load ('ABC_data.dat');
tmeas = tabledat(:,1);
ymeas = tabledat(:,2:4);

[tout,~,ic] = unique([tmeas; tplot]);
% index of times at which measurement is made
meas_ind = ic(1:numel(tmeas));
model.tout = tout;
% interploate measurement onto new grid
y_noisy = interp1(tmeas, ymeas, tout, 'previous');
y_noisy(isnan(y_noisy)) = 0.;
model.lsq_ind = meas_ind'; % only use index of measurement times in objective function

% Create a paresto instance
pe = paresto(model);

%% parameters; initial guesses and bounds;
thetaic.k1 = 0.5;
thetaic.k2 = 3;
thetaic.ca = ca0;
thetaic.cb = cb0;
thetaic.cc = cc0;

lb.k1 = 1e-4;
lb.k2 = 1e-4;
lb.ca = ca0;
lb.cb = cb0;
lb.cc = cc0;

ub.k1 = 10;
ub.k2 = 10;
ub.ca = ca0;
ub.cb = cb0;
ub.cc = cc0;

est_ind = 1:2;
%%est_ind = 1:5;

%% estimate the parameters
est = pe.optimize(y_noisy', thetaic, lb, ub);

% Also calculate confidence intervals with 95 % confidence
theta_conf = pe.confidence(est, est_ind, 0.95);

disp('Estimated Parameters and Bounding Box')
[est.theta(est_ind)  theta_conf]

%%plot the model fit to the noisy measurements

if (~ strcmp (getenv ('OMIT_PLOTS'), 'true')) %% PLOTTING
  figure(1)
  plot(model.tout, est.x(nmeas,:)', tmeas, ymeas(:,nmeas), 'o');
  if (ncase != 1)
    figure(2)
    plot(model.tout, est.x', tmeas, ymeas, 'o');
  endif
  %% TITLE
endif %% PLOTTING

tableest = [model.tout, est.x'];
save ABC.dat tableest tabledat
