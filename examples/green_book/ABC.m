% [depends] ABC_data.dat
%%
%% Estimate parameters for the A->B->C model using paresto.m
%% jbr,  4/2018
%%

%%Model
model = struct;
model.nlp_solver_options.ipopt.linear_solver = 'ma27';
model.x = {'ca', 'cb', 'cc'};
model.p = {'k1', 'k2'};
model.d = {'m_ca', 'm_cb', 'm_cc'};

function xdot = massbal_ode(t, x, p)
  r1 = p.k1*x.ca;
  r2 = p.k2*x.cb;
  xdot = {-r1, r1-r2, r2};
end%function

model.ode = @massbal_ode;

model.lsq = @(t, y, p) {y.ca-y.m_ca, y.cb-y.m_cb, y.cc-y.m_cc};

tfinal = 6;
nplot = 75;
tplot = linspace(0, tfinal, nplot)';

k1 = 2; k2 = 1;
p_ac = [k1; k2];

ca0 = 1; cb0 = 0; cc0 = 0;
x_ac = [ca0; cb0; cc0];

thetaic.k1 = 0.5;
thetaic.k2 = 3;
thetaic.ca = ca0
thetaic.cb = cb0;
thetaic.cc = cc0;

lb.k1 = 1e-4;
lb.k2 = 1e-4;

ub.k1 = 10;
ub.k2 = 10;

est_ind = 1:2;

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

%% estimate the parameters
est = pe.optimize(y_noisy', thetaic, lb, ub);

% Also calculate confidence intervals with 95 % confidence
theta_conf = pe.confidence(est, est_ind, 0.95);

disp('Estimated Parameters and Bounding Box')
[est.theta(est_ind)  theta_conf]

%%plot the model fit to the noisy measurements
figure(1);
plot(model.tout, est.x, tmeas, ymeas', 'o');

tableest = [model.tout, est.x'];
save ABC.dat tableest tabledat
