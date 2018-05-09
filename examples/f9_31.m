% We have the ODE
% dot(VR) = Qf
% dot(nA) = -k1*nA*nB/VR
% dot(nB) = Qf*cBf - nB*(k1*nA + k2*nC)/VR
% dot(nC) = nB*(k1*nA - k2*nC)/VR
% dot(nD) = k2*nC*nB/VR
%
%
% with:
% States: x = [VR, nA, nB, nC, nD]
% Qf: Volumetric flowrate of base
% cBf: Feed concentration of B
%
% Initial conditions: x(0) = [VR0, nA0, 0, 0, 0]
% Unknown parameters: k1, k2
% Output function: y = nC/(cC + 2*nD)
%

% Model
model = struct;
model.transcription = 'shooting';
model.x = {'VR', 'nA', 'nB', 'nC', 'nD'};
model.p = {'k1', 'k2', 'cBf'};

% Dependent variables with definitions
model.y = {'lc'};
model.h = @(t, v, p) { 1 / (1 + 2*v.nD/max(v.nC, 1e-6))}; % avoid divide-by-zero

% Data and measurements
model.d = {'Qf', 'lc_m'};

% ODE right-hand-side
model.ode = @(t, v, p) {v.Qf,...
                        -p.k1*v.nA*v.nB/v.VR,...
                        v.Qf*p.cBf - v.nB*(p.k1*v.nA + p.k2*v.nC)/v.VR,...
                        v.nB*(p.k1*v.nA - p.k2*v.nC)/v.VR,...
                        p.k2*v.nC*v.nB/v.VR};

% Relative least squares objective function
model.lsq = @(t, v, p) {v.lc_m/v.lc - 1};

% Load data
teaf      = 0.00721;
teaden    = 0.728;
flow = load('flow.dat');
lc = load('lc.dat');
tQf  = [0. ; flow(:,1)];
Qf = [0. ; flow(:,2)./teaden];
tlc    = lc(:,1);
lc = lc(:,2);

% Get all time points occuring in either tlc or tflow
[tout,~,ic] = unique([tQf; tlc]);
Qf_ind = ic(1:numel(tQf));
lc_ind = ic(numel(tQf)+1:end);

% Interpolate lcmeas and Qf to this new grid
Qf = interp1(tQf, Qf, tout, 'previous');
lc_m = interp1(tlc, lc, tout, 'previous');

% Replace NaNs with zeros
Qf(isnan(Qf)) = 0.;
lc_m(isnan(lc_m)) = 0.;

% Initial volume
VR0 = 2370;

% Options
model.tout = tout';
model.lsq_ind = lc_ind'; % only include lc_ind in objective

% Create a paresto instance
pe = paresto(model);

% Solution in the book
nA0 = 2.35;
k1 = 2500;
k2 = 1250;

% Initial guess, upper and lower bounds for the estimated parameters
p0 = [k1; 0.9*k2; teaf];
lbp = [k1; 0.5*k2; teaf];
ubp = [k1; 1.5*k2; teaf];
ic0 = [VR0; 1.1*nA0; 0; 0; 0];
lbic = [VR0; 0.5*nA0; 0; 0; 0];
ubic = [VR0; 1.5*nA0; 0; 0; 0];
lbtheta = [lbp; lbic];
ubtheta = [ubp; ubic];
theta0 = [p0; ic0];

% Estimate parameters
[est,v,p] = pe.optimize([Qf'; lc_m'], theta0, lbtheta, ubtheta);
est.theta
est.d2f_dtheta2

% Also calculate confidence intervals with 95 % confidence
est_ind = find(lbtheta~=ubtheta); % estimate free parameters only
theta_conf = pe.confidence(est, est_ind, 0.95)

disp('Confidence interval')
theta_conf

clf
subplot(2,2,1)
hold on
plot(model.tout, v.nA)
plot(model.tout, v.nC)
plot(model.tout, v.nD)
legend({'n_A', 'n_C', 'n_D'});
xlabel('time (min)')
ylabel('Amount of substance (kmol)')
title('Amount of substance of species A, C and D versus time')

subplot(2,2,2)
hold on
plot(model.tout, v.nB)
legend({'n_B'});
xlabel('time (min)')
ylabel('Amount of substance (kmol)')
title('Amount of substance of species B versus time')

subplot(2,2,3)
stairs(model.tout, v.Qf)
xlabel('time (min)')
ylabel('flowrate (kg/min)')
title('Base addition rate')

subplot(2,2,4)
hold on
plot(model.tout, v.lc)
plot(tlc, lc, 'o')
ylim([0, 2*max(lc)])
legend({'model', 'measurement'});
xlabel('time (min)')
title('LC measurement')
