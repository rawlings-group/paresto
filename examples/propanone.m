%% Parameter estimation for propanone mechanism

%% convert to use paresto.m
%% jbr, 2/6/2021


%% load casadi and paresto
addpath('~/src/casadi/casadi-octave')
addpath ("~/src/paresto");
pkg ('load',  'statistics');

%% parameter estimation with paresto
model = struct();
model.print_level = 1;
model.transcription = 'shooting';
## model.transcription = 'simultaneous';
## model.ord = 1;

model.nlp_solver_options.ipopt.linear_solver = 'ma27';
%model.nlp_solver_options.ipopt.mumps_scaling = 0;
% set eps to zero for algebraic model
model.nlp_solver_options.sens_linsol_options.eps = 0;
%% variable list: differential and algebraic states, parameters, and
%% measurements and constants
model.x = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};
model.z = {'X', 'Y', 'Z'};
model.p = {'k1', 'K_1', 'K4', 'K5', 'K6', 'K8'};
model.d = {'mA', 'mB', 'mC', 'mD', 'mE', 'mF', 'mG', 'mH'};

% Create suitable data such that:
% Major: A,B,C,D,H
% Minor: E,F,G
% QSSA: X,Y,Z,
% Scaling: Xbar = k2 X, Ybar = k3 Y, Zbar = k6 Z
maj = [1:4, 8];
mino = [5:7];
int = [9:11];
%% Kinetic parameters for full model
pf.k1 = 10^-2;
pf.k_1 = 10^7;
pf.k2 = 10^3;
pf.k3 = 10^2;
pf.k4 = 0.5*10^5;
pf.k5 = 10^3;
pf.k6 = 10^1;
pf.k7 = 10^-1;
pf.k8 = 10^9;
%% scaled parameters for reduced DAE model
p.k1 = pf.k1;
p.K_1 = pf.k_1/(pf.k2*pf.k3);
p.K4 = pf.k4/(pf.k3^2);
p.K5 = pf.k5/(pf.k3*pf.k7);
p.K6 = pf.k6/(pf.k7^2);
p.K8 = pf.k8/(pf.k2^2);

% Initial Condition:
nx = numel(model.x);
nz = numel(model.z);
x0 = zeros(nx + nz, 1);
x0(1) = 1;
%%
function xdot = massbal(t, x, p)
  A = x(1);
  B = x(2);
  C = x(3);
  D = x(4);
  E = x(5);
  F = x(6);
  G = x(7);
  H = x(8);
  X = x(9);
  Y = x(10);
  Z = x(11);
  %%
  r= [p.k1*A - p.k_1*X*Y;
    p.k2*X;
    p.k3*Y*A;
    p.k4*Y^2;
    p.k5*Y*Z;
    p.k6*Z^2;
    p.k7*Z;
    p.k8*X^2];
  %%
  R = [-r(1)-r(3), r(2), r(4), r(3), r(5), r(6), r(7), r(8), ...
       r(1)-r(2)-2*r(8), r(1)+r(2)-r(3)-2*r(4)-r(5)+r(7), r(3)-r(5)-2*r(6)-r(7)]';
  xdot = R;
end%function

% Solve full model to make data
tfinal = 200;
nplot = 101;
tplot = linspace(0, tfinal, nplot)';
[tplot, xf] = ode15s(@(t, x) massbal(t, x, pf), tplot, x0);

%%%% DAE Model %%%%
function res = daemodel(t, x, xdot, p)
  A = x(1);
  B = x(2);
  C = x(3);
  D = x(4);
  E = x(5);
  F = x(6);
  G = x(7);
  H = x(8);
  X = x(9);
  Y = x(10);
  Z = x(11);
  %%
  r = [p.k1*A - p.K_1*X*Y;
       X;
       Y*A;
       p.K4*Y^2;
       p.K5*Y*Z;
       p.K6*Z^2;
       Z;
       p.K8*X^2];
  %%
  R = [-r(1)-r(3), r(2), r(4), r(3), r(5), r(6), r(7), r(8), ...
       r(1)-r(2)-2*r(8), r(1)+r(2)-r(3)-2*r(4)-r(5)+r(7), r(3)-r(5)-2*r(6)-r(7)]';
  res = xdot - R;
  %% Fix up algebraic states
  int = 9:11;
  res(int) = R(int);
endfunction

%% Use decic to get consistent ICs
%% fix ODE states; free algebraic states
x0g = x0;
%% make good guess on scaled intermediate (X,Y,Z) so decic finds positive IC
x0g(int) = xf(2,int).*[pf.k2, pf.k3, pf.k7];
xp0g = massbal(0, x0, pf);
fixed_x0 = [ones(8,1); zeros(3,1)];
fixed_xp0 = zeros(size(xp0g));
options.AbsTol = 0.01*sqrt(eps);
options.RelTol = 0.01*sqrt(eps);
[x0new, xp0new, resnorm] = decic (@(t, x, xp) daemodel(t, x, xp, p), ...
				  0, x0g, fixed_x0, xp0g, fixed_xp0, options);
%% solve the DAEs
[tplot, x] = ode15i (@(t, x, xp) daemodel(t, x, xp, p), tplot, x0new, xp0new);
%% unscale the QSSA species
xunscaled = x;
xunscaled(:,int) = xunscaled(:,int)./[pf.k2; pf.k3; pf.k7]';
%% Plot full model (no noise) and dae model solutions (no noise)
figure(1);
subplot(1,3,1);
plot(tplot, xf(:,maj), '-', tplot, xunscaled(:,maj), '.');
xlabel('time');
ylabel('concentration');
title('Major Species');
legend('A','B','C','D','H');
%%
subplot(1,3,2);
plot(tplot,xf(:,mino), '-', tplot, xunscaled(:, mino), '.');
xlabel('time');
ylabel('concentration');
title('Minor Species');
legend('E','F','G');
%%
subplot(1,3,3);
%semilogy(tplot,xf(:, int), '-', tplot, xunscaled(:,int), '.');
plot(tplot,xf(:, int), '-', tplot, xunscaled(:,int), '.');
xlabel('time');
ylabel('concentration');
title('Reactive intermediates');
legend('X','Y','Z');
%% Add noise
meas = sort([maj, mino]);
randn('seed',0);
y = xf(:, meas);
stddev = sqrt(0.00001);
ymeas = y+stddev*randn(size(y));
%% Plot full model (without noise) and full model (with noise)
figure(2);
subplot(1,3,1);
plot(tplot,xf(:, maj), '-', tplot, ymeas(:, maj), '-.');
xlabel('time');
ylabel('concentration');
title('Major Species');
legend('A','B','C','D','H');
%%
subplot(1,3,2);
plot(tplot,xf(:, mino), '-', tplot, ymeas(:, mino), '-.');
xlabel('time');
ylabel('concentration');
title('Minor Species');
legend('E','F','G');
%%
subplot(1,3,3);
plot(tplot,xf(:, int), '-', tplot, 0*xf(:, int), '-.');
xlabel('time');
ylabel('concentration');
title('Reactive intermediates');
legend('X','Y','Z');

myfile = fopen('propanone_data.dat', 'w');
for i = 1:numel(tplot)
  fprintf(myfile, '%8.3f', tplot(i), ymeas(i,:));
  fprintf(myfile, '\n');
end
fclose(myfile);

function r = rxrates(t, x, p)
  r = [p.k1*x.A - p.K_1*x.X*x.Y;
       x.X;
       x.Y*x.A;
       p.K4*x.Y^2;
       p.K5*x.Y*x.Z;
       p.K6*x.Z^2;
       x.Z;
       p.K8*x.X^2];
end%function

function xdot = ode_part(t, x, p)
  r = rxrates(t, x, p);
  R = {-r(1)-r(3), r(2), r(4), r(3), r(5), r(6), r(7), r(8)};, ...
  xdot =  R;
end%function
function res = alg_part(t, x, p)
  r = rxrates(t, x, p);
  R = {r(1)-r(2)-2*r(8), r(1)+r(2)-r(3)-2*r(4)-r(5)+r(7), r(3)-r(5)-2*r(6)-r(7)};
  res =  R;
end%function
model.ode = @ode_part;
model.alg = @alg_part;
model.lsq = @(t, x, p) {x.mA-x.A, x.mB-x.B, x.mC-x.C, x.mD-x.D, x.mE-x.E, ...
  			x.mF-x.F, x.mG-x.G, x.mH-x.H, 1e-10*x.X, 1e-10*x.Y, 1e-10*x.Z};
## model.lsq = @(t, x, p) {x.mA-x.A, x.mB-x.B, x.mC-x.C, x.mD-x.D, x.mE-x.E, ...
##   			x.mF-x.F, x.mG-x.G, x.mH-x.H};
model.tout = tplot;
pe = paresto(model);

small = 0.1;
big = 10;
%% parameters; initial guesses and bounds;
for i = 1:numel(model.p)
  fni = model.p{i};
  pp = p.(fni);
  thetaic.(fni) = pp;
  lb.(fni) = small*pp;
  ub.(fni) = big*pp;
end%for

for i = 1: numel(model.x)
  fnx = model.x{i};
  xx = x0(i);
  thetaic.(fnx) = xx;
  lb.(fnx) = xx;
  ub.(fnx) = xx;
end%for

for i = 1: numel(model.z)
  fnz = model.z{i};
  thetaic.(fnz) = 1e-3;
  lb.(fnz) = 0;
  ub.(fnz) = 10;
end%for

%% check reduced model
ysim = pe.simulate(zeros(8,1), x0new(1:8), p, x0new(9:11));

%% estimate the parameters
[est, y, p] = pe.optimize(ymeas', thetaic, lb, ub);

% Also calculate confidence intervals with 95% confidence

conf = pe.confidence(est, 0.95);

disp('Estimated parameters')
disp(est.theta)
disp('Bounding box intervals')
disp(conf.bbox)

