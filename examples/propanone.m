%% Parameter estimation for propanone mechanism

%% convert to use paresto.m
%% jbr, 2/6/2021


%% load casadi and paresto
addpath('~/src/casadi/casadi-octave')
addpath ("~/src/paresto");
pkg ('load',  'statistics');

% Create suitable data such that:
% Major: A,B,C,D,H
% Minor: E,F,G
% QSSA: X,Y,Z,
% Scaling: Xbar = k2 X, Ybar = k3 Y, Zbar = k6 Z
maj = [1:4, 8];
mino = [5:7];
int = [9:11];
p.int = int;
%% Kinetic parameters
p.k1 = 10^-2;
p.k_1 = 10^7;
p.k2 = 10^3;
p.k3 = 10^2;
p.k4 = 0.5*10^5;
p.k5 = 10^3;
p.k6 = 10^1;
p.k7 = 10^-1;
p.k8 = 10^9;
%% scaled parameters for DAE model
p.K_1 = p.k_1/(p.k2*p.k3);
p.K4 = p.k4/(p.k3^2);
p.K5 = p.k5/(p.k3*p.k7);
p.K6 = p.k6/(p.k7^2);
p.K8 = p.k8/(p.k2^2);
% Initial Condition:
B0 = 0;
C0 = 0;
A0 = 1;
D0 = 0;
H0 = 0;
E0 = 0;
F0 = 0;
G0 = 0;
X0 = 0;
Y0 = 0;
Z0 = 0;
%%
x0 = [A0; B0; C0; D0; E0; F0; G0; H0; X0; Y0; Z0];
%%
p.nu = [-1 +0 +0 +0 +0 +0 +0 +0 +1 +1 +0;
	+0 +1 +0 +0 +0 +0 +0 +0 -1 +1 +0;
	-1 +0 +0 +1 +0 +0 +0 +0 +0 -1 +1;
	+0 +0 +1 +0 +0 +0 +0 +0 +0 -2 +0;
	+0 +0 +0 +0 +1 +0 +0 +0 +0 -1 -1;
	+0 +0 +0 +0 +0 +1 +0 +0 +0 +0 -2;
	+0 +0 +0 +0 +0 +0 +1 +0 +0 +1 -1;
	+0 +0 +0 +0 +0 +0 +0 +1 -2 +0 +0];
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
  xdot = p.nu'*r;
end%function

% Solve full model to make data
tfinal = 200;
nplot = 101;
tplot = linspace(0, tfinal, nplot)';
[tplot, xf] = ode15s(@(t, x) massbal(t, x, p), tplot, x0);

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
  R = p.nu'*r;
  res = xdot - R;
  %% Fix up algebraic states
  res(p.int,:) = R(p.int,:);
endfunction

%% Use decic to get consistent ICs
%% fix ODE states; free algebraic states
x0g = x0;
%% make good guess on scaled intermediate (X,Y,Z) so decic finds positive IC
x0g(int) = xf(2,int).*[p.k2, p.k3, p.k7];
xp0g = massbal(0, x0, p);
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
xunscaled(:,int) = xunscaled(:,int)./[p.k2;p.k3;p.k7]';
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
for i = 1: size(tplot)
  fprintf(myfile, '%8.3f', tplot(i), ymeas(i,:));
  fprintf(myfile, '\n');
end
fclose(myfile);


%% parameter estimation with paresto
model = struct();
model.print_level = 1;
model.transcription = 'shooting';
## model.transcription = 'simultaneous';
## model.ord = 1;

model.nlp_solver_options.ipopt.linear_solver = 'ma27';
%model.nlp_solver_options.ipopt.mumps_scaling = 0;
% set eps to zero for algebraic model
algmodel.nlp_solver_options.sens_linsol_options.eps = 0;
%% variable list: differential and algebraic states, parameters, and
%% measurements and constants
model.x = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};
model.z = {'X', 'Y', 'Z'};
model.p = {'k1', 'K_1', 'K4', 'K5', 'K6', 'K8'};
model.d = {'mA', 'mB', 'mC', 'mD', 'mE', 'mF', 'mG', 'mH'};
function xdot = ode_part(t, x, p)
  r = [p.k1*x.A - p.K_1*x.X*x.Y;
       x.X;
       x.Y*x.A;
       p.K4*x.Y^2;
       p.K5*x.Y*x.Z;
       p.K6*x.Z^2;
       x.Z;
       p.K8*x.X^2];
  R = {-r(1)-r(3), r(2), r(4), r(3), r(5), r(6), r(7), r(8)};, ...
  xdot =  R;
endfunction
function res = alg_part(t, x, p)
  r = [p.k1*x.A - p.K_1*x.X*x.Y;
       x.X;
       x.Y*x.A;
       p.K4*x.Y^2;
       p.K5*x.Y*x.Z;
       p.K6*x.Z^2;
       x.Z;
       p.K8*x.X^2];
  R = {r(1)-r(2)-2*r(8), r(1)+r(2)-r(3)-2*r(4)-r(5)+r(7), r(3)-r(5)-2*r(6)-r(7)};
  res =  R;
endfunction
model.ode = @ode_part;
model.alg = @alg_part;
model.lsq = @(t, x, p) {x.mA-x.A, x.mB-x.B, x.mC-x.C, x.mD-x.D, x.mE-x.E, ...
  			x.mF-x.F, x.mG-x.G, x.mH-x.H, 1e-10*x.X, 1e-10*x.Y, 1e-10*x.Z};
## model.lsq = @(t, x, p) {x.mA-x.A, x.mB-x.B, x.mC-x.C, x.mD-x.D, x.mE-x.E, ...
##   			x.mF-x.F, x.mG-x.G, x.mH-x.H};
model.tout = tplot;
pe = paresto(model);

%% parameters; initial guesses and bounds;
thetaic.k1 = p.k1;
thetaic.K_1 = p.K_1;
thetaic.K4 = p.K4;
thetaic.K5 = p.K5;
thetaic.K6 = p.K6;
thetaic.K8 = p.K8;

## thetaic.A = x0new(1);
## thetaic.B = x0new(2);
## thetaic.C = x0new(3);
## thetaic.D = x0new(4);
## thetaic.E = x0new(5);
## thetaic.F = x0new(6);
## thetaic.G = x0new(7);
## thetaic.H = x0new(8);
## thetaic.X = x0new(9);
## thetaic.Y = x0new(10);
## thetaic.Z = x0new(11);

thetaic.A = x0(1);
thetaic.B = x0(2);
thetaic.C = x0(3);
thetaic.D = x0(4);
thetaic.E = x0(5);
thetaic.F = x0(6);
thetaic.G = x0(7);
thetaic.H = x0(8);
thetaic.X = 1e-3;
thetaic.Y = 1e-3;
thetaic.Z = 1e-3;

lb = thetaic;
ub = thetaic;

small = 0.1;
%small = 1e-3;
lb.k1  = small*thetaic.k1;
lb.K_1 = small*thetaic.K_1;
lb.K4  = small*thetaic.K4;
lb.K5  = small*thetaic.K5;
lb.K6  = small*thetaic.K6;
lb.K8  = small*thetaic.K8;
big = 10;
%big = 1e3;
ub.k1  = big*thetaic.k1;
ub.K_1 = big*thetaic.K_1;
ub.K4  = big*thetaic.K4;
ub.K5  = big*thetaic.K5;
ub.K6  = big*thetaic.K6;
ub.K8  = big*thetaic.K8;

lb.X = 0;
lb.Y = 0;
lb.Z = 0;
ub.X = 10;
ub.Y = 10;
ub.Z = 10;


%% check model
%% need a vector of parameters for simulate function
pp = [p.k1; p.K_1; p.K4; p.K5; p.K6; p.K8];
ysim = pe.simulate(zeros(8,1), x0new(1:8), pp, x0new(9:11));

%% estimate the parameters
[est, y, p] = pe.optimize(ymeas', thetaic, lb, ub);

% Also calculate confidence intervals with 95% confidence

conf = pe.confidence(est, 0.95);

disp('Estimated parameters')
disp(est.theta)
disp('Bounding box intervals')
disp(conf.bbox)

