%%
%% jbr,  3/16/21
%%
%%  simple reaction equilibrium problem in CSTR
%%  A <-> B  <-> C,   second reaction at equilibrium
%%
%%  model:  dA/dt = (Af - A)/tau -r1,
%%          dB/dt = (Bf - B)/tau + r1 - r2,
%%          dC/dt = (Cf - C)/tau + r2
%%  r1 = k1 A - k_1 B, , r2 = k2 B - k_2 C
%%  measurement: A,B,C

daemodel = struct;
%daemodel.transcription = 'simultaneous';
%daemodel.ord = 1;
daemodel.transcription = 'shooting';
%% daemodel.nlp_solver_options.ipopt.linear_solver = 'ma27';
daemodel.nlp_solver_options.ipopt.mumps_scaling = 0;
% set eps to zero for daeebraic model
%daemodel.nlp_solver_options.sens_linsol_options.eps = 0;
daemodel.print_level = 1;
%
daemodel.x = {'A', 'Y'};  % Y = B + C 
daemodel.z = {'B'};
daemodel.p = {'k1', 'k_1', 'K2', 'tau', 'Af', 'Bf', 'Cf'};
daemodel.d = {'Ameas', 'Bmeas', 'Cmeas'};

nplot = 31;
tplot = linspace(0, 5, nplot);
daemodel.tout = tplot;

daemodel.ode = @(t, y, p){ (p.Af-y.A)/p.tau - (p.k1*y.A - p.k_1*y.B), ...
			  (p.Bf + p.Cf - y.Y)/p.tau + (p.k1*y.A - p.k_1*y.B) };
daemodel.alg = @(t, y, p) {p.K2*y.B - (y.Y-y.B)};  % K2*C_B - C_C
daemodel.lsq = @(t, y, p) {y.Ameas - y.A, y.Bmeas - y.B, y.Cmeas - (y.Y-y.B)};
pe = paresto(daemodel);

%% create measurements from full model

function rhs = fullmodel(t, x, p)
  A = x(1);
  B = x(2);
  C = x(3);
  r = [p.k1*A - p.k_1*B; p.k2*B - p.k_2*C];
  flow = [p.Af - A; p.Bf - B; p.Cf - C] / p.tau;
  rhs = flow + [-r(1); r(1) - r(2); r(2)];
end%function
%% page 214 in green book
pf = struct();
pf.k1 = 0.5;
pf.k_1 = 0.5;
pf.k2 = 40;
pf.k_2 = 20;
pf.tau = 1;
pf.Af = 1;
pf.Bf = 0.5;
pf.Cf = 0;
Af0 = 1;
Bf0 = 0.4;
Cf0 = 0;
xf0 = [Af0; Bf0; Cf0];

[tplot, xf] = ode15s(@(t, x) fullmodel(t, x, pf), tplot, xf0);
var = 0.001;
randn('seed',0);
meas = xf + sqrt(var)*randn(size(xf));

p = struct();
p.k1 = pf.k1;
p.k_1 = pf.k_1;
p.K2 = pf.k2/pf.k_2;
p.tau = 1;
p.Af = pf.Af;
p.Bf = pf.Bf;
p.Cf = pf.Cf;

x0 = [Af0; Bf0 + Cf0];

%% check model

x0 = [Af0; Bf0+Cf0];
%%z0 = (Bf0+Cf0)/(1+p.K2);
z0 = 0;
[xsim, zsim] = pe.simulate(zeros(3,1), x0, p, z0);

%% set up parameter initial guesses and bounds
theta0 = p;
theta0.k1 = 2*p.k1;
theta0.k_1 = 0.5*p.k_1;
theta0.K2 = 2*p.K2;
theta0.A = Af0;
theta0.Y = [Bf0+Cf0];
%%theta0.B = theta0.Y/(1+theta0.K2);
theta0.B = 0;

lb = p;
lb.k1 = 0;
lb.k_1 = 0;
lb.K2 = 0;
## lb.A = theta0.A;
## lb.Y = theta0.Y;
lb.A = theta0.A;
lb.Y = theta0.Y;

ub = p;
ub.k1 = 10;
ub.k_1 = 10;
ub.K2 = 100;
## ub.A = theta0.A;
## ub.Y = theta0.Y;
ub.A = theta0.A;
ub.Y = theta0.Y;

## %% estimate the parameters

est = pe.optimize(meas', theta0, lb, ub);

yest = [est.x(1,:); est.z; est.x(2,:) - est.z];
figure(1)
plot(tplot, meas', 'o', tplot, yest)

conf =  pe.confidence(est, 0.95);

disp('Estimated parameters')
disp(est.theta)
disp('Bounding box intervals')
disp(conf.bbox)






