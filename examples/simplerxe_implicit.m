%%
%% jbr,  3/14/21
%%
%%  simple reaction equilibrium problem
%%  A <-> B  <-> C,   second reaction at equilibrium
%%
%%  model:  dA/dt = -r1, dB/dt = r1 - r2, dC/dt = r2
%%  r1 = k1 A - k_1 B, , r2 = k2 B - k_2 C
%%  measurement: A,B,C

daemodel = struct;
%daemodel.transcription = 'simultaneous';
%daemodel.ord = 1;
daemodel.transcription = 'shooting';
%% daemodel.nlp_solver_options.ipopt.linear_solver = 'ma27';
daemodel.nlp_solver_options.ipopt.mumps_scaling = 0;
% set eps to zero for daeebraic model
daemodel.nlp_solver_options.sens_linsol_options.eps = 0;
daemodel.print_level = 1;
%
daemodel.x = {'A', 'B', 'C'};
daemodel.z = {'Adot', 'Bdot', 'Cdot'};
daemodel.p = {'k1', 'k_1', 'K2'};
daemodel.d = {'Ameas', 'Bmeas', 'Cmeas'};

nplot = 31;
tplot = linspace(0, 3, nplot);
daemodel.tout = tplot;

function retval = alg_part(t, y, p)
  r1 = (p.k1*y.A - p.k_1*y.B);
  retval = {y.Adot + r1, y.Bdot + y.Cdot - r1, p.K2*y.B - y.C};
end%function

daemodel.ode = @(t, y, p) {y.Adot, y.Bdot, y.Cdot};
daemodel.alg = @alg_part;
%daemodel.alg = @(t, y, p) {p.K2*y.Bdot - y.Cdot, p.K2*y.B - y.C};
daemodel.lsq = @(t, y, p) {y.Ameas - y.A, y.Bmeas - y.B, y.Cmeas - y.C};

pe = paresto(daemodel);

%% create measurements from full model

function rhs = fullmodel(t, x, p)
  A = x(1);
  B = x(2);
  C = x(3);
  r = [p.k1*A - p.k_1*B; p.k2*B - p.k_2*C];
  rhs = [-r(1); r(1) - r(2); r(2)];
end%function
%% page 214 in green book
pf = struct();
pf.k1 = 0.5;
pf.k_1 = 0.5;
pf.k2 = 40;
pf.k_2 = 20;
A0 = 1;
B0 = 0.4;
C0 = 0;
xf0 = [A0; B0; C0];

[tplot, xf] = ode15s(@(t, x) fullmodel(t, x, pf), tplot, xf0);
var = 0.001;
randn('seed',0);
meas = xf + sqrt(var)*randn(size(xf));

p = struct();
p.k1 = pf.k1;
p.k_1 = pf.k_1;
p.K2 = pf.k2/pf.k_2;

x0 = xf0;
z0 = [0;0;0];

%% check model
x0 = [A0; (B0+C0)/(1+p.K2); (B0+C0)*p.K2/(1+p.K2)];
[xsim, zsim] = pe.simulate(zeros(3,1), x0, p, z0);

# ysim = [p.A0, p.B0, p.C0] + [xsim', zsim']*stoi;


%% set up parameter initial guesses and bounds
theta0 = struct();
theta0.k1 = 2*p.k1;
theta0.k_1 = 0.5*p.k_1;
theta0.K2 = 2*p.K2;
theta0.A = x0(1);
theta0.B = x0(2);
theta0.C = x0(3);
theta0.Adot = 0;
theta0.Bdot = 0;
theta0.Cdot = 0;

lb = struct();
lb.k1 = 0;
lb.k_1 = 0;
lb.K2 = 0;
lb.A = 0.1*theta0.A;
lb.B = 0.1*theta0.B;
lb.C = 0.1*theta0.C;

ub = struct();
ub.k1 = 10;
ub.k_1 = 10;
ub.K2 = 100;
ub.A = 10*theta0.A;
ub.B = 10*theta0.B;
ub.C = 10*theta0.C;

## %% estimate the parameters

est = pe.optimize(meas', theta0, lb, ub);

figure(1)
plot(tplot, meas', 'o', tplot, est.x)
figure(2)
plot(tplot, est.z)

conf =  pe.confidence(est, 0.95);

disp('Estimated parameters')
disp(est.theta)
disp('Bounding box intervals')
disp(conf.bbox)
