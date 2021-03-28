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
daemodel.nlp_solver_options.ipopt.linear_solver = 'ma27';
%daemodel.nlp_solver_options.ipopt.mumps_scaling = 0;
% set eps to zero for daeebraic model
daemodel.nlp_solver_options.sens_linsol_options.eps = 0;
daemodel.print_level = 1;
%
daemodel.x = {'ex1'};
daemodel.z = {'ex2'};
%daemodel.p = {'k', 'b'};
daemodel.p = {'k1', 'k_1', 'K2', 'A0', 'B0', 'C0'};
daemodel.d = {'Ameas', 'Bmeas', 'Cmeas'};

nplot = 31;
tplot = linspace(0, 3, nplot);
daemodel.tout = tplot;

daemodel.ode = @(t, y, p) {p.k1*(p.A0-y.ex1)-p.k_1*(p.B0+y.ex1-y.ex2)};
daemodel.alg = @(t, y, p) {p.K2*(p.B0+y.ex1-y.ex2) - (p.C0+y.ex2)};
daemodel.lsq = @(t, y, p) {y.Ameas - (p.A0-y.ex1), y.Bmeas - (p.B0+y.ex1-y.ex2),...
			   y.Cmeas - (p.C0+y.ex2), 1e-10*y.ex1, 1e-10*y.ex2};
pe = paresto(daemodel);

%% create measurements from full model

function rhs = fullmodel(t, x, p)
  A = p.A0 - x(1);
  B = p.B0 + x(1) - x(2);
  C = p.C0 + x(2);
  rhs = [p.k1*A - p.k_1*B; p.k2*B - p.k_2*C];
end%function
%% page 214 in green book
pf.k1 = 0.5;
pf.k_1 = 0.5;
pf.k2 = 40;
pf.k_2 = 40;
pf.A0 = 1;
pf.B0 = 0.4;
pf.C0 = 0;
xf0 = zeros(2,1);

[tplot, xf] = ode15s(@(t, x) fullmodel(t, x, pf), tplot, xf0);
stoi = [-1, 1, 0; 0, -1, 1];
yf = [pf.A0, pf.B0, pf.C0] + xf*stoi;
var = 0.001;
meas = yf + sqrt(var)*randn(size(yf));

p.k1 = pf.k1;
p.k_1 = pf.k_1;
p.K2 = pf.k2/pf.k_2;
p.A0 = pf.A0;
p.B0 = pf.B0;
p.C0 = pf.C0;

%%consistent initial conditions
x0 = 0;
%% set r_2(x0, z0) = 0 to find z0; K2(B0 + x0 - z0) - (C0 + z0)=0
z0 = (p.K2*(pf.B0 + x0) - pf.C0 ) / (1 + p.K2);


%% solve reduced model with dae solver
function res = reducedmodel(t, x, xdot, p)
  A = p.A0 - x(1);
  B = p.B0 + x(1) - x(2);
  C = p.C0 + x(2);
  r = [p.k1*A - p.k_1*B; p.K2*B - C];
  res(1) = xdot(1) - r(1);
  res(2) = r(2);
end%function

xr0 = [x0; z0];
xp0g = fullmodel(0, xr0, pf);
## fixed_x0 = [1;0];
## fixed_xp0 = [0;0];
## options.AbsTol = 0.01*sqrt(eps);
## options.RelTol = 0.01*sqrt(eps);
## [x0new, xp0new, resnorm] = decic (@(t, x, xp) reducedmodel(t, x, xp, p), ...
## 				  0, x0g, fixed_x0, xp0g, fixed_xp0, options);
[tplot, xr] = ode15i(@(t, x, xdot) reducedmodel(t, x, xdot, p), tplot, xr0, xp0g);
  

%% check model
[xsim, zsim] = pe.simulate(zeros(3,1), x0, p, z0);

ysim = [p.A0, p.B0, p.C0] + [xsim', zsim']*stoi;


%% set up parameter initial guesses and bounds
theta0 = p;
theta0.ex1 = 0;
theta0.ex2 = 0;

lb = struct();
lb.k1 = 0;
lb.k_1 = 0;
lb.K2 = 0;
lb.A0 = 0;
lb.B0 = 0;
lb.C0 = 0;
lb.ex1 = 0;
lb.ex2 = 0;

ub = struct();
ub.k1 = 10;
ub.k_1 = 10;
ub.K2 = 100;
ub.A0 = 10;
ub.B0 = 10;
ub.C0 = 10;
ub.ex1 = 0;
ub.ex2 = 0;


## %% estimate the parameters

est = pe.optimize(meas', theta0, lb, ub);

yest = [est.theta.A0, est.theta.B0, est.theta.C0] + [est.x', est.z']*stoi;
plot(tplot, meas', 'o', tplot, yest)

conf =  pe.confidence(est, 0.95);

disp('Estimated parameters')
disp(est.theta)
disp('Bounding box intervals')
disp(conf.bbox)






