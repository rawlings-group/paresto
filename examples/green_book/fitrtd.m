% [depends] ch9selecrtd.dat

%% Copyright (C) 2001, James B. Rawlings and John G. Ekerdt
%%
%% This program is free software; you can redistribute it and/or
%% modify it under the terms of the GNU General Public License as
%% published by the Free Software Foundation; either version 2, or (at
%% your option) any later version.
%%
%% This program is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; see the file COPYING.  If not, write to
%% the Free Software Foundation, 59 Temple Place - Suite 330, Boston,
%% MA 02111-1307, USA.

%% Revised 7/12/02
%% Created 7/12/02
%%
%% actual rho = 0.62

%% Converted to use parest.m file; 
%% 3/31/2007

%% converted to use paresto.m;
%% jbr, 4/23/2018
%% jbr, new paresto struct for parameters, 5/24/2019
%%

%tmp = load ('../../ch8/figures/selecrtd.dat');
tmp = load ('ch9selecrtd.dat');
if (isstruct (tmp))
  table2 = tmp.table2;
else
  table2 = tmp;
endif

model = struct;
model.nlp_solver_options.ipopt.linear_solver = 'ma27';
model.x = {'c1', 'c2'};
model.p = {'rho', 'alpha', 'cf', 'tau1', 'tau2'};
model.d = {'m_c2'};

model.ode = @(t, x, p) {(p.alpha*p.cf-(p.alpha+p.rho)*x.c1+p.rho*x.c2)/p.tau1, ...
		      ((p.alpha+p.rho)*x.c1-(1+p.alpha+p.rho)*x.c2)/p.tau2 };
model.lsq = @(t,y,p) {y.c2 - y.m_c2};

c10 = 0;
c20 = 0;
rho0 = 0.5;
alpha = 0.1;
cf = 1;
tau1 = 1;
tau2 = 2;
parac = [rho0; alpha; cf; tau1; tau2];

tplot = linspace(0,32,75)';
tmeas = table2(:,1);
ymeas = table2(:,2)';

[tout,~,ic] = unique([tmeas; tplot]);
% index of times at which measurement is made
meas_ind = ic(1:numel(tmeas));
model.tout = tout;
% interploate measurement onto new grid
y_noisy = interp1(tmeas, ymeas, tout, 'previous');
y_noisy(isnan(y_noisy)) = 0.;
model.lsq_ind = meas_ind'; % only use index of measurement times in objective function

pe = paresto(model);

est_ind = 1;

par = struct();
par.rho = 0.5
par.alpha = 0.1;
par.cf = 1;
par.tau1 = 1;
par.tau2 = 2;
par.c1 = 0;
par.c2 = 0;

lb = struct();
lb.rho = 0.1;
lb.alpha = par.alpha;
lb.cf = par.cf;
lb.tau1 = par.tau1;
lb.tau2 = par.tau2;

ub = lb;
ub.rho = 0.95;

% Estimate parameters
[est, y, p] = pe.optimize(y_noisy', par, lb, ub);

% Also calculate confidence intervals with 95 % confidence
theta_conf = pe.confidence(est, numel(est_ind), 0.95);


disp('Estimated parameters and confidence intervals')
[est.theta(est_ind), theta_conf]

c2dim = y.c2/cf;
gnuplotsave('fitrtd.dat', struct('data', [tmeas(:), ymeas(:)], ...
                                 'fit', [model.tout(:), c2dim(:)]));

if (~ strcmp (getenv ('OMIT_PLOTS'), 'true')) %% PLOTTING
plot (tmeas,  ymeas, 'o', model.tout, c2dim');
%% TITLE
endif %% PLOTTING

%% Now solve for the steady-state conversion and yield with the estimated rho
par = struct;
par.rho = est.theta(est_ind);
par.k1 = 1;
par.k2 = 2;
par.caf = 1;
par.cbf = 1;
par.n   = 2;
par.tau1 = tau1;
par.tau2 = tau2;
par.alpha = alpha;
function res = tworeac (x, p)
  ca1 = x(1);
  cb1 = x(2);
  ca2 = x(3);
  cb2 = x(4);
  r11 = p.k1*ca1*cb1;
  r21 = p.k2*ca1^p.n;
  r12 = p.k1*ca2*cb2;
  r22 = p.k2*ca2^p.n;
  R = [ -(r11+r21)*p.tau1;
        -r11*p.tau1;
        -(r12+r22)*p.tau2;
        -r12*p.tau2];
  flow = [ p.alpha*p.caf-(p.alpha+p.rho)*ca1+p.rho*ca2;
           -(p.alpha+p.rho)*cb1+p.rho*cb2;
           -(1+p.alpha)*ca2+(p.alpha+p.rho)*ca1-p.rho*ca2;
           p.cbf+(p.alpha+p.rho)*cb1-p.rho*cb2-(1+p.alpha)*cb2 ];
  res = flow + R;
endfunction
x0 = [0; par.cbf-par.caf; 0; par.cbf-par.caf];
[x, fval, info] = fsolve(@(x) tworeac(x,par), x0);
conv  = (par.alpha*par.caf-(1+par.alpha)*x(3))/(par.alpha*par.caf)
yield = (par.cbf-(1+par.alpha)*x(4))/(par.alpha*par.caf-(1+par.alpha)*x(3))
