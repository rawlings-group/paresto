% [depends] batch_data.dat

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

%% Revised 7/22/02
%% Created 7/12/02

%% converted to use paresto, jbr, 4/27/2018
%%


tmp = load ('batch_data.dat');
if (isstruct (tmp))
  table = tmp.table;
else
  table = tmp;
endif
time = table(:,1);
ymeas = table(:,2);
pts = length(ymeas);


%% PART B
%% First estimate the parameters from the linearized transformation via
%% least squares
rate = (ymeas(2:pts)-ymeas(1:pts-1))./(time(2:pts)-time(1:pts-1));
bigA = [log(ymeas(1:pts-1)) ones(pts-1,1)];
x = inv(bigA'*bigA)*bigA'*log(-rate);
n_linear = x(1);
k_linear = exp(x(2));
theta_linear = [n_linear; k_linear];

function xdot = rxrate(t, x, p)
  xdot = -p(2)*x^p(1);
end%function  

[tout, ca_linear] = ode15s(@(t, x) rxrate(t, x, theta_linear), time, ymeas(1));

%% PART C
%% nonlinear parameter estimation

model = struct;
model.transcription = 'shooting';
model.x = {'ca'};
model.p = {'n','k'};
model.d = {'m_ca'};
model.tout = time;

model.ode = @(t, y, p) {-p.k*y.ca^p.n};

model.lsq = @(t, y, p) {y.ca - y.m_ca};

est_ind = 1:2;
small = 1e-3; 
large  = 3;

pe = paresto(model);

objective.parlb = [small; small];
objective.parub = [large; large];
% true parameters are [n; k; ca0] = [1.45; .05; 0.83]; see batch_data.m file

%theta0 = [3; .1; 0.8];
theta0 = struct;
theta0.n = 3;
theta0.k = 0.1;
theta0.ca = 0.8;

lbtheta = theta0;
lbtheta.n = small;
lbtheta.k = small;
lbtheta.ca = 0.79;

ubtheta = theta0;
ubtheta.n = large;
ubtheta.k = large;
ubtheta.ca = 0.81;

% Estimate parameters
[est, y, p] = pe.optimize(ymeas', theta0, lbtheta, ubtheta);

% Also calculate confidence intervals with 95 % confidence
theta_conf = pe.confidence(est, est_ind, 0.95);
disp('Estimated parameters and confidence intervals')
[est.theta(est_ind), theta_conf]


ca = y.ca;

%%echo results
n_linear
k_linear
n = est.theta(1)
k = est.theta(2)

data = [time ymeas ca_linear(:) ca(:)];
save -ascii batch_data_solution.dat data;

if (~ strcmp (getenv ('OMIT_PLOTS'), 'true')) %% PLOTTING
  plot (time, ymeas, '+', time, [ca_linear(:), ca(:)]);
%% TITLE
endif %% PLOTTING
