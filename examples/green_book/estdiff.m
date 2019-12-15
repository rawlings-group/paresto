%% [depends] batch_noise.dat
%%
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

%%
%% jbr,  11/6/03, taken largely from nthorder.m
%%
%% converted to paresto.m
%% jbr, 4/27/2018
%%

model = struct;
model.transcription = 'shooting';
model.x = {'ca'};
model.p = {'k', 'n'};
model.d = {'m_ca'};

model.ode = @(t, y, p) {-p.k * y.ca^p.n};

model.lsq = @(t, y, p) {y.ca - y.m_ca};

data = load ('batch_noise.dat');
tdata = data(:,1);
ntimes = length(tdata);
y_noisy = data(:,2:4);
model.tout = tdata;

kac   = 0.05;
ca0ac = 0.83;
nac   = 1.45;

pe = paresto(model);

thetaac = [kac; nac; ca0ac];
x_ac = ca0ac;
p_ac = [kac; nac];

small = 1e-3; 
large  = 3;
est_ind = 1:3;

% Initial guess, upper and lower bounds for the estimated parameters
theta0 = struct();   
theta0.k = kac;      
theta0.n = nac;      
theta0.ca = ca0ac;

ubtheta = struct(); 
ubtheta.k = large;	 
ubtheta.n = large;
ubtheta.ca = large;

lbtheta = struct(); 
lbtheta.k = small;  
lbtheta.n = small;  
lbtheta.ca = small;

nsets = size(y_noisy,2);
thetaest = NaN(numel(thetaac), nsets);
caest = NaN(size(y_noisy));
bbox = thetaest;

for j = 1:nsets
  % Estimate parameters
  [est, y, p] = pe.optimize(y_noisy(:,j)', theta0, lbtheta, ubtheta);
  % Also calculate confidence intervals with 95 % confidence
  theta_conf = pe.confidence(est, est_ind, 0.95);
  disp('Estimated parameters and confidence intervals')
  [est.theta(est_ind), theta_conf]
  thetaest(:,j) = est.theta(est_ind);
  caest(:,j) = y.ca;
  bbox(:,j) = theta_conf;
endfor

table  = [model.tout, caest];
save estdiff.dat table data;
if (~ strcmp (getenv ('OMIT_PLOTS'), 'true')) %% PLOTTING
subplot(3,1,1)
plot (table(:,1), table(:,2), data(:,1), data(:,2), '+')
%% TITLE estdiff_1

subplot(3,1,2)
plot (table(:,1), table(:,3), data(:,1), data(:,3), '+')
%% TITLE estdiff_2

subplot(3,1,3)
plot (table(:,1), table(:,4), data(:,1), data(:,4), '+')
%% TITLE estdiff_3
endif %% PLOTTING
