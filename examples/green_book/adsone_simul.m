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
%% jbr,  11/28/01
%%
%  modified to use paresto, jbr, 6/8/18
%

Hgas1 = [2.659, 27.590, 141.038, ...
	332.613, 546.857, 608.583, 875.353, 991.072, 1080.134, ...
	1142.782, 1203.178, 1256.244]';

Hads1 = [32.31, 40.02, 45.40, 48.49 ...
	50.19, 53.86, 52.77, 53.19, 53.67, 54.15, 54.51, 54.73]';


Hgas2 = [3.755, 36.331, 189.904, ...
	 413.881, 624.665, 817.592]';

Hads2 = [33.17, 40.90, 45.59, 47.72, ...
	 48.74, 49.43]';

Hgas3 = [12.794, 119.331, 339.917, ...
	 545.227, 735.746, 890.642]';

Hads3 = [36.91, 42.90, 45.65, 47.00, ...
	 47.82, 48.60]';

Hgas4 = [7.424, 73.312, 251.177, ...
	 491.691, 707.806, 873.926]';

Hads4 = [39.11, 46.06, 49.87, 51.57, ...
	 52.48, 53.18]';

Hgas5 = [9.176, 72.065, 259.002, ...
	 475.250, 670.034, 818.669]';

Hads5 = [39.52, 46.06, 49.92, 51.86, ...
	 52.96, 53.88]';

Hgas6 = [2.962, 34.350, 176.860, 405.233, ...
	 607.476, 799.723, 946.765]';

Hads6 = [32.70, 40.11, 44.92, 47.40, 48.62, ...
	 49.50, 50.12]';

Hgas7 = [3.159, 38.528, 200.877, 433.924, ...
	 647.630, 820.673]';

Hads7 = [32.42, 39.81, 44.51, 46.79, 48.14, ...
	 48.95]';

Hgas8 = [34.127, 172.316, ...
	 385.633, 616.149, 812.932, 961.191, 1061.967, 1155.829, 1219.299]';

Hads8 = [43.79, 48.97, 51.57, ...
	 53.09, 54.05, 54.67, 55.20, 55.60, 55.93]';


%% choose a dataset or datasets to fit
cadsmeas = [Hads4];
cgmeas = [Hgas4];

Ks0 = 2;
cms0 = 40;
cmp0 = 0;

algmodel = struct;
%% daemodel.transcription = 'simultaneous';
%% daemodel.ord = 1;
algmodel.transcription = 'simultaneous';
%% algmodel.nlp_solver_options.ipopt.linear_solver = 'ma27';
algmodel.nlp_solver_options.ipopt.mumps_scaling = 0;
%% set eps to zero for algebraic model
algmodel.nlp_solver_options.sens_linsol_options.eps = 0;
algmodel.print_level = 1;
algmodel.z = {'cg', 'cads'};
algmodel.p = {'sqKs', 'cms', 'cmp'};
algmodel.d = {'cgmeas', 'cadsmeas'};

%% least squares
algmodel.alg = @(t, y, p) {y.cads - p.cmp - ...
			   (p.cms*p.sqKs*sqrt(y.cg))/(1 + p.sqKs*sqrt(y.cg)), ...
			   y.cg - y.cgmeas};
algmodel.lsq = @(t, y, p) {y.cadsmeas-y.cads};
%% error in variables, both cg and cads have error
%% algmodel.alg = @(t, y, p) {y.cads - p.cmp - ...
%%    			   (p.cms*p.sqKs*sqrt(y.cg))/(1 + p.sqKs*sqrt(y.cg))};
%% algmodel.lsq = @(t, y, p) {(y.cgmeas-y.cg)*10, (y.cadsmeas-y.cads)};

algmodel.tout = 1:numel(cgmeas);

theta0 = struct;
theta0.sqKs = sqrt(Ks0);
theta0.cms = cms0;
theta0.cmp = cmp0;
theta0.cg = 100;
theta0.cads = 30;

lb = struct;
lb.sqKs = sqrt(1E-3);
lb.cms = 10;
lb.cmp = 0;

ub = struct;
ub.sqKs = sqrt(5);
ub.cms = 200;
ub.cmp = 100;

pe = paresto(algmodel);

est = pe.optimize([cgmeas'; cadsmeas'], theta0, lb, ub);

conf =  pe.confidence(est, 0.95);

disp('Estimated parameters')
disp(est.theta)
disp('Bounding box intervals')
disp(conf.bbox)


nplot = 100;
Xplot = linspace(0,1.05*max(cgmeas),nplot)';
sqKs = est.theta.sqKs;
cms = est.theta.cms;
cmp = est.theta.cmp;
Yplot = cms*sqKs*sqrt(Xplot)./(1+sqKs*sqrt(Xplot)) + cmp;

table1 = [Xplot, Yplot];
table2 = [cgmeas, cadsmeas];
save adsone.dat table1 table2;

if (~ strcmp (getenv ('OMIT_PLOTS'), 'true')) %% PLOTTING
plot (table1(:,1), table1(:,2), table2(:,1), table2(:,2), '+');
%% TITLE
endif %% PLOTTING
