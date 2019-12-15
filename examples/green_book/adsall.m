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
cgmeas = [Hgas1; Hgas2; Hgas3; Hgas4; Hgas5; Hgas6; Hgas7; Hgas8];
cadsmeas = [Hads1; Hads2; Hads3; Hads4; Hads5; Hads6; Hads7; Hads8];

Ks0 = 2;
cms0 = 40;
cmp0 = 0;
par0 = [Ks0; cms0; cmp0];

theta0 = struct;
theta0.Ks = Ks0;
theta0.cms = cms0;
theta0.cmp = cmp0;

theta0.time = 0;

theta0.cg = 100;
theta0.cads = 30;

thetalb = struct;
thetalb.Ks = 1E-3;
thetalb.cms = 10;
thetalb.cmp = 0;

thetalb.time = 0;

thetalb.cg = 0;
thetalb.cads = 0;

thetaub = struct;
thetaub.Ks = 5;
thetaub.cms = 200;
thetaub.cmp = 100;

est_ind = 1:3;

% use a dummy differential state to measure time
algmodel = struct;
algmodel.nlp_solver_options.ipopt.linear_solver = 'ma27';
algmodel.x = {'time'};
algmodel.z = {'cg', 'cads'};
algmodel.p = {'Ks', 'cms', 'cmp'};
algmodel.d = {'cgmeas', 'cadsmeas'};

algmodel.ode = @(t, y, p) {1};
algmodel.alg = @(t, y, p) {y.cads - p.cmp - ...
			   (p.cms*sqrt(p.Ks)*sqrt(y.cgmeas))/(1 + sqrt(p.Ks)*sqrt(y.cgmeas)), ...
			   y.cg - y.cgmeas};
##algmodel.lsq = @(t, y, p) {(y.cgmeas-y.cg)/10, y.cadsmeas-y.cads};
algmodel.lsq = @(t, y, p) {y.cadsmeas-y.cads};

algmodel.tout = 1:numel(cgmeas);
		      
pe = paresto(algmodel);

est = pe.optimize([cgmeas'; cadsmeas'], theta0, thetalb, thetaub);
theta_opt =  est.theta(est_ind)

theta_conf =  pe.confidence(est, est_ind, 0.95)

table1 = [Hgas1 Hads1];
table2 = [Hgas2 Hads2];
table3 = [Hgas3 Hads3];
table4 = [Hgas4 Hads4];
table5 = [Hgas5 Hads5];
table6 = [Hgas6 Hads6];
table7 = [Hgas7 Hads7];
table8 = [Hgas8 Hads8];

nplot = 100;
Xplot = linspace(0,1.05*max(cgmeas),nplot)';
Ks = est.theta(1);
cms = est.theta(2);
cmx = est.theta(3);
Yplot = cms*sqrt(Ks)*sqrt(Xplot)./(1+sqrt(Ks)*sqrt(Xplot)) + cmx;

table_all = [Xplot, Yplot];

save adsall.dat table_all table1 table2  table3  table4  table5  table6  table7  table8;

if (~ strcmp (getenv ('OMIT_PLOTS'), 'true')) %% PLOTTING
plot (table_all(:,1), table_all(:,2), ...
       table1(:,1), table1(:,2), '+', ...
       table2(:,1), table2(:,2), 'x', ...
       table3(:,1), table3(:,2), '*', ...
       table4(:,1), table4(:,2), 'o', ...
       table5(:,1), table5(:,2), '+', ...
       table6(:,1), table6(:,2), 'x', ...
       table7(:,1), table7(:,2), '*', ...
       table8(:,1), table8(:,2), 'o');
%% TITLE
endif %% PLOTTING
