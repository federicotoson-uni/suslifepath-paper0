%%   SIMPLIFIED ALGORITHM FOR COLLISION RISK ASSESSMENT - Paper 0          %%
%%   Figure generator: Fig.1 altitude distribution, Fig.2 asymmetry        %%
% ----------------------------------------------------------------------- %
% Author: Federico Toson                                                  %
% Produces the two data-driven figures of Paper 0 as vector PDFs.         %
% Reuses the SAME model/inputs as paper0_casestudies.m, so the figures are %
% guaranteed consistent with Table 2 (results).                           %
% ----------------------------------------------------------------------- %
clear; close all; clc;

% --- Auto path discovery (relative to this script's location) --------- %
SCRIPT_DIR = fileparts(mfilename('fullpath'));
REPO_ROOT  = fileparts(SCRIPT_DIR);
DATA_DIR   = fullfile(REPO_ROOT, 'data');
FIG_DIR    = fullfile(REPO_ROOT, 'figures');
addpath(SCRIPT_DIR);
if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

% --- Real catalogue (active + tracked debris) -------------------------- %
CSV    = fullfile(DATA_DIR, 'celestrak_active.csv');
DEBCSV = fullfile(DATA_DIR, 'celestrak_debris.csv');
PopulationData = databasecreator_real(CSV, Inf, DEBCSV);

% --- Case definitions (identical to the driver) ------------------------ %
cases = struct();
cases(1).name = 'ENVISAT';
cases(1).mission = struct('r',7148,'lt',10, 'dt',200,'cost',2.3e9);
cases(1).sc      = struct('ex_surf',100,'tot_surf',250);
cases(2).name = 'Sentinel-6';
cases(2).mission = struct('r',7714,'lt',5.5,'dt',25, 'cost',4.0e8);
cases(2).sc      = struct('ex_surf',10, 'tot_surf',25);
cases(3).name = 'Starlink V2';
cases(3).mission = struct('r',6903,'lt',5,  'dt',5,  'cost',1.0e6);
cases(3).sc      = struct('ex_surf',12, 'tot_surf',30);

% --- Fixed reference mission (700 km SSO smallsat) --------------------- %
ref_m  = struct('r',7078,'lt',7,'dt',25,'cost',1e8);
ref_sc = struct('ex_surf',5,'tot_surf',15);
S_ref  = ref_sc.tot_surf;

[ref_Pind,~] = individual_probability_flux(PopulationData, ref_m, ref_sc);
[ref_Pcol,~] = collective_probability(PopulationData, ref_m, ref_sc, S_ref);
R_ref_ind = ref_Pind*ref_m.cost; if R_ref_ind<=0, R_ref_ind = eps; end
R_ref_col = ref_Pcol;            if R_ref_col<=0, R_ref_col = eps; end

% --- Normalised indices (these reproduce Table 2) ---------------------- %
N = numel(cases); Rt_ind = zeros(1,N); Rt_col = zeros(1,N);
for i = 1:N
    m = cases(i).mission; sc = cases(i).sc;
    [P_ind,~] = individual_probability_flux(PopulationData, m, sc);
    [P_col,~] = collective_probability(PopulationData, m, sc, S_ref);
    Rt_ind(i) = (P_ind*m.cost)/R_ref_ind;
    Rt_col(i) = P_col/R_ref_col;
end

% --- Colours (colour-blind friendly) ----------------------------------- %
col_blue = [0.00 0.45 0.74];     % active / individual
col_oran = [0.85 0.33 0.10];     % debris / collective
col_grey = [0.45 0.45 0.45];

% ====================================================================== %
%  FIGURE 1 - Perigee-altitude distribution
% ====================================================================== %
rearth = 6378;
ALLact = [PopulationData.LEO; PopulationData.GEO; PopulationData.OTH];
h_act  = ALLact(:,1) - rearth;
h_deb  = PopulationData.DEB(:,1) - rearth;
edges  = 0:25:2000;

f1 = figure('Visible','off','Units','centimeters','Position',[2 2 12 7]);
histogram(h_act, edges, 'FaceColor',col_blue,'EdgeColor','none','FaceAlpha',0.85); hold on;
histogram(h_deb, edges, 'FaceColor',col_oran,'EdgeColor','none','FaceAlpha',0.85);
set(gca,'YScale','log');           % populations span orders of magnitude
xlim([0 2000]);
% Extend ylim upward to leave headroom for the case-study labels above the
% data peaks, so the labels never overlap either the histogram bars or
% the legend.
ylim([0.7, 5e4]);

mk    = [525 700 770 1336];
lbl   = {'Starlink (525)','Reference (700)','ENVISAT (770)','Sentinel-6 (1336)'};
% Place all four labels above the data peaks (max bar ~4.5e3) on a single
% horizontal row at y = 1.5e4, rotated 90 deg. With 'VerticalAlignment'
% set to 'bottom', the text grows UPWARD from this baseline, staying in
% the headroom region [1.5e4, 5e4]. Bars never exceed 5e3 in this band,
% so no overlap with data is possible.
y_label = 1.5e4;
x_off   = 22;   % horizontal offset to the right of each dashed line
for k = 1:numel(mk)
    xline(mk(k),'LineStyle','--','Color',col_grey,'LineWidth',0.8);
    text(mk(k)+x_off, y_label, lbl{k}, 'Rotation',90, ...
         'FontSize',14,'Color',[0.25 0.25 0.25], ...
         'VerticalAlignment','bottom','HorizontalAlignment','left');
end
xlabel('Perigee altitude (km)','FontSize',14);
ylabel('Objects per 25 km bin','FontSize',14);
legend({'Active satellites','Tracked debris'},'Box','off', ...
       'Location','southeast','FontSize',13);
set(gca,'FontSize',13,'Box','on','Layer','top');
exportgraphics(f1, fullfile(FIG_DIR,'fig1_altitude_distribution.pdf'),'ContentType','vector');

% ====================================================================== %
%  FIGURE 2 - Individual vs collective asymmetry
% ====================================================================== %
names = {'ENVISAT','Sentinel-6','Starlink V2'};
f2 = figure('Visible','off','Units','centimeters','Position',[2 2 12 7]);
Y  = [Rt_ind; Rt_col]';            % 3x2
b  = bar(Y,'grouped','BarWidth',0.9); hold on;
b(1).FaceColor = col_blue; b(1).EdgeColor = 'none';
b(2).FaceColor = col_oran; b(2).EdgeColor = 'none';
set(gca,'YScale','log');
thr   = [1e-2 1e-1 1e0 1e1];
clab  = {'L','M','H','VH'};
for t = 1:numel(thr)
    yline(thr(t),'LineStyle',':','Color',col_grey,'LineWidth',0.6);
    text(3.55, thr(t), clab{t}, 'FontSize',14,'Color',col_grey, ...
         'VerticalAlignment','middle');
end
ylim([1e-3 1e4]);
set(gca,'XTickLabel',names,'FontSize',13,'Box','on','Layer','top');
xlabel('Case study','FontSize',14);
ylabel('Normalised risk index $\tilde{R}$','Interpreter','latex','FontSize',14);
legend({'Individual','Collective'},'Box','off','Location','northwest','FontSize',13);
exportgraphics(f2, fullfile(FIG_DIR,'fig2_asymmetry.pdf'),'ContentType','vector');

% --- Console echo (sanity check vs Table 2) ---------------------------- %
fprintf('\nFigures written to:\n  %s\n', FIG_DIR);
fprintf('Rt_ind = %8.3f %8.3f %8.3f\n', Rt_ind);
fprintf('Rt_col = %8.3f %8.3f %8.3f\n', Rt_col);
fprintf('(expected ~ ind 1018 / 0.02 / 0.09 ; col 903 / 0.05 / 29.95)\n');
