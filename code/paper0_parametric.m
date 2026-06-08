%%   SIMPLIFIED ALGORITHM FOR COLLISION RISK ASSESSMENT - Paper 0          %%
%%   Parametric extension: cost-driven individual-collective asymmetry     %%
% ----------------------------------------------------------------------- %
% Author: Federico Toson                                                  %
% Generates Figure 4 and the cost-threshold table for paper Section 4.2.  %
%                                                                         %
% Mathematical core (cf. Section 4.2):                                    %
%   With sigma, A, T_op, T_res held at the reference values, the          %
%   asymmetry between normalised collective and individual indices is     %
%       R_col_tilde / R_ind_tilde (h, C) = K(h) / C                       %
%   K(h) = C_ref * [P_col(h)/P_col_ref] * [P_ind_ref/P_ind(h)]            %
%   so a satellite at altitude h with replacement cost below              %
%       C*(h, A) = K(h) / A                                               %
%   is in the regime where the individual-only assessment understates the %
%   collective burden by at least a factor A (e.g. A = 100 -> two orders).%
% ----------------------------------------------------------------------- %
clear; close all; clc;

% Paths auto-discovered relative to this script's location ------------- %
SCRIPT_DIR = fileparts(mfilename('fullpath'));
REPO_ROOT  = fileparts(SCRIPT_DIR);
DATA_DIR   = fullfile(REPO_ROOT, 'data');
FIG_DIR    = fullfile(REPO_ROOT, 'figures');
addpath(SCRIPT_DIR);
if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

% Population (same as the driver) --------------------------------------- %
CSV    = fullfile(DATA_DIR,'celestrak_active.csv');
DEBCSV = fullfile(DATA_DIR,'celestrak_debris.csv');
PopulationData = databasecreator_real(CSV, Inf, DEBCSV);

% Reference mission (700 km SSO, identical to paper) -------------------- %
ref_m  = struct('r',7078,'lt',7,'dt',25,'cost',1e8);
ref_sc = struct('ex_surf',5,'tot_surf',15);
S_ref  = ref_sc.tot_surf;
[ref_Pind,~] = individual_probability_flux(PopulationData, ref_m, ref_sc);
[ref_Pcol,~] = collective_probability(PopulationData, ref_m, ref_sc, S_ref);
R_ref_ind = ref_Pind * ref_m.cost;
R_ref_col = ref_Pcol;

rearth = 6378;
alts   = [525, 700, 770, 1336];
labels = {'525 km (Starlink shell)', '700 km (reference SSO)', ...
          '770 km (ENVISAT / debris band)', '1336 km (Sentinel-6 band)'};

% P_ind, P_col at each altitude, with REFERENCE geometry & lifetime ----- %
P_ind_h = zeros(numel(alts),1);
P_col_h = zeros(numel(alts),1);
for k = 1:numel(alts)
    m = struct('r', alts(k)+rearth, 'lt', ref_m.lt, 'dt', ref_m.dt, 'cost', 0);
    [P_ind_h(k), ~] = individual_probability_flux(PopulationData, m, ref_sc);
    [P_col_h(k), ~] = collective_probability(PopulationData, m, ref_sc, S_ref);
end

% K(h) = C_ref * (P_col(h)/P_col_ref) * (P_ind_ref/P_ind(h))
K_h = ref_m.cost .* (P_col_h ./ ref_Pcol) ./ (P_ind_h ./ ref_Pind);

% Cost sweep and asymmetry surface
costs = logspace(4, 10, 200);
asym  = K_h ./ costs;          % size: numel(alts) x numel(costs)

% Thresholds
C_thresh_100 = K_h / 100;      % 2 orders of magnitude
C_thresh_10  = K_h / 10;       % 1 order
C_thresh_1   = K_h / 1;        % crossover (col index = ind index)

% Console echo (for paper table) ---------------------------------------- %
fprintf('\n========== K(h) and cost thresholds ==========\n');
fprintf('Altitude (km) |    K(h)     | C*(asym=10) | C*(asym=100)\n');
fprintf('--------------+-------------+-------------+-------------\n');
for k = 1:numel(alts)
    fprintf('  %4d        |  %.3e  |  %.2e   |  %.2e\n', ...
            alts(k), K_h(k), C_thresh_10(k), C_thresh_100(k));
end
fprintf('==============================================\n');
fprintf('R_ref_ind = %.3e   R_ref_col = %.3e\n', R_ref_ind, R_ref_col);
fprintf('P_ind_h = %s\n', mat2str(P_ind_h',3));
fprintf('P_col_h = %s\n', mat2str(P_col_h',3));

% ====================================================================== %
% PLOT - Figure 4
% ====================================================================== %
col_blue = [0.00 0.45 0.74];
col_yell = [0.93 0.69 0.13];
col_oran = [0.85 0.33 0.10];
col_purp = [0.49 0.18 0.56];
cmap = [col_blue; col_yell; col_oran; col_purp];

f = figure('Visible','off','Units','centimeters','Position',[2 2 16 10]);

% Distinct line styles + markers so that visually overlapping pairs
% (525/1336 and 700/770 km) remain individually identifiable.
linestyles = {'-', '--', '-.', ':'};
markers    = {'o', 's', '^', 'd'};
nC         = numel(costs);
% Stagger marker x-positions so they don't pile up at the same point
mkri       = {15:30:nC, 22:30:nC, 8:30:nC, 28:30:nC};

% main curves
for k = 1:numel(alts)
    loglog(costs, asym(k,:), linestyles{k}, 'Color', cmap(k,:), 'LineWidth', 1.5, ...
           'Marker', markers{k}, 'MarkerIndices', mkri{k}, ...
           'MarkerFaceColor', cmap(k,:), 'MarkerEdgeColor', 'none', ...
           'MarkerSize', 4.5, 'DisplayName', labels{k}); hold on;
end

% asymmetry guide lines
yline(100, 'LineStyle','--','Color',[0.30 0.30 0.30],'LineWidth',0.8, ...
      'HandleVisibility','off');
yline(10,  'LineStyle',':', 'Color',[0.45 0.45 0.45],'LineWidth',0.6, ...
      'HandleVisibility','off');
yline(1,   'LineStyle','-', 'Color',[0.65 0.65 0.65],'LineWidth',0.4, ...
      'HandleVisibility','off');

text(5e9, 140, 'asymmetry = 100 (two orders)', ...
     'FontSize',13,'Color',[0.30 0.30 0.30], ...
     'HorizontalAlignment','right');
text(5e9, 14,  'asymmetry = 10', ...
     'FontSize',13,'Color',[0.45 0.45 0.45], ...
     'HorizontalAlignment','right');

xlabel('Replacement cost $C$ (USD)','Interpreter','latex','FontSize',13);
ylabel('Asymmetry $\tilde{R}^{col}/\tilde{R}^{ind}$','Interpreter','latex','FontSize',13);
legend('show','Location','northeast','FontSize',11,'Box','off');
ylim([1e-3, 1e7]);
xlim([1e4, 1e10]);
grid on;
set(gca,'FontSize',12,'Box','on','GridAlpha',0.15,'Layer','top');

exportgraphics(f, fullfile(FIG_DIR,'fig4_parametric_asymmetry.pdf'), ...
               'ContentType','vector');
fprintf('Figure written to: %s\n', ...
        fullfile(FIG_DIR,'fig4_parametric_asymmetry.pdf'));
