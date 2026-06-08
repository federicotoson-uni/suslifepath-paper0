%%   SIMPLIFIED ALGORITHM FOR COLLISION RISK ASSESSMENT - Paper 0         %%
%%   SNAPSHOT-ECOB PROXY BENCHMARK on case studies + N=100 sample         %%
% ----------------------------------------------------------------------- %
% Author: Federico Toson                                                  %
% Driver for Figure 7 + Table 5 of Paper 0 (revision after first review)  %
% ----------------------------------------------------------------------- %
%{
  Computes the snapshot-ECOB proxy of ecob_proxy.m on (i) the three case
  studies of Section 4.1 and (ii) the random N=100 sample of Section 4.2,
  then ranks both indicators (our normalised collective $\tilde{R}^{col}$
  and the snapshot-ECOB proxy) and reports the Spearman rank correlation
  on the joint sample. The goal is INTERNAL BENCHMARKING: showing that
  the rapid screening produced by P_col reproduces the ECOB-style ranking
  for the catalogue snapshot of May 2026, without claiming numerical
  equivalence with the full ECOB toolchain (which requires 200-yr
  debris-evolution simulation and is out of scope for this paper).

  MASS ASSIGNMENT:
    * Case studies use published / engineering-estimate values:
        ENVISAT          : 8211 kg  (ESA fact sheet)
        Sentinel-6       : 1192 kg  (NASA JPL fact sheet)
        Starlink V2 Mini :  800 kg  (SpaceX engineering value)
        Reference smallsat: 500 kg
    * N=100 sample: stratified log-normal mass prior matching the cost prior
      tiers (constellation 500 kg, EO/science 1500 kg, specialty 3000 kg).

  Same RNG seed as paper0_validation.m (rng(20260605, 'twister')) for
  byte-identical reproduction.

  OUTPUT:
    figures/Figure_7_ecob_benchmark.pdf  - scatter Rcol_norm vs ECOB-proxy_norm
                                           with Spearman rho on log-log axes
    Console: 3-case study Tab 5 + Spearman over the joint sample (3 cases + 100)
%}
% ----------------------------------------------------------------------- %
clear; close all; clc;

% --- Auto path discovery ---------------------------------------------- %
THIS_DIR = fileparts(mfilename('fullpath'));
DATA_DIR = fullfile(THIS_DIR, '..', 'data');
FIG_DIR  = fullfile(THIS_DIR, '..', 'figures');
if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end
addpath(THIS_DIR);

CSV    = fullfile(DATA_DIR, 'celestrak_active.csv');
DEBCSV = fullfile(DATA_DIR, 'celestrak_debris.csv');

% --- Reproducibility seed (SAME as paper0_validation.m) --------------- %
rng(20260605, 'twister');

% --- Build full population once --------------------------------------- %
PopulationData = databasecreator_real(CSV, Inf, DEBCSV);

% --- Re-read names + elements ----------------------------------------- %
T       = readtable(CSV);
mu      = 398600.4418;
conv    = pi/180;
n_rev   = T.MEAN_MOTION;
ecc_all = T.ECCENTRICITY;
inc_all = T.INCLINATION;
w_all   = T.ARG_OF_PERICENTER * conv;
ra_all  = T.RA_OF_ASC_NODE    * conv;
n       = n_rev * 2*pi / 86400;
a_all   = (mu ./ n.^2).^(1/3);
rp_all  = a_all .* (1 - ecc_all);
alt_all = rp_all - 6378;
name_all= string(T.OBJECT_NAME);
isLEO   = isfinite(alt_all) & alt_all > 0 & alt_all < 2000;
idx_LEO = find(isLEO);

% --- Reference mission (SAME as everywhere else) ---------------------- %
ref_m  = struct('yol',2020,'lt',7,'dt',25,'it',0.2,'r',7078, ...
    'ecc',0.001,'w',0,'ra',0,'inc',98,'pl','G','net',0,'man','U','cost',1e8);
ref_sc = struct('ex_surf',5,'tot_surf',15,'rho',500);
S_ref  = ref_sc.tot_surf;
M_ref  = 500;                          % [kg] reference smallsat
[ref_Pind, ~] = individual_probability_flux(PopulationData, ref_m, ref_sc);
[ref_Pcol, ~] = collective_probability(PopulationData, ref_m, ref_sc, S_ref);
[ref_Peco, ~] = ecob_proxy(PopulationData, ref_m, ref_sc, M_ref, M_ref);
R_ref_ind = ref_Pind * ref_m.cost;   if R_ref_ind<=0, R_ref_ind = eps; end
R_ref_col = ref_Pcol;                if R_ref_col<=0, R_ref_col = eps; end
R_ref_eco = ref_Peco;                if R_ref_eco<=0, R_ref_eco = eps; end
fprintf('R_ref_ind = %.3e | R_ref_col = %.3e | R_ref_eco = %.3e\n\n', ...
        R_ref_ind, R_ref_col, R_ref_eco);

% ====================================================================== %
%  THREE CASE STUDIES (with published / engineering-estimate masses)
% ====================================================================== %
case_names = {'ENVISAT','Sentinel-6','Starlink V2 Mini'};
case_alts  = [770, 1336, 525];
case_rps   = case_alts + 6378;
case_costs = [2.3e9, 4.0e8, 1.0e6];
case_lt    = [10, 5.5, 5];
case_dt    = [200, 25, 5];
case_ex    = [100, 10, 12];
case_tot   = [250, 25, 30];
case_inc   = [98.4, 66, 53];
case_mass  = [8211, 1192, 800];       % [kg] published / SpaceX engineering value

case_Cind  = zeros(1,3);
case_Ccol  = zeros(1,3);
case_Ceco  = zeros(1,3);
for c = 1:3
    m  = struct('yol',2020,'lt',case_lt(c),'dt',case_dt(c),'it',0.2, ...
                'r',case_rps(c),'ecc',0.001,'w',0,'ra',0, ...
                'inc',case_inc(c),'pl','G','net',0,'man','U','cost',case_costs(c));
    sc = struct('ex_surf',case_ex(c),'tot_surf',case_tot(c),'rho',500);
    [Pi, ~] = individual_probability_flux(PopulationData, m, sc);
    [Pc, ~] = collective_probability(PopulationData, m, sc, S_ref);
    [Pe, ~] = ecob_proxy(PopulationData, m, sc, M_ref, case_mass(c));
    case_Cind(c) = (Pi * case_costs(c)) / R_ref_ind;
    case_Ccol(c) = Pc / R_ref_col;
    case_Ceco(c) = Pe / R_ref_eco;
end

fprintf('===== ECOB-PROXY BENCHMARK on three case studies (Table 5) =====\n');
fprintf('%-18s %8s %12s %12s %12s\n','Case','M [kg]','C_ind','C_col','C_eco');
fprintf('%s\n', repmat('-',1,68));
for c = 1:3
    fprintf('%-18s %8.0f %12.3f %12.3f %12.3f\n', ...
            case_names{c}, case_mass(c), case_Cind(c), case_Ccol(c), case_Ceco(c));
end
fprintf('\n');

% Rank-order check on the three cases (Letizia 2016 places ENVISAT as the
% upper baseline; Letizia 2019 places sparse-band operational missions and
% per-element constellation contributions in the lower classes).
[~, rk_col] = sort(-case_Ccol);   [~, rk_col_pos(rk_col)] = sort(1:3);  %#ok<NASGU>
[~, rk_eco] = sort(-case_Ceco);   [~, rk_eco_pos(rk_eco)] = sort(1:3);  %#ok<NASGU>
fprintf('Ranking by C_col : '); fprintf('%-18s ', case_names{rk_col}); fprintf('\n');
fprintf('Ranking by C_eco : '); fprintf('%-18s ', case_names{rk_eco}); fprintf('\n');
fprintf('ECOB literature  : %-18s %-18s %-18s (Letizia 2016 ENVISAT baseline)\n\n', ...
        case_names{1}, case_names{3}, case_names{2});

% ====================================================================== %
%  N=100 RANDOM SAMPLE (same seed and sample as paper0_validation.m)
% ====================================================================== %
N = 100;
samp_idx = randperm(numel(idx_LEO), N);
sel      = idx_LEO(samp_idx);

names   = strings(N,1);
alts    = zeros(N,1);
costs   = zeros(N,1);
masses  = zeros(N,1);
exsurfs = zeros(N,1);
tsurfs  = zeros(N,1);
Pind_a  = zeros(N,1);
Pcol_a  = zeros(N,1);
Peco_a  = zeros(N,1);
Cind_a  = zeros(N,1);
Ccol_a  = zeros(N,1);
Ceco_a  = zeros(N,1);
tiers   = strings(N,1);

for i = 1:N
    k   = sel(i);
    h_k = alt_all(k);
    [cost, ex_surf, tot_surf, lt, dt, mass_kg] = tier_params(h_k);

    m  = struct('yol',2020,'lt',lt,'dt',dt,'it',0.2, ...
                'r',rp_all(k),'ecc',ecc_all(k),'w',w_all(k),'ra',ra_all(k), ...
                'inc',inc_all(k),'pl','G','net',0,'man','U','cost',cost);
    sc = struct('ex_surf',ex_surf,'tot_surf',tot_surf,'rho',500);

    [Pind, ~] = individual_probability_flux(PopulationData, m, sc);
    [Pcol, ~] = collective_probability(PopulationData, m, sc, S_ref);
    [Peco, ~] = ecob_proxy(PopulationData, m, sc, M_ref, mass_kg);

    Cind_a(i) = (Pind * cost) / R_ref_ind;
    Ccol_a(i) = Pcol         / R_ref_col;
    Ceco_a(i) = Peco         / R_ref_eco;

    if h_k < 600
        tier = "constellation";
    elseif h_k < 1500
        tier = "EO/science";
    else
        tier = "specialty";
    end

    names(i)   = name_all(k);
    alts(i)    = h_k;
    costs(i)   = cost;
    masses(i)  = mass_kg;
    exsurfs(i) = ex_surf;
    tsurfs(i)  = tot_surf;
    Pind_a(i)  = Pind;
    Pcol_a(i)  = Pcol;
    Peco_a(i)  = Peco;
    tiers(i)   = tier;
end

% --- Spearman rank correlation ---------------------------------------- %
% Joint vector: 3 case studies + 100 sample = 103 points.
joint_col = [case_Ccol(:); Ccol_a(:)];
joint_eco = [case_Ceco(:); Ceco_a(:)];
[rho_S, pval] = corr(joint_col, joint_eco, 'type','Spearman');
[rho_P, ~]    = corr(log10(joint_col), log10(joint_eco), 'type','Pearson');
fprintf('===== RANK CORRELATION C_col vs C_eco (joint sample n=%d) =====\n', ...
        numel(joint_col));
fprintf('  Spearman rho = %.4f  (p = %.2e)\n', rho_S, pval);
fprintf('  Pearson  rho = %.4f  on log10 scale\n', rho_P);
fprintf('\n');

% Per-tier Spearman on the sample alone
for tl = ["constellation","EO/science"]
    mask = tiers == tl;
    if sum(mask) >= 3
        [rho_t, ~] = corr(Ccol_a(mask), Ceco_a(mask), 'type','Spearman');
        fprintf('  Spearman (tier %s, n=%d) = %.4f\n', tl, sum(mask), rho_t);
    end
end
fprintf('\n');

% --- LaTeX-ready Tab 5 snippet ---------------------------------------- %
fprintf('===== LaTeX SNIPPET (Table 5 - case studies) =====\n');
fprintf('\\begin{table}[ht]\n');
fprintf('\\centering\n');
fprintf('\\footnotesize\n');
fprintf('\\setlength{\\tabcolsep}{6pt}\n');
fprintf('\\caption{Snapshot-ECOB proxy benchmark on the three case studies. $\\tilde{R}^{col}$ is the normalised collective index of this paper (area-based fragmentation weight, operational+residual horizon); $\\tilde{R}^{eco}$ is the snapshot-ECOB proxy defined in Section~\\ref{sec:ecob-bench} (mass-based NASA-Standard-Breakup fragmentation weight, 200-year horizon, same catalogue snapshot). Both indices are normalised against the reference mission (\\SI{500}{kg}, \\SI{700}{km} sun-synchronous). Object masses are published / engineering-estimate values: ENVISAT \\SI{8211}{kg} (ESA), Sentinel-6 \\SI{1192}{kg} (NASA JPL), Starlink V2 Mini \\SI{800}{kg} (SpaceX). The two indices produce the same qualitative ranking as the ECOB literature baseline (Letizia 2016, with ENVISAT as the upper reference).}\n');
fprintf('\\label{tab:ecob-cases}\n');
fprintf('\\begin{tabular}{lrrrr}\n');
fprintf('\\toprule\n');
fprintf('Case study & $M$ (kg) & $\\tilde{R}^{col}$ & $\\tilde{R}^{eco}$ & ratio $\\tilde{R}^{eco}/\\tilde{R}^{col}$ \\\\ \n');
fprintf('\\midrule\n');
for c = 1:3
    fprintf('%s & %.0f & %.2f & %.2f & %.2f \\\\ \n', ...
            case_names{c}, case_mass(c), case_Ccol(c), case_Ceco(c), ...
            case_Ceco(c)/case_Ccol(c));
end
fprintf('\\bottomrule\n');
fprintf('\\end{tabular}\n');
fprintf('\\end{table}\n\n');

% ====================================================================== %
%  FIGURE 7: scatter R_col vs R_eco (log-log), Spearman in legend
% ====================================================================== %
fig = figure('Color','w','Units','centimeters','Position',[2 2 13 11]);
hold on; box on; grid on;
set(gca,'XScale','log','YScale','log');

% Sample, coloured by tier
tier_colors = containers.Map( ...
    {'constellation','EO/science','specialty'}, ...
    {[0.30 0.45 0.75], [0.10 0.60 0.30], [0.85 0.50 0.00]});
h_legend = gobjects(0);
legend_labels = {};
for tl = ["constellation","EO/science"]
    mask = tiers == tl;
    if any(mask)
        cc = tier_colors(char(tl));
        h = scatter(Ccol_a(mask), Ceco_a(mask), 42, ...
                    'MarkerFaceColor', cc, 'MarkerEdgeColor','w', ...
                    'MarkerFaceAlpha', 0.85);
        h_legend(end+1) = h;
        legend_labels{end+1} = char(tl);
    end
end

% Case studies, diamonds
case_colors = [0.7 0.1 0.1; 0.2 0.6 0.2; 0.85 0.5 0.0];
for c = 1:3
    hc = plot(case_Ccol(c), case_Ceco(c), 'd', ...
              'MarkerFaceColor', case_colors(c,:), 'MarkerEdgeColor','k', ...
              'MarkerSize', 13, 'LineWidth', 1.2);
    text(case_Ccol(c)*1.25, case_Ceco(c), case_names{c}, ...
         'FontSize',14,'FontWeight','bold','Color',case_colors(c,:));
    if c == 1
        h_legend(end+1) = hc;
        legend_labels{end+1} = 'case studies';
    end
end

% Identity reference line (Rcol = Reco)
xl = [1e-3, 1e4];
plot(xl, xl, 'k:', 'LineWidth', 1.1, 'HandleVisibility','off');
text(1e3, 5e2, '$\tilde{R}^{eco}=\tilde{R}^{col}$', ...
     'Interpreter','latex','FontSize',13,'Color',[.3 .3 .3]);

xlim(xl);
ylim([1e-3, 1e5]);
xlabel('$\tilde{R}^{col}$ (our normalised collective index)','Interpreter','latex','FontSize',14);
ylabel('$\tilde{R}^{eco}$ (snapshot-ECOB proxy)','Interpreter','latex','FontSize',14);
title(sprintf('Internal benchmark: $\\rho_{Spearman} = %.2f$ on $n = %d$ joint sample', ...
              rho_S, numel(joint_col)), 'Interpreter','latex','FontSize',14);
legend(h_legend, legend_labels, 'Location','northwest','FontSize',13);
set(gca,'FontSize',13);

% --- Save vector PDF -------------------------------------------------- %
out_pdf = fullfile(FIG_DIR, 'Figure_7_ecob_benchmark.pdf');
exportgraphics(fig, out_pdf, 'ContentType','vector');
fprintf('Saved: %s\n', out_pdf);

% ----------------------------------------------------------------------- %
%  LOCAL FUNCTIONS (extends tier_params of paper0_validation.m with mass)
% ----------------------------------------------------------------------- %
function [cost, ex_surf, tot_surf, lt, dt, mass_kg] = tier_params(h_km)
% Stratified prior for non-orbital parameters by altitude tier (with mass).
    if h_km < 600
        cost     = lognrnd(log(5e6),  0.5);
        ex_surf  = max(3, 30 + 10*randn);
        tot_surf = ex_surf * 2.0;
        mass_kg  = lognrnd(log(500),  0.4);
        lt       = 5;
        dt       = 5;
    elseif h_km < 1500
        cost     = lognrnd(log(150e6), 0.7);
        ex_surf  = max(5, 50 + 20*randn);
        tot_surf = ex_surf * 2.2;
        mass_kg  = lognrnd(log(1500), 0.6);
        lt       = 7;
        dt       = 25;
    else
        cost     = lognrnd(log(400e6), 0.5);
        ex_surf  = max(8, 80 + 30*randn);
        tot_surf = ex_surf * 2.5;
        mass_kg  = lognrnd(log(3000), 0.5);
        lt       = 10;
        dt       = 25;
    end
end
