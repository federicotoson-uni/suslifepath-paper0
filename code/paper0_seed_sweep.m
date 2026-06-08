%%   SIMPLIFIED ALGORITHM FOR COLLISION RISK ASSESSMENT - Paper 0         %%
%%   SEED SWEEP for statistical stability of validation metrics          %%
% ----------------------------------------------------------------------- %
% Author: Federico Toson                                                  %
% Defuses reviewer attack "single seed, no sweep"                         %
% ----------------------------------------------------------------------- %
%{
  Reruns paper0_validation.m logic across K independent random seeds
  (seeds 20260605 + 0..K-1) and reports the range of:
    - median R (asymmetry)
    - p99 R
    - P(R>1)
    - Spearman rho between Ccol and Ceco on the joint sample
  This addresses the methodology reviewer attack that all reported numbers
  are conditional on a single rng draw.

  OUTPUT: console table + LaTeX snippet for an extra row in Tab 6.
%}
% ----------------------------------------------------------------------- %
clear; close all; clc;

THIS_DIR = fileparts(mfilename('fullpath'));
DATA_DIR = fullfile(THIS_DIR, '..', 'data');
addpath(THIS_DIR);

CSV    = fullfile(DATA_DIR, 'celestrak_active.csv');
DEBCSV = fullfile(DATA_DIR, 'celestrak_debris.csv');

PopulationData = databasecreator_real(CSV, Inf, DEBCSV);

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
isLEO   = isfinite(alt_all) & alt_all > 0 & alt_all < 2000;
idx_LEO = find(isLEO);

% Reference
ref_m  = struct('yol',2020,'lt',7,'dt',25,'it',0.2,'r',7078, ...
    'ecc',0.001,'w',0,'ra',0,'inc',98,'pl','G','net',0,'man','U','cost',1e8);
ref_sc = struct('ex_surf',5,'tot_surf',15,'rho',500);
S_ref  = ref_sc.tot_surf;
M_ref  = 500;
[ref_Pind, ~] = individual_probability_flux(PopulationData, ref_m, ref_sc);
[ref_Pcol, ~] = collective_probability(PopulationData, ref_m, ref_sc, S_ref);
[ref_Peco, ~] = ecob_proxy(PopulationData, ref_m, ref_sc, M_ref, M_ref);
R_ref_ind = max(ref_Pind * ref_m.cost, eps);
R_ref_col = max(ref_Pcol, eps);
R_ref_eco = max(ref_Peco, eps);

% Case studies for Spearman joint sample
case_alts  = [770, 1336, 525];
case_rps   = case_alts + 6378;
case_costs = [2.3e9, 4.0e8, 1.0e6];
case_lt    = [10, 5.5, 5];
case_dt    = [200, 25, 5];
case_ex    = [100, 10, 12];
case_tot   = [250, 25, 30];
case_mass  = [8211, 1192, 800];
case_inc   = [98.4, 66, 53];
case_Ccol_fixed = zeros(1,3);
case_Ceco_fixed = zeros(1,3);
for c = 1:3
    m  = struct('yol',2020,'lt',case_lt(c),'dt',case_dt(c),'it',0.2, ...
                'r',case_rps(c),'ecc',0.001,'w',0,'ra',0, ...
                'inc',case_inc(c),'pl','G','net',0,'man','U','cost',case_costs(c));
    sc = struct('ex_surf',case_ex(c),'tot_surf',case_tot(c),'rho',500);
    [Pc, ~] = collective_probability(PopulationData, m, sc, S_ref);
    [Pe, ~] = ecob_proxy(PopulationData, m, sc, M_ref, case_mass(c));
    case_Ccol_fixed(c) = Pc / R_ref_col;
    case_Ceco_fixed(c) = Pe / R_ref_eco;
end

% --- Seed sweep -------------------------------------------------------- %
K = 100;
seeds = 20260605 + (0:K-1);

med_R   = zeros(K,1);
p99_R   = zeros(K,1);
PR_gt_1 = zeros(K,1);
rho_S   = zeros(K,1);
N = 100;

for s = 1:K
    rng(seeds(s), 'twister');
    samp = randperm(numel(idx_LEO), N);
    sel  = idx_LEO(samp);

    Cind = zeros(N,1);  Ccol = zeros(N,1);  Ceco = zeros(N,1);
    for i = 1:N
        k = sel(i);
        h_k = alt_all(k);
        [cost, ex_surf, tot_surf, lt, dt, mass_kg] = tier_params(h_k);
        m  = struct('yol',2020,'lt',lt,'dt',dt,'it',0.2, ...
                    'r',rp_all(k),'ecc',ecc_all(k),'w',w_all(k),'ra',ra_all(k), ...
                    'inc',inc_all(k),'pl','G','net',0,'man','U','cost',cost);
        sc = struct('ex_surf',ex_surf,'tot_surf',tot_surf,'rho',500);
        [Pind, ~] = individual_probability_flux(PopulationData, m, sc);
        [Pcol, ~] = collective_probability(PopulationData, m, sc, S_ref);
        [Peco, ~] = ecob_proxy(PopulationData, m, sc, M_ref, mass_kg);
        Cind(i) = (Pind * cost) / R_ref_ind;
        Ccol(i) = Pcol / R_ref_col;
        Ceco(i) = Peco / R_ref_eco;
    end

    asym = Ccol ./ max(Cind, eps);
    med_R(s)   = median(asym);
    p99_R(s)   = prctile(asym, 99);
    PR_gt_1(s) = mean(asym > 1);

    joint_col = [case_Ccol_fixed(:); Ccol(:)];
    joint_eco = [case_Ceco_fixed(:); Ceco(:)];
    rho_S(s)  = corr(joint_col, joint_eco, 'type','Spearman');

    if mod(s, 10) == 0
        fprintf('  seed %d/%d done\n', s, K);
    end
end

% --- Summary ---------------------------------------------------------- %
fprintf('\n===== SEED SWEEP RESULTS (K=%d, N=%d per seed) =====\n', K, N);
fprintf('%-18s %10s %10s %10s %10s\n', 'Statistic','median','p05','p95','SD');
fprintf('%s\n', repmat('-', 1, 64));
fprintf('%-18s %10.2f %10.2f %10.2f %10.2f\n', 'median R',     median(med_R),   prctile(med_R,5),   prctile(med_R,95),   std(med_R));
fprintf('%-18s %10.2f %10.2f %10.2f %10.2f\n', 'p99 R',        median(p99_R),   prctile(p99_R,5),   prctile(p99_R,95),   std(p99_R));
fprintf('%-18s %10.3f %10.3f %10.3f %10.3f\n', 'P(R>1)',       median(PR_gt_1), prctile(PR_gt_1,5), prctile(PR_gt_1,95), std(PR_gt_1));
fprintf('%-18s %10.3f %10.3f %10.3f %10.3f\n', 'Spearman rho', median(rho_S),   prctile(rho_S,5),   prctile(rho_S,95),   std(rho_S));

% Min and max for honest range
fprintf('\nMin/max across the %d seeds:\n', K);
fprintf('  median R   : [%.1f, %.1f]\n', min(med_R), max(med_R));
fprintf('  p99 R      : [%.1f, %.1f]\n', min(p99_R), max(p99_R));
fprintf('  P(R>1)     : [%.0f%%, %.0f%%]\n', 100*min(PR_gt_1), 100*max(PR_gt_1));
fprintf('  Spearman   : [%.3f, %.3f]\n', min(rho_S), max(rho_S));

% Clopper-Pearson 95% lower bound for binomial P(R>1) on the original seed
% (paper0_validation seed 20260605, 100/100 successes)
% For x=100, n=100: lower = (alpha/2)^(1/n), alpha = 0.05
lower_cp = (0.025)^(1/100);
fprintf('\nClopper-Pearson 95%% lower bound for 100/100: P(R>1) >= %.4f (%.2f%%)\n', ...
        lower_cp, 100*lower_cp);

% ----------------------------------------------------------------------- %
function [cost, ex_surf, tot_surf, lt, dt, mass_kg] = tier_params(h_km)
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
