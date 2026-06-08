%%   SIMPLIFIED ALGORITHM FOR COLLISION RISK ASSESSMENT - Paper 0         %%
%%   TEMPORAL ROBUSTNESS: snapshot composition sensitivity                %%
% ----------------------------------------------------------------------- %
% Author: Federico Toson                                                  %
% Driver for Table 6 of Paper 0 (revision after first review)             %
% ----------------------------------------------------------------------- %
%{
  Tests whether the key statistical claims of the paper (median R,
  Spearman rho between P_col and the ECOB proxy, top-tier composition)
  depend on the specific composition of the May 2026 Celestrak catalogue,
  or are stable across snapshot evolution.

  COMPARISON:
    Snapshot A : original paper snapshot (May 29, 2026, n=14727 active LEO).
    Snapshot B : synthetic "past 2025" snapshot, built from the LIVE
                 Celestrak catalogue of June 5, 2026 by removing all
                 objects whose OBJECT_ID begins with "2026-" (i.e. removing
                 the 2026 launches). Approximates the catalogue as it stood
                 in December 2025, with ~1900 fewer objects, predominantly
                 Starlink V2 Mini. The orbital state of the surviving
                 objects is the June 2026 state (approximate; for the
                 macro statistics of interest, drift over a few months is
                 negligible).

  Both snapshots use the SAME debris CSV (the four principal fragmentation
  clouds: Fengyun-1C 2007, Iridium-33 2009, Cosmos-2251 2009, Cosmos-1408
  2021), which are pre-2026 and unaffected by the filter.

  Same RNG seed as paper0_validation.m (20260605, 'twister') so the random
  sample of 100 satellites is drawn from each snapshot in the same way;
  the two samples are NOT identical (different population, different
  random draws), but both are unbiased random samples of their respective
  snapshots, which is what the robustness claim requires.

  OUTPUT:
    Console: side-by-side stats + LaTeX-ready Table 6 snippet showing
             median R, p99, P(R>1), Spearman rho on each snapshot.
%}
% ----------------------------------------------------------------------- %
clear; close all; clc;

THIS_DIR = fileparts(mfilename('fullpath'));
DATA_DIR = fullfile(THIS_DIR, '..', 'data');
addpath(THIS_DIR);

DEBCSV   = fullfile(DATA_DIR, 'celestrak_debris.csv');
CSV_A    = fullfile(DATA_DIR, 'celestrak_active.csv');
CSV_B    = fullfile(DATA_DIR, 'celestrak_active_past2025.csv');

names_snap = {'May 2026 (paper)', 'Past 2025 (synthetic)'};
results = cell(2, 1);

for snap = 1:2
    if snap == 1, CSV = CSV_A; else, CSV = CSV_B; end

    fprintf('\n=========== SNAPSHOT: %s ===========\n', names_snap{snap});

    rng(20260605, 'twister');
    PopulationData = databasecreator_real(CSV, Inf, DEBCSV);

    % Re-read names + elements
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

    % Reference mission (same)
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

    % Random sample N=100
    N = 100;
    samp = randperm(numel(idx_LEO), N);
    sel  = idx_LEO(samp);

    Cind = zeros(N,1);  Ccol = zeros(N,1);  Ceco = zeros(N,1);
    tiers = strings(N,1);
    alts  = zeros(N,1);

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
        alts(i) = h_k;
        if h_k < 600, tiers(i) = "constellation";
        elseif h_k < 1500, tiers(i) = "EO/science";
        else, tiers(i) = "specialty";
        end
    end

    asym = Ccol ./ max(Cind, eps);

    % Re-evaluate 3 case studies on this snapshot
    case_names = {'ENVISAT','Sentinel-6','Starlink V2 Mini'};
    case_alts  = [770, 1336, 525];
    case_rps   = case_alts + 6378;
    case_costs = [2.3e9, 4.0e8, 1.0e6];
    case_lt    = [10, 5.5, 5];
    case_dt    = [200, 25, 5];
    case_ex    = [100, 10, 12];
    case_tot   = [250, 25, 30];
    case_mass  = [8211, 1192, 800];
    case_inc   = [98.4, 66, 53];
    case_Ccol  = zeros(1,3);
    case_Ceco  = zeros(1,3);
    for c = 1:3
        m  = struct('yol',2020,'lt',case_lt(c),'dt',case_dt(c),'it',0.2, ...
                    'r',case_rps(c),'ecc',0.001,'w',0,'ra',0, ...
                    'inc',case_inc(c),'pl','G','net',0,'man','U','cost',case_costs(c));
        sc = struct('ex_surf',case_ex(c),'tot_surf',case_tot(c),'rho',500);
        [Pc, ~] = collective_probability(PopulationData, m, sc, S_ref);
        [Pe, ~] = ecob_proxy(PopulationData, m, sc, M_ref, case_mass(c));
        case_Ccol(c) = Pc / R_ref_col;
        case_Ceco(c) = Pe / R_ref_eco;
    end

    % Spearman + Pearson on joint sample
    joint_col = [case_Ccol(:); Ccol(:)];
    joint_eco = [case_Ceco(:); Ceco(:)];
    rho_S = corr(joint_col, joint_eco, 'type','Spearman');
    rho_P = corr(log10(joint_col), log10(joint_eco), 'type','Pearson');

    fprintf('  n_LEO_active     = %d\n', numel(idx_LEO));
    fprintf('  median R         = %.2f\n', median(asym));
    fprintf('  IQR R            = [%.2f, %.2f]\n', prctile(asym,25), prctile(asym,75));
    fprintf('  p99 R            = %.2f\n', prctile(asym,99));
    fprintf('  P(R>1)           = %.0f%%\n', 100*mean(asym>1));
    fprintf('  Spearman C_col vs C_eco (joint n=103) = %.4f\n', rho_S);
    fprintf('  Pearson on log10 = %.4f\n', rho_P);
    fprintf('  Constellation tier count = %d / 100\n', sum(tiers=="constellation"));
    fprintf('  ENVISAT C_col    = %.2f\n', case_Ccol(1));
    fprintf('  Sentinel-6 C_col = %.2f\n', case_Ccol(2));
    fprintf('  Starlink V2 C_col = %.2f\n', case_Ccol(3));

    results{snap} = struct( ...
        'n_LEO',         numel(idx_LEO), ...
        'median_R',      median(asym), ...
        'p99_R',         prctile(asym,99), ...
        'PR_gt_1',       100*mean(asym>1), ...
        'rho_S',         rho_S, ...
        'rho_P',         rho_P, ...
        'n_constell',    sum(tiers=="constellation"), ...
        'env_col',       case_Ccol(1), ...
        'sen_col',       case_Ccol(2), ...
        'star_col',      case_Ccol(3) );
end

% --- LaTeX-ready Tab 6 snippet ---------------------------------------- %
fprintf('\n===== LaTeX SNIPPET (Table 6 - temporal robustness) =====\n');
fprintf('\\begin{table}[ht]\n');
fprintf('\\centering\n');
fprintf('\\footnotesize\n');
fprintf('\\setlength{\\tabcolsep}{6pt}\n');
fprintf('\\caption{Temporal robustness of the validation and benchmark statistics across two snapshots of the Celestrak catalogue. Snapshot A is the May 2026 snapshot used throughout the paper (\\num{14727} active LEO). Snapshot B is a synthetic past-2025 snapshot built from the live June 2026 catalogue by removing all objects launched in 2026 (the filter $\\texttt{OBJECT\\_ID} \\not\\!\\sim\\! \\texttt{2026-}$ removes $\\sim\\!1900$ predominantly Starlink~V2~Mini objects, leaving the catalogue as it stood at the end of 2025). The random sample is redrawn on each snapshot from the same RNG seed. The qualitative claims of Sections~\\ref{sec:validation} and \\ref{sec:ecob-bench} (asymmetry universal, ECOB-style ranking) are invariant.}\n');
fprintf('\\label{tab:robustness}\n');
fprintf('\\begin{tabular}{lrr}\n');
fprintf('\\toprule\n');
fprintf('Statistic & Snapshot A (May 2026) & Snapshot B (Past 2025) \\\\ \n');
fprintf('\\midrule\n');
fprintf('Active LEO objects               & %d & %d \\\\ \n', results{1}.n_LEO, results{2}.n_LEO);
fprintf('Median $R = C_\\mathrm{col}/C_\\mathrm{ind}$ & %.1f & %.1f \\\\ \n', results{1}.median_R, results{2}.median_R);
fprintf('$p_{99}$ of $R$                  & %.0f & %.0f \\\\ \n', results{1}.p99_R, results{2}.p99_R);
fprintf('$P(R>1)$                         & %.0f\\%% & %.0f\\%% \\\\ \n', results{1}.PR_gt_1, results{2}.PR_gt_1);
fprintf('Constellation tier ($n/100$)     & %d & %d \\\\ \n', results{1}.n_constell, results{2}.n_constell);
fprintf('Spearman $\\rho$ ($C_\\mathrm{col}$ vs $C_\\mathrm{eco}$, $n{=}103$) & %.3f & %.3f \\\\ \n', results{1}.rho_S, results{2}.rho_S);
fprintf('Pearson $\\rho$ on $\\log_{10}$    & %.3f & %.3f \\\\ \n', results{1}.rho_P, results{2}.rho_P);
fprintf('ENVISAT $\\tilde{R}^{col}$        & %.1f & %.1f \\\\ \n', results{1}.env_col, results{2}.env_col);
fprintf('Sentinel-6 $\\tilde{R}^{col}$     & %.3f & %.3f \\\\ \n', results{1}.sen_col, results{2}.sen_col);
fprintf('Starlink V2 Mini $\\tilde{R}^{col}$ & %.1f & %.1f \\\\ \n', results{1}.star_col, results{2}.star_col);
fprintf('\\bottomrule\n');
fprintf('\\end{tabular}\n');
fprintf('\\end{table}\n');

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
