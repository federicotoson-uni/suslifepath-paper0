%%   SIMPLIFIED ALGORITHM FOR COLLISION RISK ASSESSMENT - Paper 0         %%
%%   Educational smallsat corollary: CubeSat / PocketCube at 525 km       %%
% ----------------------------------------------------------------------- %
% Author: Federico Toson                                                  %
% Driver script for the closing paragraph of Section 4.3      (new script) %
% ----------------------------------------------------------------------- %
%{
  Quantifies the educational-smallsat corollary of eq. (8): sigma_m
  cancels in the ratio R = R_col/R_ind, but f_frag enters the collective
  index only, so a kilogram-scale spacecraft does not inherit the
  constellation asymmetry unchanged. The smaller fragmentation footprint
  attenuates R; the lower replacement cost amplifies it.

  Cases (Starlink shell, 525 km):
    - Starlink V2 Mini  (control: must reproduce R ~ 317 of Table 2)
    - 3U CubeSat, 1.0 MUSD   (university-class upper budget)
    - 2P PocketCube, 50 kUSD (typical educational budget)

  Also reports the geometry-specific cost threshold C* (asymmetry >= 100)
  for each educational geometry, the analogue of Table 3 with K(h)
  rescaled by the case f_frag and lifetimes.

  Numbers quoted in the paper: R ~ 43 (3U at 1 MUSD), C* ~ 0.43 MUSD,
  R ~ 460 (2P at 50 kUSD).
%}
% ----------------------------------------------------------------------- %
clear; close all; clc;

% --- Auto-discover paths relative to this script's location ----------- %
SCRIPT_DIR = fileparts(mfilename('fullpath'));   % .../code
REPO_ROOT  = fileparts(SCRIPT_DIR);              % .../repo root
DATA_DIR   = fullfile(REPO_ROOT, 'data');
addpath(SCRIPT_DIR);
assert(exist('risk_index','file')==2, ...
       'risk_index.m not found on path; check SCRIPT_DIR.');

CSV    = fullfile(DATA_DIR, 'celestrak_active.csv');
DEBCSV = fullfile(DATA_DIR, 'celestrak_debris.csv');
PopulationData = databasecreator_real(CSV, Inf, DEBCSV);

% --- Reference mission (identical to paper0_casestudies.m) ------------ %
ref_m  = struct('yol',2020,'lt',7,'dt',25,'it',0.2,'r',7078, ...
    'ecc',0.001,'w',0,'ra',0,'inc',98,'pl','G','net',0,'man','U','cost',1e8);
ref_m.my = ref_m.yol + (1:ref_m.lt);
ref_sc = struct('ex_surf',5,'tot_surf',15,'rho',500);
S_ref  = ref_sc.tot_surf;

[ref_Pind, ~] = individual_probability_flux(PopulationData, ref_m, ref_sc);
[ref_Pcol, ~] = collective_probability(PopulationData, ref_m, ref_sc, S_ref);
R_ref_ind = ref_Pind * ref_m.cost;
R_ref_col = ref_Pcol;

% --- Cases at the Starlink shell (525 km) ------------------------------ %
tc = struct();

% 1) Starlink V2 Mini (control, parameters of paper0_casestudies.m)
tc(1).name = 'Starlink V2 Mini (control)';
tc(1).m  = struct('yol',2023,'lt',5,'dt',5,'it',0.2,'r',6903,'ecc',0.0001, ...
    'w',0,'ra',0,'inc',53,'pl','C','net',1,'man','U','cost',1.0e6);
tc(1).sc = struct('ex_surf',12,'tot_surf',30,'rho',400);

% 2) 3U CubeSat, university-class: 10x10x30 cm, ~4 kg, SSO rideshare.
%    ex_surf = largest face 0.1x0.3 m; tot_surf = 2*(0.1*0.1)+4*(0.1*0.3).
tc(2).name = '3U CubeSat, 1.0 MUSD';
tc(2).m  = struct('yol',2026,'lt',2,'dt',6,'it',0.2,'r',6903,'ecc',0.0001, ...
    'w',0,'ra',0,'inc',97.5,'pl','G','net',0,'man','U','cost',1.0e6);
tc(2).sc = struct('ex_surf',0.03,'tot_surf',0.14,'rho',1300);

% 3) 2P PocketCube (RedPill-class): 5x5x10 cm, ~0.8 kg.
tc(3).name = '2P PocketCube, 50 kUSD';
tc(3).m  = struct('yol',2026,'lt',1,'dt',4,'it',0.2,'r',6903,'ecc',0.0001, ...
    'w',0,'ra',0,'inc',97.5,'pl','G','net',0,'man','U','cost',5.0e4);
tc(3).sc = struct('ex_surf',0.005,'tot_surf',0.025,'rho',1300);

% --- Run ---------------------------------------------------------------- %
fprintf('\n========== Educational smallsat corollary (Sec. 4.3) ==========\n');
fprintf('%-28s %12s %12s %10s %10s %9s %12s\n', ...
    'Case','R_ind~','R_col~','Cls_ind','Cls_col','R','C* (USD)');
for i = 1:numel(tc)
    m = tc(i).m; sc = tc(i).sc;
    m.my = m.yol + (1:max(1,round(m.lt)));
    [Pind,~] = individual_probability_flux(PopulationData, m, sc);
    [Pcol,~] = collective_probability(PopulationData, m, sc, S_ref);
    Rind = (Pind*m.cost)/R_ref_ind;
    Rcol = Pcol/R_ref_col;
    [~,~,li] = risk_index(Pind*m.cost,1,R_ref_ind,'individual');
    [~,~,lc] = risk_index(Pcol,1,R_ref_col,'collective');
    % Geometry-specific C*: R(C) = K_case/C with K_case = Rcol*R_ref_ind/Pind
    K_case = Rcol * R_ref_ind / Pind;
    fprintf('%-28s %12.4g %12.4g %10s %10s %9.1f %12.3g\n', ...
        tc(i).name, Rind, Rcol, li, lc, Rcol/Rind, K_case/100);
end
fprintf('\nf_frag: 3U = %.4f | 2P = %.4f | Starlink V2 Mini = %.4f\n', ...
    sqrt(0.14/S_ref), sqrt(0.025/S_ref), sqrt(30/S_ref));
fprintf(['\nsigma cancels in R; f_frag does not (collective only).\n' ...
    'Attenuation (f_frag) vs amplification (cost): at typical educational\n' ...
    'budgets the cost term dominates and R exceeds the constellation value.\n']);
