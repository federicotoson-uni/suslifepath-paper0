%%   SIMPLIFIED ALGORITHM FOR COLLISION RISK ASSESSMENT - Paper 0         %%
%%   Case studies driver: ENVISAT / Sentinel-6 / Starlink V2 Mini         %%
% ----------------------------------------------------------------------- %
% Author: Federico Toson                                                  %
% Driver script for Section 4 of Paper 0                      (new script) %
% ----------------------------------------------------------------------- %
%{
  Uses the flux modules:
    databasecreator_real.m, individual_probability_flux.m,
    collective_probability.m, risk_index.m
  to produce Table 2 of the paper (individual + collective risk index and
  classification for three representative LEO objects).

  Paths are auto-discovered relative to this script's location.

  The population is built from the real 2026 Celestrak GP catalogue via
  databasecreator_real.m (active + tracked debris). The collective index
  uses the kinetic-flux formulation of collective_probability.m
  (operational density rho_op x sigma_m x v_rel x T_total x f_frag); see
  Section 3.2 of the manuscript for the formal definition of P_col.
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

% Load the resident population (real Celestrak GP catalogue snapshot).
CSV    = fullfile(DATA_DIR, 'celestrak_active.csv');
DEBCSV = fullfile(DATA_DIR, 'celestrak_debris.csv');
PopulationData = databasecreator_real(CSV, Inf, DEBCSV);   % full real pop + debris
% Fallback to synthetic random population (thesis 2022):
% PopulationData = databasecreator("y","n","n","n","n");

% ====================================================================== %
%  CASE DEFINITIONS  (perigee radius r = altitude + 6378 km)
%  Cross-sections (ex_surf/tot_surf) are engineering estimates from public
%  fact sheets (ESA for ENVISAT and Sentinel-6; FCC SAT-MOD-20230207-00022
%  for Starlink V2 Mini). Replacement costs are first-order class
%  estimates (build cost for ENVISAT; programme cost for Sentinel-6; unit
%  marginal cost for Starlink V2 Mini); the IDA 2023 LEO satellite cost
%  report (D-33436) is used as the order-of-magnitude reference for
%  mid-size operational small satellites.
% ====================================================================== %
cases = struct();

% --- 1) ENVISAT : large defunct platform, ~770 km polar -------------- %
cases(1).name        = 'ENVISAT';
cases(1).mission.yol = 2002;  cases(1).mission.lt = 10;  cases(1).mission.dt = 200;
cases(1).mission.it  = 0.2;   cases(1).mission.r  = 7148;            % 770+6378
cases(1).mission.ecc = 0.001; cases(1).mission.w  = 0;   cases(1).mission.ra = 0;
cases(1).mission.inc = 98.4;  cases(1).mission.pl = 'G'; cases(1).mission.net = 0;
cases(1).mission.man = 'U';   cases(1).mission.cost = 2.3e9;         % original build cost
cases(1).sc.ex_surf  = 100;   cases(1).sc.tot_surf = 250; cases(1).sc.rho = 750; % large RCS

% --- 2) Sentinel-6 Michael Freilich : operational, ~1336 km ---------- %
cases(2).name        = 'Sentinel-6';
cases(2).mission.yol = 2020;  cases(2).mission.lt = 5.5; cases(2).mission.dt = 25;
cases(2).mission.it  = 0.2;   cases(2).mission.r  = 7714;            % 1336+6378
cases(2).mission.ecc = 0.0001;cases(2).mission.w  = 0;   cases(2).mission.ra = 0;
cases(2).mission.inc = 66;    cases(2).mission.pl = 'G'; cases(2).mission.net = 0;
cases(2).mission.man = 'U';   cases(2).mission.cost = 4.0e8;
cases(2).sc.ex_surf  = 10;    cases(2).sc.tot_surf = 25;  cases(2).sc.rho = 480;

% --- 3) Starlink V2 Mini : constellation element, ~525 km ------------ %
cases(3).name        = 'Starlink V2 Mini';
cases(3).mission.yol = 2023;  cases(3).mission.lt = 5;   cases(3).mission.dt = 5;
cases(3).mission.it  = 0.2;   cases(3).mission.r  = 6903;            % 525+6378
cases(3).mission.ecc = 0.0001;cases(3).mission.w  = 0;   cases(3).mission.ra = 0;
cases(3).mission.inc = 53;    cases(3).mission.pl = 'C'; cases(3).mission.net = 1;
cases(3).mission.man = 'U';   cases(3).mission.cost = 1.0e6;
cases(3).sc.ex_surf  = 12;    cases(3).sc.tot_surf = 30;  cases(3).sc.rho = 400;

% ====================================================================== %
%  REFERENCE MISSION (fixed, external) - stable normaliser
%  A "typical operational small satellite": ~500 kg class, 700 km SSO,
%  compliant 25-yr disposal. Evaluated in BOTH scenarios. Using a fixed
%  external reference (not the median of the 3 cases) makes the normalisation
%  stable and SSCI interpretable as "how many times a typical reference".
% ====================================================================== %
ref_m  = struct('yol',2020,'lt',7,'dt',25,'it',0.2,'r',7078, ...   % 700+6378
    'ecc',0.001,'w',0,'ra',0,'inc',98,'pl','G','net',0,'man','U','cost',1e8);
ref_m.my = ref_m.yol + (1:ref_m.lt);
ref_sc = struct('ex_surf',5,'tot_surf',15,'rho',500);
S_ref  = ref_sc.tot_surf;   % fragmentation reference -> reference has f_frag = 1

% Flux-based individual + collective (consistent, robust in sparse bands).
% Severity differs by scenario:
%   INDIVIDUAL severity = mission replacement cost (value of m lost if hit).
%   COLLECTIVE severity = homogeneous: the harm to OTHERS depends on how many
%     operational satellites m threatens (already in P_col), NOT on m's own
%     value. So R_col = P_col (the constant per-target value cancels in the
%     normalisation). This is why a cheap constellation element can still be a
%     high collective risk.
[ref_Pind, ~] = individual_probability_flux(PopulationData, ref_m, ref_sc);
[ref_Pcol, ~] = collective_probability(PopulationData, ref_m, ref_sc, S_ref);
R_ref_ind = ref_Pind * ref_m.cost;   if R_ref_ind<=0, R_ref_ind = eps; end
R_ref_col = ref_Pcol;                if R_ref_col<=0, R_ref_col = eps; end

% ====================================================================== %
%  RUN
% ====================================================================== %
N = numel(cases);
R_raw_ind = zeros(1,N);
R_raw_col = zeros(1,N);

for i = 1:N
    m  = cases(i).mission;
    sc = cases(i).sc;
    m.my = m.yol + (1:max(1,round(m.lt)));

    % --- Individual risk: P_ind (kinetic flux) x cost ----------------- %
    [P_ind, inddet]    = individual_probability_flux(PopulationData, m, sc);
    R_raw_ind(i)       = P_ind * m.cost;

    % --- Collective risk: P_col (kinetic flux); homogeneous severity --- %
    [P_col, coldet]    = collective_probability(PopulationData, m, sc, S_ref);
    R_raw_col(i)       = P_col;

    cases(i).inddetails = inddet;   % inspect: rho_local, sigma_m, T_op, ...
    cases(i).coldetails = coldet;   % inspect: rho_op, sigma_m, f_frag, ...
end

% ====================================================================== %
%  TABLE 1
% ====================================================================== %
fprintf('\n================ Paper 0 - Table 1 ================\n');
fprintf('R_ref individual = %.3e | R_ref collective = %.3e\n\n', ...
        R_ref_ind, R_ref_col);
fprintf('%-18s | %-22s | %-22s\n','Case','Individual','Collective');
fprintf('%s\n', repmat('-',1,66));
for i = 1:N
    % risk_index expects (probtot,sev,R_ref); we already have R_raw, so we
    % pass R_raw as probtot and sev=1 to reuse the normaliser cleanly.
    [ci, cls_i, lbl_i] = risk_index(R_raw_ind(i), 1, R_ref_ind, 'individual'); %#ok<ASGLU>
    [cc, cls_c, lbl_c] = risk_index(R_raw_col(i), 1, R_ref_col, 'collective'); %#ok<ASGLU>
    fprintf('%-18s | %5.2f  class %d (%-9s) | %5.2f  class %d (%-9s)\n', ...
        cases(i).name, ci, cls_i, lbl_i, cc, cls_c, lbl_c);
end
fprintf('===================================================\n');

%{
  EXPECTED QUALITATIVE PATTERN (to validate, not to hardcode):
   - ENVISAT     : individual High, collective Very High (large RCS, debris)
   - Sentinel-6  : individual/collective Low (protected 1336 km band)
   - Starlink V2 : individual Medium, collective High (dense 525 km band)
  If the run contradicts this, check: database currency, cross-sections,
  collective_factor placeholder, and the conjunction tolerance (toll=25).
%}
