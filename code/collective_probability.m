%%   SIMPLIFIED ALGORITHM FOR COLLISION RISK ASSESSMENT - Paper 0 module   %%
function [P_col, details] = collective_probability(PopulationData, mission, sc, S_ref)
% ----------------------------------------------------------------------- %
% Author: Federico Toson                                                  %
% Collective collision-risk contribution                      (new module)%
% ----------------------------------------------------------------------- %
%{
  Computes the COLLECTIVE contribution of mission m: how much the mission
  (during its operational life + residual orbital lifetime after disposal)
  adds to the collision risk experienced by the OPERATIONAL population at
  its altitude band.

  POSITIONING vs ECOB (read this before citing the module):
  This is a SIMPLIFIED, REAL-TIME ENGINEERING PROXY for early-design
  screening. It does NOT replace and does NOT compete with the full ECOB
  index [Letizia, Colombo, Lewis & Krag], which is the gold standard and
  requires long-term (200-year) debris-evolution modelling. Our module
  gives a first-order estimate computable in seconds on a laptop, using
  only the catalogued population snapshot. For the multi-domain SSCI
  framework (Paper 1) the ORBITAL input should be the full ECOB index,
  not this proxy. This proxy belongs to Paper 0 (rapid classification);
  ECOB belongs to Paper 1 (rigorous orbital component). Always cite ECOB.

  MODEL  (kinetic theory of the debris flux, Kessler & Cour-Palais 1978,
  the NASA/ESA standard that also underpins ECOB; cf. Paper 0 Sec. 3.2):

      collision rate  =  rho_op * sigma_m * v_rel          [events/s]
      P_col_base      =  rate * T_total                     [exp. events]
      P_col           =  P_col_base * f_frag                [ECOB-weighted]

  where
      rho_op  - spatial number density of operational objects in a shell
                of half-width dh around the mission altitude   [obj/km^3]
      sigma_m - mission collision cross-section                [km^2]
      v_rel   - mean relative velocity in LEO                  [km/s]
      T_total - operational + residual orbital lifetime        [s]
      f_frag  - fragmentation-severity weight (large/massive objects
                generate larger fragment clouds; ECOB-inspired)  [-]

  Dimensional check:
      [obj/km^3]*[km^2]*[km/s] = [obj/s]  ->  *[s] = [obj] (dimensionless)

  INPUT:
    PopulationData - resident population (from databasecreator.m); fields
                     LEO/GEO/OTH/DEB, col 1 = perigee radius [km],
                     col 7 = cross-section/surface [m^2]
    mission        - struct: r [km perigee], lt [y], dt [y residual], ...
    sc             - struct: ex_surf [m^2], tot_surf [m^2], rho [kg/m^3]

  OUTPUT:
    P_col   - collective contribution (expected fragment-weighted events)
    details - struct with all intermediate quantities for inspection
%}
% ----------------------------------------------------------------------- %

% --- Constants --------------------------------------------------------- %
rearth = 6378;               % [km]
v_rel  = 10;                 % [km/s] mean relative velocity in LEO (typical;
                             % consistent with Rossi and Farinella 1992,
                             % 9.65 +/- 0.88 km/s)
dh     = 25;                 % [km] shell half-width around mission altitude
sec_per_year = 3.15576e7;    % [s/y]
if nargin < 4 || isempty(S_ref)
    S_ref = 25;              % [m^2] default reference cross-section (typical small sat)
end

% --- Mission altitude (near-circular: use perigee radius) -------------- %
r_m  = mission.r;            % [km] perigee radius
h_m  = r_m - rearth;         % [km] altitude

% --- Operational population in the altitude shell ---------------------- %
% Operational objects = LEO catalogue (exclude DEB = debris, which are the
% hazard *source*, not the population to be protected). Adjust if needed.
LEO = PopulationData.LEO;
r_pop = LEO(:,1);                                  % [km] perigee radii
in_shell = (r_pop >= r_m - dh) & (r_pop <= r_m + dh);
n_op = sum(in_shell);

% Shell volume between (rearth+h-dh) and (rearth+h+dh) ------------------- %
r_lo = rearth + h_m - dh;
r_hi = rearth + h_m + dh;
V_shell = (4/3)*pi*(r_hi^3 - r_lo^3);              % [km^3]

rho_op = n_op / V_shell;                           % [obj/km^3]

% --- Mission collision cross-section ----------------------------------- %
% Convert exposed surface [m^2] -> [km^2]; ex_surf is the proxy used in the
% thesis severity model, kept here for consistency.
sigma_m = sc.ex_surf * 1e-6;                       % [km^2]

% --- Total exposure time: operational + residual orbital lifetime ------ %
% mission.dt is the post-mission residual/decommissioning time [y].
T_total = (mission.lt + mission.dt) * sec_per_year; % [s]

% --- Base collective rate * time (expected events) --------------------- %
rate       = rho_op * sigma_m * v_rel;             % [events/s]
P_col_base = rate * T_total;                       % [expected events]

% --- Fragmentation-severity weight (ECOB-inspired) --------------------- %
% Larger objects produce larger fragment clouds. We use the total surface as
% a size proxy, normalised to the reference cross-section S_ref (passed in:
% the reference mission's tot_surf) so the reference mission has f_frag = 1.
% A 0.5 exponent is used (more conservative than the 0.75 fragment-count
% scaling of NASA breakup models) to avoid over-amplifying very large objects
% such as defunct platforms. The area-based proxy is the design choice
% of this paper; a mass-based alternative (M^0.75, NASA Standard Breakup
% Model) is implemented separately in ecob_proxy.m and benchmarked in
% Section 5.1 of the manuscript.
f_frag = (sc.tot_surf / S_ref)^0.5;                % [-]

% --- Collective contribution ------------------------------------------- %
P_col = P_col_base * f_frag;

% --- Diagnostics ------------------------------------------------------- %
details = struct( ...
    'altitude_km',   h_m, ...
    'n_op_in_shell', n_op, ...
    'V_shell_km3',   V_shell, ...
    'rho_op',        rho_op, ...
    'sigma_m_km2',   sigma_m, ...
    'v_rel_kms',     v_rel, ...
    'T_total_s',     T_total, ...
    'rate_per_s',    rate, ...
    'P_col_base',    P_col_base, ...
    'f_frag',        f_frag, ...
    'P_col',         P_col );
end
