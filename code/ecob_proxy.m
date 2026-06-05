%%   SIMPLIFIED ALGORITHM FOR COLLISION RISK ASSESSMENT - Paper 0 module   %%
function [P_eco, details] = ecob_proxy(PopulationData, mission, sc, M_ref, M_m)
% ----------------------------------------------------------------------- %
% Author: Federico Toson                                      (new module) %
% Snapshot-ECOB proxy for internal benchmarking against P_col              %
% ----------------------------------------------------------------------- %
%{
  Snapshot-only proxy of the ECOB index of Letizia, Colombo, Lewis and Krag
  (Letizia 2016, Letizia 2017, Letizia 2023), implemented as a SECOND
  collective-risk indicator on the same Celestrak snapshot used by
  collective_probability.m. It is NOT a re-implementation of the full ECOB,
  which requires long-term (200-year) debris-evolution modelling via
  MASTER / MOCAT-MC or equivalent. It is a transparent, reproducible
  approximation whose only purpose is INTERNAL BENCHMARKING of the
  P_col ranking produced by this paper against the documented ECOB
  literature ranking.

  THREE DIFFERENCES FROM collective_probability.m (the structural choices
  that align this proxy with ECOB):
    (1) Mass-based fragmentation weight following the NASA Standard Breakup
        Model fragment-count scaling adopted by ECOB:
            f_frag_eco = (M_m / M_ref)^0.75
        whereas collective_probability.m uses A^0.5 (area-based, conservative
        relative to ECOB).
    (2) Fixed long-term horizon T_eco = 200 yr replacing the operational +
        residual horizon T_op + T_res used in our index. This emulates the
        ECOB convention without simulating the actual debris evolution.
    (3) Operational-population density rho_op is the SAME as in
        collective_probability.m and is held FIXED at the snapshot value
        (no forward propagation). This is the explicit simplification of
        this proxy: long-term evolution would shift rho_op upward through
        debris generation and downward through atmospheric decay; the net
        effect for screening on a 200-yr horizon is documented in Letizia
        (2017) and is preserved only QUALITATIVELY here.

  Honest framing for the paper:
    * 'ECOB-proxy ranking matches our P_col ranking' is a defendable claim.
    * 'Our P_col equals ECOB up to a constant' is NOT defendable.
    * 'The ECOB-proxy is computable from a catalogue snapshot in seconds,
       and reproduces the published ENVISAT-very-high / sparse-mission-low
       qualitative ranking' is the design intent.

  INPUT:
    PopulationData - struct with LEO/GEO/OTH/DEB (Nx8); col 1 = perigee radius [km]
    mission        - struct: r [km perigee]
    sc             - struct: ex_surf [m^2]
    M_ref          - reference mission mass  [kg]
    M_m            - mission mass            [kg]

  OUTPUT:
    P_eco   - snapshot-ECOB-proxy index (expected ECOB-weighted events
              over 200 yr, normalised by mass^0.75 ratio)
    details - struct with intermediate quantities
%}
% ----------------------------------------------------------------------- %

% --- Constants --------------------------------------------------------- %
rearth = 6378;               % [km]
v_rel  = 10;                 % [km/s] mean relative velocity in LEO
dh     = 25;                 % [km] shell half-width (same as P_col)
sec_per_year = 3.15576e7;    % [s/y]
T_eco  = 200 * sec_per_year; % [s] ECOB long-term horizon

% --- Mission altitude (near-circular: use perigee radius) -------------- %
r_m = mission.r;             % [km] perigee radius
h_m = r_m - rearth;          % [km] altitude

% --- Operational density (SAME as collective_probability.m) ------------ %
LEO = PopulationData.LEO;
r_pop = LEO(:,1);                                  % [km] perigee radii
in_shell = (r_pop >= r_m - dh) & (r_pop <= r_m + dh);
n_op = sum(in_shell);

r_lo = rearth + h_m - dh;
r_hi = rearth + h_m + dh;
V_shell = (4/3)*pi*(r_hi^3 - r_lo^3);              % [km^3]
rho_op = n_op / V_shell;                           % [obj/km^3]

% --- Mission collision cross-section (SAME convention as P_col) -------- %
sigma_m = sc.ex_surf * 1e-6;                       % [km^2]

% --- ECOB-aligned mass-based fragmentation weight ---------------------- %
% NASA Standard Breakup Model: cumulative fragment count scales with target
% mass as M^0.75 (Johnson 2001). ECOB inherits this scaling.
f_frag_eco = (M_m / M_ref)^0.75;                   % [-]

% --- Snapshot ECOB-proxy ----------------------------------------------- %
rate  = rho_op * sigma_m * v_rel;                  % [events/s]
P_eco = rate * T_eco * f_frag_eco;                 % [exp. ECOB-weighted events]

% --- Diagnostics ------------------------------------------------------- %
details = struct( ...
    'altitude_km',   h_m, ...
    'n_op_in_shell', n_op, ...
    'V_shell_km3',   V_shell, ...
    'rho_op',        rho_op, ...
    'sigma_m_km2',   sigma_m, ...
    'v_rel_kms',     v_rel, ...
    'T_eco_s',       T_eco, ...
    'M_m_kg',        M_m, ...
    'M_ref_kg',      M_ref, ...
    'f_frag_eco',    f_frag_eco, ...
    'rate_per_s',    rate, ...
    'P_eco',         P_eco );
end
