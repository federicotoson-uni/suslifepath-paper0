%%   SIMPLIFIED ALGORITHM FOR COLLISION RISK ASSESSMENT - Paper 0 module   %%
function [P_ind, details] = individual_probability_flux(PopulationData, mission, sc)
% ----------------------------------------------------------------------- %
% Author: Federico Toson                                      (new module) %
% Individual collision risk via kinetic flux                              %
% ----------------------------------------------------------------------- %
%{
  Individual collision risk of mission m: the expected number of impacts
  ON m during its operational life. Computed with the SAME kinetic-flux
  approach as collective_probability.m, so the two scenarios are
  methodologically consistent and both robust in sparse altitude bands
  (unlike the discrete geometric intersection of the thesis probability.m,
  which under-counts conjunctions where the population is thin).

  MODEL (kinetic theory of the debris flux; cf. Paper 0 Section 3.2):

      P_ind = rho_local * sigma_m * v_rel * T_op            [expected impacts]

  where
      rho_local - spatial number density of ALL catalogued objects in the
                  mission's altitude shell (anything can hit m)  [obj/km^3]
      sigma_m   - mission collision cross-section                [km^2]
      v_rel     - mean relative velocity in LEO                  [km/s]
      T_op      - operational lifetime = exposure window         [s]

  Difference from collective_probability.m (the asymmetry that the paper
  highlights):
    INDIVIDUAL  -> "is m hit?"      : density of ALL objects, T_op only,
                                      NO fragmentation weight.
    COLLECTIVE  -> "does m harm others?" : density of OPERATIONAL objects,
                                      T_op + residual lifetime, x f_frag.

  Note vs ECOB: as for the collective module, this is a fast engineering
  proxy for screening, NOT a replacement for ECOB. See Paper 0 Section 3.2.

  INPUT:
    PopulationData - struct with LEO/GEO/OTH/DEB (Nx8), col 1 = perigee radius [km]
    mission        - struct: r [km perigee], lt [y operational], ...
    sc             - struct: ex_surf [m^2], ...

  OUTPUT:
    P_ind   - expected impacts on m during operations (dimensionless)
    details - struct with intermediate quantities for inspection
%}
% ----------------------------------------------------------------------- %

% --- Constants --------------------------------------------------------- %
rearth = 6378;               % [km]
v_rel  = 10;                 % [km/s] mean relative velocity in LEO
dh     = 25;                 % [km] shell half-width
sec_per_year = 3.15576e7;    % [s/y]

% --- Mission altitude (near-circular: use perigee radius) -------------- %
r_m = mission.r;             % [km] perigee radius
h_m = r_m - rearth;          % [km] altitude

% --- Local density: ALL catalogued objects can hit m ------------------- %
ALL = [PopulationData.LEO; PopulationData.GEO; PopulationData.OTH; PopulationData.DEB];
r_pop = ALL(:,1);                                  % [km] perigee radii
in_shell = (r_pop >= r_m - dh) & (r_pop <= r_m + dh);
n_local = sum(in_shell);

r_lo = rearth + h_m - dh;
r_hi = rearth + h_m + dh;
V_shell = (4/3)*pi*(r_hi^3 - r_lo^3);              % [km^3]
rho_local = n_local / V_shell;                     % [obj/km^3]

% --- Mission collision cross-section ----------------------------------- %
sigma_m = sc.ex_surf * 1e-6;                       % [km^2]

% --- Exposure time = operational lifetime ------------------------------ %
T_op = mission.lt * sec_per_year;                  % [s]

% --- Expected impacts on m --------------------------------------------- %
rate  = rho_local * sigma_m * v_rel;               % [impacts/s]
P_ind = rate * T_op;                               % [expected impacts]

% --- Diagnostics ------------------------------------------------------- %
details = struct( ...
    'altitude_km',     h_m, ...
    'n_local_in_shell',n_local, ...
    'V_shell_km3',     V_shell, ...
    'rho_local',       rho_local, ...
    'sigma_m_km2',     sigma_m, ...
    'v_rel_kms',       v_rel, ...
    'T_op_s',          T_op, ...
    'rate_per_s',      rate, ...
    'P_ind',           P_ind );
end
