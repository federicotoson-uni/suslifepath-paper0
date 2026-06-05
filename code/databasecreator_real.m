%%   SIMPLIFIED ALGORITHM FOR COLLISION RISK ASSESSMENT - Paper 0         %%
function [PopulationData] = databasecreator_real(csvpath, max_per_class, debris_csv)
% ----------------------------------------------------------------------- %
% Author: Federico Toson                                      (new module) %
% Real-population database from Celestrak GP CSV (active + debris)         %
% ----------------------------------------------------------------------- %
%{
  Drop-in replacement for databasecreator.m that builds PopulationData from
  REAL catalogues (Celestrak GP, FORMAT=csv) instead of a synthetic random
  population. Fix #1: gives the true altitude distribution, so the model
  sees the real congestion (the ~9,600 active satellites in the 450-550 km
  Starlink band, the debris clusters at ~770-800 km, the near-empty 1336 km
  Sentinel-6 band, etc.).

  Output structure (compatible with the flux modules), each row Nx8:
    col 1  perigee radius [km]        col 5  life duration [y]
    col 2  argument of perigee [rad]  col 6  RAAN [rad]
    col 3  eccentricity               col 7  exposed surface [m^2]
    col 4  inclination [rad]          col 8  cost [USD/g]

  Orbital elements from the TLE/GP record:
    a  = (mu / n^2)^(1/3),  n = MEAN_MOTION * 2*pi / 86400   [rad/s]
    rp = a * (1 - ecc)                                       [km]
  Fields absent from a TLE use documented defaults (active: lt=5, surf=5,
  cost=10; debris: lt=50, surf=1, cost=0). TODO: refine surface from SATCAT
  RCS size if needed.

  INPUT:
    csvpath        - Celestrak GP CSV of ACTIVE satellites (celestrak_active.csv)
    max_per_class  - (optional) cap per class for fast tests; default Inf (all).
                     Random sample -> preserves the altitude-distribution shape.
    debris_csv     - (optional) Celestrak GP CSV of DEBRIS (celestrak_debris.csv).
                     If empty/missing, DEB is left empty.

  OUTPUT:
    PopulationData - struct with fields LEO, GEO, OTH, DEB (Nx8 each).
%}
% ----------------------------------------------------------------------- %
if nargin < 2 || isempty(max_per_class), max_per_class = Inf; end
if nargin < 3, debris_csv = ''; end

rearth = 6378;          % [km]

% --- Active satellites -> LEO / GEO / OTH ------------------------------ %
[Ma, alt_a, incdeg_a] = parse_gp(csvpath, 5, 5, 10);   % lt, surf, cost defaults
valid = isfinite(alt_a) & (alt_a > 0);
isLEO = valid & (alt_a < 2000);
isGEO = valid & (abs(alt_a - 35786) < 500) & (incdeg_a < 10);
isOTH = valid & ~isLEO & ~isGEO;

LEO = Ma(isLEO,:);
GEO = Ma(isGEO,:);
OTH = Ma(isOTH,:);

% --- Debris -> DEB ----------------------------------------------------- %
if ~isempty(debris_csv) && isfile(debris_csv)
    [Md, alt_d, ~] = parse_gp(debris_csv, 50, 1, 0);   % debris: lt=50, surf=1, cost=0
    DEB = Md(isfinite(alt_d) & (alt_d > 0), :);
else
    DEB = zeros(0,8);
end

% --- Optional random cap for fast tests -------------------------------- %
if isfinite(max_per_class)
    LEO = subsample(LEO, max_per_class);
    GEO = subsample(GEO, max_per_class);
    OTH = subsample(OTH, max_per_class);
    DEB = subsample(DEB, max_per_class);
end

PopulationData.LEO = LEO;
PopulationData.GEO = GEO;
PopulationData.OTH = OTH;
PopulationData.DEB = DEB;

fprintf(['databasecreator_real: LEO=%d  GEO=%d  OTH=%d  DEB=%d\n'], ...
        size(LEO,1), size(GEO,1), size(OTH,1), size(DEB,1));
end

% ----------------------------------------------------------------------- %
function [M, alt, incdeg] = parse_gp(csvpath, lt_def, surf_def, cost_def)
% Parse a Celestrak GP CSV into the Nx8 matrix + perigee altitude [km].
    mu   = 398600.4418;   % [km^3/s^2]
    conv = pi/180;
    T = readtable(csvpath);

    n_rev  = T.MEAN_MOTION;                 % [rev/day]
    ecc    = T.ECCENTRICITY;
    incdeg = T.INCLINATION;                 % [deg]
    inc    = incdeg * conv;                 % [rad]
    w      = T.ARG_OF_PERICENTER * conv;    % [rad]
    ra     = T.RA_OF_ASC_NODE    * conv;    % [rad]

    n  = n_rev * 2*pi / 86400;              % [rad/s]
    a  = (mu ./ n.^2).^(1/3);               % [km]
    rp = a .* (1 - ecc);                    % [km] perigee radius
    alt = rp - 6378;                        % [km] perigee altitude

    N = height(T);
    M = [rp, w, ecc, inc, lt_def*ones(N,1), ra, surf_def*ones(N,1), cost_def*ones(N,1)];
end

% ----------------------------------------------------------------------- %
function Y = subsample(X, k)
    if size(X,1) <= k
        Y = X;
    else
        idx = randperm(size(X,1), k);
        Y = X(idx,:);
    end
end
