%%   SIMPLIFIED ALGORITHM FOR COLLISION RISK ASSESSMENT - Paper 0 module   %%
function [R_cont, R_class, R_label] = risk_index(probtot, sev, R_ref, mode)
% ----------------------------------------------------------------------- %
% Author: Federico Toson                                                  %
% Risk Index normaliser & classifier                          (new module)%
% ----------------------------------------------------------------------- %
%{
  This module extends the thesis formulation RISK = probtot * sev
  (Toson, 2022) into a CONTINUOUS, normalised risk index suitable for
  classification AND for aggregation across missions, bridging the
  categorical ECSS scheme with continuous orbital-sustainability indices.

  It is the building-block cited by Paper 1 (SSCI) as the orbital risk
  component of the multi-domain Space Sustainability Composite Indicator.

  INPUT:
    probtot - collision probability over the integration horizon
              (from probability.m, individual or collective variant)
    sev     - severity (from severity.m); economic + environmental
    R_ref   - reference-mission risk used for normalisation
              (median of the catalogued LEO population in the band)
    mode    - 'individual' or 'collective' (affects only labelling/log)

  OUTPUT:
    R_cont  - continuous, dimensionless risk index (multiplier of R_ref)
    R_class - discrete class 1..5 (compatible with ECSS-M-ST-80C)
    R_label - text label of the class
%}
% ----------------------------------------------------------------------- %

if nargin < 4
    mode = 'individual';
end

% --- Raw risk: thesis formulation R = P x S ---------------------------- %
R_raw = probtot * sev;

% --- Normalisation against the reference mission ----------------------- %
if R_ref <= 0
    error('risk_index:badRef','R_ref must be positive.');
end
R_cont = R_raw / R_ref;

% --- Discrete classification (logarithmic, ECSS-compatible 5 classes) -- %
% Boundaries are one order of magnitude apart, centred on the reference:
%   class 1 (Very Low):  R_cont <  1e-2
%   class 2 (Low):       1e-2 <= R_cont < 1e-1
%   class 3 (Medium):    1e-1 <= R_cont < 1e0
%   class 4 (High):      1e0  <= R_cont < 1e1
%   class 5 (Very High): R_cont >= 1e1
edges   = [1e-2, 1e-1, 1e0, 1e1];
R_class = 1 + sum(R_cont >= edges);

labels  = {'Very Low','Low','Medium','High','Very High'};
R_label = labels{R_class};

% --- Console log ------------------------------------------------------- %
fprintf('[%s] R_cont = %.3f  ->  class %d (%s)\n', ...
        mode, R_cont, R_class, R_label);
end
