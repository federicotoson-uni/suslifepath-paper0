%% paper0_sentinel6_eol.m  Revision (ASR): Sentinel-6 collective index and end of life (R2-2)
%
% The case study charges the residual lifetime at the density of the mission shell
% (screening convention, T_res = 25 yr). Two physical end-of-life options instead:
%   graveyard: the satellite stays in the 1336 km band, as its predecessor Jason-2 did,
%              over the 200-year horizon used for ENVISAT and the ECOB proxy;
%   descent:   a 25-year compliant disposal, which needs a perigee inside the dense
%              shells; the residual phase is charged at the operational density
%              averaged over a fixed disposal ellipse (1336 km x h_p), weighted by
%              the time spent at each altitude. A fixed apogee understates the time
%              spent low as the apogee decays, so these values are lower bounds.
% The individual index depends on T_op only and does not change.

clear;
THIS_DIR = fileparts(mfilename('fullpath'));
addpath(THIS_DIR);
Pop = databasecreator_real(fullfile(THIS_DIR,'..','data','celestrak_active.csv'), Inf, ...
                           fullfile(THIS_DIR,'..','data','celestrak_debris.csv'));
RE = 6378; DH = 25;
V  = @(r) 4/3 * pi * ((r + DH).^3 - (r - DH).^3);
r_op = Pop.LEO(:,1);
r_s6 = 7714;
rho_op = @(r) arrayfun(@(x) sum(abs(r_op - x) <= DH) - (abs(x - r_s6) <= DH), r) ./ V(r);   % without Sentinel-6A

ref_m  = struct('r',7078,'lt',7,'dt',25,'cost',1e8);
ref_sc = struct('ex_surf',5,'tot_surf',15);
R_ref_col = collective_probability(Pop, ref_m, ref_sc, ref_sc.tot_surf);
s6 = struct('r',r_s6,'lt',5.5,'dt',25,'cost',4.0e8,'in_catalogue',true);
s6_sc = struct('ex_surf',10,'tot_surf',25);
[P_col, d] = collective_probability(Pop, s6, s6_sc, ref_sc.tot_surf);
per_anno = P_col / (s6.lt + s6.dt);                   % collective per year at the 1336 km density
classe = @(x) 1 + sum(x >= 10.^(-2:1));
etich  = {'VL','L','M','H','VH'};
stampa = @(nome, x) fprintf('%-50s %9.3f  %s\n', nome, x, etich{classe(x)});

fprintf('\nSentinel-6, normalised collective index (operational objects in the shell without 6A: %d)\n', d.n_op_in_shell);
stampa('convention: 25 yr at the 1336 km density', P_col / R_ref_col);
stampa('graveyard in the band, 200-yr horizon', per_anno * (s6.lt + 200) / R_ref_col);
M = linspace(0, 2*pi, 3601); M(end) = [];
for hp = [350 450 550]
    rp = RE + hp; ra = r_s6; a = (rp + ra) / 2; e = (ra - rp) / (2 * a);
    E = M; for it = 1:30, E = E - (E - e*sin(E) - M) ./ (1 - e*cos(E)); end
    fattore = mean(rho_op(a * (1 - e*cos(E)))) / rho_op(r_s6);
    x = per_anno * (s6.lt + s6.dt * fattore) / R_ref_col;
    stampa(sprintf('25-yr descent, perigee %d km (density x%.0f)', hp, fattore), x);
end

% Descent with the collision-avoidance model of the revised manuscript: only Starlink, OneWeb and Kuiper avoid a derelict
T  = readtable(fullfile(THIS_DIR,'..','data','celestrak_active.csv'));
n  = T.MEAN_MOTION * 2*pi / 86400;  a_ = (398600.4418 ./ n.^2).^(1/3);  rp_ = a_ .* (1 - T.ECCENTRICITY);
r_man = rp_(startsWith(upper(string(T.OBJECT_NAME)), ["STARLINK", "ONEWEB", "KUIPER"]) & rp_ - RE > 0 & rp_ - RE < 2000);
fman = @(r) arrayfun(@(x) sum(abs(r_man - x) <= DH) / max(sum(abs(r_op - x) <= DH), 1), r);
for eta = [0.9 0.99]
    for hp = [350 450 550]
        rp = RE + hp; ra = r_s6; a = (rp + ra) / 2; e = (ra - rp) / (2 * a);
        E = M; for it = 1:30, E = E - (E - e*sin(E) - M) ./ (1 - e*cos(E)); end
        r = a * (1 - e*cos(E));
        fattore = mean(rho_op(r) .* (1 - eta * fman(r))) / rho_op(r_s6);
        stampa(sprintf('descent with manoeuvres %.0f%%, perigee %d km (x%.0f)', 100*eta, hp, fattore), per_anno * (s6.lt + s6.dt * fattore) / R_ref_col);
    end
end

% The reference under the convention questioned in R2-2: its 25 yr along a disposal ellipse
rho_rif = @(r) arrayfun(@(x) sum(abs(r_op - x) <= DH), r) ./ V(r);
for hp = [550 600]
    rp = RE + hp; ra = ref_m.r; a = (rp + ra) / 2; e = (ra - rp) / (2 * a);
    E = M; for it = 1:30, E = E - (E - e*sin(E) - M) ./ (1 - e*cos(E)); end
    fattore = mean(rho_rif(a * (1 - e*cos(E)))) / rho_rif(ref_m.r);
    k = (ref_m.lt + ref_m.dt * fattore) / (ref_m.lt + ref_m.dt);
    fprintf('Reference disposed on 700 x %d km: density x%.1f, reference collective x%.1f, Starlink R %.0f\n', hp, fattore, k, 317.3 / k);
end
