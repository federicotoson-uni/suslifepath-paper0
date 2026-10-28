%% paper0_catalogue.m  Revision (ASR): every active LEO object in place of the N = 100 sample
%
% Evaluates the individual and collective indices for all active LEO objects of the
% May 2026 Celestrak snapshot (R1-4, R2-6), with D = 100 independent draws of the
% non-orbital priors of paper0_validation.m. The formulas are those of
% individual_probability_flux.m and collective_probability.m, vectorised: only the
% shell counts depend on the catalogue, the priors enter multiplicatively. Each object
% is removed from its own shell count (R2-7.1). On the same draws it also computes:
%   - the snapshot-ECOB proxy of ecob_proxy.m and its Spearman correlation with R_col;
%   - A2: collective severity with the fragments of the struck target (R2-3);
%   - A4: dwell-averaged densities for eccentric orbits (R2-7.3);
%   - A5: collision avoidance with efficiency eta (R3-2, R3-16);
%   - the past-2025 snapshot of Section 4.2.
% Outputs: console statistics, figures/Figure_6_validation.pdf,
% figures/Figure_7_ecob_benchmark.pdf.

clear; close all;
THIS_DIR = fileparts(mfilename('fullpath'));
DATA_DIR = fullfile(THIS_DIR, '..', 'data');
FIG_DIR  = fullfile(THIS_DIR, '..', 'figures');
addpath(THIS_DIR);
rng(20260605, 'twister');
D = 100;                                            % draws of the priors

% Reference mission and the three cases, as in paper0_casestudies.m
ref_m  = struct('r',7078,'lt',7,'dt',25,'cost',1e8);
ref_sc = struct('ex_surf',5,'tot_surf',15);
casi = struct('nome', {'ENVISAT','Sentinel-6','Starlink V2 Mini'}, 'r', {7148, 7714, 6903}, ...
              'lt', {10, 5.5, 5}, 'dt', {200, 25, 5}, 'cost', {2.3e9, 4.0e8, 1.0e6}, ...
              'ex', {100, 10, 12}, 'tot', {250, 25, 30}, 'massa', {8211, 1192, 800}, ...
              'in_catalogue', {false, true, true});

cat = carica(fullfile(DATA_DIR,'celestrak_active.csv'), fullfile(DATA_DIR,'celestrak_debris.csv'));
N = numel(cat.rp);
pri = priors(cat.h, cat.costellazione, D);
out = indici(cat, pri, ref_m, ref_sc);
fprintf('\nActive LEO objects: %d (constellations by name %d, others below 600 km %d, 600-1500 km %d, above 1500 km %d)\n', N, ...
        sum(pri.gruppo == 1), sum(pri.gruppo == 2), sum(pri.gruppo == 3), sum(pri.gruppo == 4));
fprintf('Alone in their shell (zero indices, excluded): %d\n', sum(~out.ok));

% ---- A3: statistics over the whole catalogue ---------------------------------------------
nomi_t = {'all', 'constellations (by name)', 'others, below 600 km', 'others, 600-1500 km', 'others, above 1500 km'};
fprintf('\n%-26s %6s | %-22s | %-22s | %-22s\n', 'group', 'n', 'median R [min max]', ...
        'share R>1 [min max]', 'share R>100 [min max]');
for t = 0:4
    sel = out.ok & (t == 0 | pri.gruppo == t);
    if ~any(sel), fprintf('%-26s %6d | empty\n', nomi_t{t+1}, 0); continue, end
    medR = median(out.R(sel,:), 1);
    f1   = mean(out.R(sel,:) > 1, 1);
    f100 = mean(out.R(sel,:) > 100, 1);
    fprintf('%-26s %6d | %6.1f [%6.1f %6.1f] | %.3f [%.3f %.3f] | %.3f [%.3f %.3f]\n', nomi_t{t+1}, ...
            sum(sel), median(medR), min(medR), max(medR), median(f1), min(f1), max(f1), ...
            median(f100), min(f100), max(f100));
end
Rmed = median(out.R, 2);                            % per object, median over the draws
q = prctile(Rmed(out.ok), [25 50 75 90 99]);
fprintf('Per-object median: p25 %.1f, p50 %.1f, p75 %.1f, p90 %.1f, p99 %.1f\n', q);
fprintf('Objects with R > 1 in every draw: %d of %d\n', sum(all(out.R(out.ok,:) > 1, 2)), sum(out.ok));
[~, io] = sort(Rmed(out.ok), 'descend'); iok = find(out.ok);
fprintf('Upper tail (top 1%%): altitudes %.0f-%.0f km\n', min(cat.h(iok(io(1:ceil(0.01*numel(io)))))), ...
        max(cat.h(iok(io(1:ceil(0.01*numel(io)))))));

% The three cases with the same formulas (check against the modules)
for c = 1:3
    k = casi(c);
    m = struct('r',k.r,'lt',k.lt,'dt',k.dt,'cost',k.cost,'in_catalogue',k.in_catalogue);
    sc = struct('ex_surf',k.ex,'tot_surf',k.tot);
    casi(c).Rind = individual_probability_flux(cat.pop, m, sc) * k.cost / out.Rref_ind;
    casi(c).Rcol = collective_probability(cat.pop, m, sc, ref_sc.tot_surf) / out.Rref_col;
    fprintf('%-17s ind %9.4f  col %9.4f  R %8.2f\n', k.nome, casi(c).Rind, casi(c).Rcol, casi(c).Rcol / casi(c).Rind);
end

% How much of the Starlink load falls on its own fleet, and how much R depends on the reference altitude
nst = sum(abs(cat.rp(startsWith(cat.nome, "STARLINK")) - 6903) <= 25) - 1;      % without the element itself
nop = conta(cat.r_op, 6903) - 1;
fprintf('525 km shell: %d operational, %d Starlink (%.0f%%); against the other operators R %.0f, collective %.1f\n', ...
        nop, nst, 100 * nst / nop, casi(3).Rcol / casi(3).Rind * (1 - nst / nop), casi(3).Rcol * (1 - nst / nop));
g = @(r) conta(cat.r_op, r) / conta(cat.r_all, r);   % share of operational objects in the reference shell
for href = [525 550 600 700 770]
    fprintf('Reference at %d km: operational %3.0f%% of the shell, Starlink R %5.0f, catalogue median %6.1f\n', href, ...
            100 * g(6378 + href), casi(3).Rcol / casi(3).Rind * g(7078) / g(6378 + href), median(median(out.R(out.ok,:),1)) * g(7078) / g(6378 + href));
end

% ---- ECOB proxy (Section 5.1): Spearman with R_col, catalogue plus the three cases -------------
rhoS = zeros(1, D);
for d = 1:D
    x = [out.Rcol(out.ok,d); [casi.Rcol]'];
    e = [out.Reco(out.ok,d); arrayfun(@(k) eco_caso(cat, k, out), casi)'];
    rhoS(d) = corr(x, e, 'Type', 'Spearman');
end
fprintf('\nSpearman R_col vs ECOB proxy (n = %d): median %.3f [%.3f %.3f]\n', sum(out.ok) + 3, median(rhoS), min(rhoS), max(rhoS));

% ---- A2: fragments of the target (target of reference size) -----------------
f_m = sqrt(pri.tot / ref_sc.tot_surf);
Rcol2 = out.Rcol .* (f_m + 1) ./ (2 * f_m);
rA2 = arrayfun(@(d) corr(out.Rcol(out.ok,d), Rcol2(out.ok,d), 'Type','Spearman'), 1:D);
R2 = Rcol2 ./ out.Rind;
fprintf('\nA2: Spearman with R_col %.4f [%.4f %.4f]; median R %.1f (was %.1f); share R>1 %.3f\n', ...
        median(rA2), min(rA2), max(rA2), median(median(R2(out.ok,:),1)), median(median(out.R(out.ok,:),1)), median(mean(R2(out.ok,:) > 1, 1)));
for c = 1:3
    fc = sqrt(casi(c).tot / ref_sc.tot_surf);
    fprintf('    %-17s col %9.4f (factor %.2f)\n', casi(c).nome, casi(c).Rcol * (fc + 1) / (2 * fc), (fc + 1) / (2 * fc));
end

% ---- A4: eccentric orbits, orbit-averaged density against perigee density -----
ecc = find(cat.ra - cat.rp > 100 & out.ok);
Mg = linspace(0, 2*pi, 721); Mg(end) = [];
rap_op = zeros(numel(ecc),1); rap_all = rap_op;
for i = 1:numel(ecc)
    k = ecc(i); a = (cat.rp(k) + cat.ra(k)) / 2; e = (cat.ra(k) - cat.rp(k)) / (2 * a);
    r = a * (1 - e * cos(kepler(Mg, e)));
    proprio = abs(r - cat.rp(k)) <= cat.DH;           % the object itself counts only near its perigee
    rap_op(i)  = mean((conta(cat.r_op,  r) - proprio) ./ cat.V(r)) / out.rho_op(k);
    rap_all(i) = mean((conta(cat.r_all, r) - proprio) ./ cat.V(r)) / out.rho_all(k);
end
rapR = rap_op ./ rap_all;
fprintf('\nA4: %d objects with apogee - perigee > 100 km\n', numel(ecc));
fprintf('    operational density, orbit-averaged/perigee: median %.2f [%.2f %.2f]; effect on R: median %.2f [%.2f %.2f]\n', ...
        median(rap_op), min(rap_op), max(rap_op), median(rapR), min(rapR), max(rapR));

% ---- A5: collision avoidance with efficiency eta -----------------------------------
fman_nome = out.fman;                                 % manoeuvrable targets by name
for eta = [0.9 0.99]
    for v = 1:2
        if v == 1, fm = fman_nome; fm_ref = out.fman_ref; lab = 'manoeuvrable targets by name';
        else,      fm = ones(N,1); fm_ref = 1;            lab = 'all targets manoeuvre'; end
        g_m   = ((1-eta) * pri.lt + (1 - eta * fm) .* pri.dt) ./ ((1-eta) * (pri.lt + pri.dt));
        g_ref = ((1-eta) * ref_m.lt + (1 - eta * fm_ref) * ref_m.dt) / ((1-eta) * (ref_m.lt + ref_m.dt));
        R5 = out.R .* g_m / g_ref;
        fs = @(k) ((1-eta) * k.lt + (1 - eta * out.fman_caso(k.r)) * k.dt) / ((1-eta) * (k.lt + k.dt));
        if v == 2, fs = @(k) 1; end
        Rst = casi(3).Rcol / casi(3).Rind * fs(casi(3)) / g_ref;
        sel = out.ok & pri.gruppo == 1;
        fprintf('A5 eta %.2f, %-30s: Starlink R %7.1f; constellation median R %7.1f, share R>100 %.3f; all, share R>1 %.3f\n', ...
                eta, lab, Rst, median(median(R5(sel,:),1)), median(mean(R5(sel,:) > 100, 1)), median(mean(R5(out.ok,:) > 1, 1)));
    end
end
fprintf('Share of manoeuvrable targets by name: 525 km %.2f, 700 km %.2f, 1336 km %.2f\n', ...
        out.fman_caso(6903), out.fman_ref, out.fman_caso(7714));

% ---- Past snapshot (without the 2026 launches) ------------------------------------------
catB = carica(fullfile(DATA_DIR,'celestrak_active_past2025.csv'), fullfile(DATA_DIR,'celestrak_debris.csv'));
rng(20260605, 'twister');
priB = priors(catB.h, catB.costellazione, D);
outB = indici(catB, priB, ref_m, ref_sc);
fprintf('\nSnapshot without the 2026 launches: %d objects; median R %.1f [%.1f %.1f]; share R>1 %.3f\n', numel(catB.rp), ...
        median(median(outB.R(outB.ok,:),1)), min(median(outB.R(outB.ok,:),1)), max(median(outB.R(outB.ok,:),1)), median(mean(outB.R(outB.ok,:) > 1, 1)));

% ---- Figure 6: whole catalogue ---------------------------------------------------------
col_t = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];
R1 = out.R(:,1);                                     % one draw, like the sample of the submitted version
f6 = figure('Visible','off','Units','centimeters','Position',[2 2 18 7]);
subplot(1,2,1);
edges = -2:0.1:4;
histogram(log10(R1(out.ok & pri.gruppo == 1)), edges, 'FaceColor', col_t(1,:), 'EdgeColor','none'); hold on;
histogram(log10(R1(out.ok & pri.gruppo > 1)), edges, 'FaceColor', col_t(2,:), 'EdgeColor','none');
set(gca,'YScale','log');
xline(0, ':', 'Color', [0.4 0.4 0.4]);
xline(log10(casi(2).Rcol / casi(2).Rind), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
xline(log10(casi(3).Rcol / casi(3).Rind), '-.', 'Color', [0 0 0], 'LineWidth', 1);
ylim([0.8 1e6]);
xlabel('$\log_{10} R$','Interpreter','latex'); ylabel('Number of objects');
h = findobj(gca, 'Type', 'Histogram'); h = flipud(h); x = findobj(gca, 'Type', 'ConstantLine'); x = flipud(x);
legend([h(1:2); x(2:3)], {'Broadband constellations','Other objects','Sentinel-6','Starlink V2 Mini'}, 'Box','off','Location','northwest','FontSize',8);
title('(a)'); set(gca, 'TitleHorizontalAlignment', 'left');
subplot(1,2,2);
for t = 1:2
    s = out.ok & ((pri.gruppo == 1) == (t == 1));
    scatter(cat.h(s), R1(s), 4, col_t(t,:), 'filled', 'MarkerFaceAlpha', 0.35); hold on;
end
set(gca,'YScale','log'); yline(1, ':', 'Color', [0.4 0.4 0.4]); yline(100, ':', 'Color', [0.4 0.4 0.4]);
xlabel('Perigee altitude (km)'); ylabel('$R$','Interpreter','latex');
title('(b)'); set(gca, 'TitleHorizontalAlignment', 'left');
exportgraphics(f6, fullfile(FIG_DIR,'Figure_6_validation.pdf'), 'ContentType','vector');

% ---- Figure 7: ECOB proxy against R_col (first draw) ------------------------------
f7 = figure('Visible','off','Units','centimeters','Position',[2 2 11 9]);
for t = 1:2
    s = out.ok & ((pri.gruppo == 1) == (t == 1));
    scatter(out.Rcol(s,1), out.Reco(s,1), 4, col_t(t,:), 'filled', 'MarkerFaceAlpha', 0.35); hold on;
end
ec = arrayfun(@(k) eco_caso(cat, k, out), casi);
scatter([casi.Rcol], ec, 40, 'k', 'd', 'filled');
text([casi.Rcol] * 1.4, ec, {casi.nome}, 'FontSize', 8);
set(gca,'XScale','log','YScale','log'); lim = [1e-4 1e4]; plot(lim, lim, ':k'); xlim(lim); ylim(lim);
xlabel('$\tilde{R}^{col}$','Interpreter','latex'); ylabel('$\tilde{R}^{eco}$','Interpreter','latex');
legend({'Broadband constellations','Other objects'}, 'Box','off','Location','northwest','FontSize',8);
exportgraphics(f7, fullfile(FIG_DIR,'Figure_7_ecob_benchmark.pdf'), 'ContentType','vector');
fprintf('\nFigures 6 and 7 written to %s\n', FIG_DIR);

% =======================================================================================
function cat = carica(csv_att, csv_deb)
    cat.pop = databasecreator_real(csv_att, Inf, csv_deb);
    T  = readtable(csv_att);
    n  = T.MEAN_MOTION * 2*pi / 86400;
    a  = (398600.4418 ./ n.^2).^(1/3);
    rp = a .* (1 - T.ECCENTRICITY);
    al = rp - 6378;
    leo = isfinite(al) & al > 0 & al < 2000;
    assert(max(abs(rp(leo) - cat.pop.LEO(:,1))) < 1e-6, 'catalogue order differs from databasecreator_real');
    cat.rp = rp(leo); cat.ra = a(leo) .* (1 + T.ECCENTRICITY(leo)); cat.h = cat.rp - 6378;
    cat.nome = upper(string(T.OBJECT_NAME(leo)));
    cat.r_op  = sort(cat.pop.LEO(:,1));
    cat.r_all = sort([cat.pop.LEO(:,1); cat.pop.GEO(:,1); cat.pop.OTH(:,1); cat.pop.DEB(:,1)]);
    cat.DH = 25;
    cat.V  = @(r) 4/3 * pi * ((r + 25).^3 - (r - 25).^3);
    man = startsWith(cat.nome, ["STARLINK", "ONEWEB", "KUIPER"]);
    cat.costellazione = startsWith(cat.nome, ["STARLINK", "ONEWEB", "KUIPER", "HULIANWANG", "QIANFAN"]);
    cat.r_man = sort(cat.rp(man));
end

function n = conta(rs, r)
    n = arrayfun(@(x) sum(rs >= x - 25 & rs <= x + 25), r);
end

function pri = priors(h, cost5, D)
    % Same priors as tier_params (paper0_validation.m) and the mass prior of paper0_ecob_benchmark.m.
    % Revision: the five large broadband constellations take the constellation prior at every altitude.
    N = numel(h); col = @(v) reshape(v, [], 1);
    tier = ones(N,1); tier(h >= 600) = 2; tier(h >= 1500) = 3; tier(cost5) = 1;
    pri.gruppo = tier + 1; pri.gruppo(cost5) = 1;        % 1 constellations by name, 2-4 others by altitude
    cm = [5e6 150e6 400e6];  cs = [0.5 0.7 0.5];
    em = [30 50 80];  es = [10 20 30];  ef = [3 5 8];  tf = [2.0 2.2 2.5];
    LT = [5 7 10];  DT = [5 25 25];
    mm = [500 1500 3000];  ms = [0.4 0.6 0.5];
    pri.tier = tier;
    pri.cost = lognrnd(col(log(cm(tier))) .* ones(1,D), col(cs(tier)) .* ones(1,D));
    ex = max(col(ef(tier)), col(em(tier)) + col(es(tier)) .* randn(N, D));
    pri.ex = ex;  pri.tot = col(tf(tier)) .* ex;
    pri.lt = col(LT(tier)) .* ones(1,D);  pri.dt = col(DT(tier)) .* ones(1,D);
    pri.massa = lognrnd(col(log(mm(tier))) .* ones(1,D), col(ms(tier)) .* ones(1,D));
end

function out = indici(cat, pri, ref_m, ref_sc)
    YR = 3.15576e7; VREL = 10;
    n_op  = conta(cat.r_op,  cat.rp) - 1;             % without the object itself
    n_all = conta(cat.r_all, cat.rp) - 1;
    out.rho_op = n_op ./ cat.V(cat.rp);  out.rho_all = n_all ./ cat.V(cat.rp);
    out.Rref_ind = individual_probability_flux(cat.pop, ref_m, ref_sc) * ref_m.cost;
    out.Rref_col = collective_probability(cat.pop, ref_m, ref_sc, ref_sc.tot_surf);
    rho_ref = conta(cat.r_op, ref_m.r) / cat.V(ref_m.r);
    out.Rref_eco = rho_ref * ref_sc.ex_surf * 1e-6 * VREL * 200 * YR;
    sig = pri.ex * 1e-6;
    out.Rind = out.rho_all .* sig * VREL .* pri.lt * YR .* pri.cost / out.Rref_ind;
    out.Rcol = out.rho_op .* sig * VREL .* (pri.lt + pri.dt) * YR .* sqrt(pri.tot / ref_sc.tot_surf) / out.Rref_col;
    out.Reco = out.rho_op .* sig * VREL * 200 * YR .* (pri.massa / 500).^0.75 / out.Rref_eco;
    out.R = out.Rcol ./ out.Rind;
    out.ok = n_all > 0 & n_op > 0;
    nm = conta(cat.r_man, cat.rp) - startsWith(cat.nome, ["STARLINK", "ONEWEB", "KUIPER"]);
    out.fman = nm ./ max(n_op, 1);
    out.fman_caso = @(r) conta(cat.r_man, r) / max(conta(cat.r_op, r), 1);
    out.fman_ref = out.fman_caso(ref_m.r);
end

function e = eco_caso(cat, k, out)
    YR = 3.15576e7;
    n = conta(cat.r_op, k.r) - k.in_catalogue;
    e = n / cat.V(k.r) * k.ex * 1e-6 * 10 * 200 * YR * (k.massa / 500)^0.75 / out.Rref_eco;
end

function E = kepler(M, e)
    E = M;
    for it = 1:30
        E = E - (E - e .* sin(E) - M) ./ (1 - e .* cos(E));
    end
end
