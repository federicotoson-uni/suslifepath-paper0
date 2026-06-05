%%   SIMPLIFIED ALGORITHM FOR COLLISION RISK ASSESSMENT - Paper 0          %%
%%   Sensitivity check: f_frag exponent                                    %%
% ----------------------------------------------------------------------- %
% Author: Federico Toson                                                  %
% Verifies how the collective indices of the three case studies vary      %
% across alpha in {0.3, 0.4, 0.5, 0.6, 0.75}, and that the discrete       %
% five-class classification is invariant in this range.                   %
% Note: the parametric thresholds C*(h) of Table 3 are by construction    %
% INDEPENDENT of alpha, because A_m = A_ref -> f_frag = 1 there.          %
% ----------------------------------------------------------------------- %
clear; close all; clc;

% Surface values used in paper (m^2)
A_ref = 15;
A_env = 250;
A_sen = 25;
A_sta = 30;

% Baseline (alpha = 0.5) collective indices from paper0_casestudies.m run
Rcol_env_05 = 903.186;
Rcol_sen_05 = 0.047;
Rcol_sta_05 = 29.945;

% Sweep exponents
alphas = [0.30, 0.40, 0.50, 0.60, 0.75];

fprintf('\n=== Sensitivity of R_col_tilde to f_frag exponent alpha ===\n');
fprintf('alpha | ENVISAT (R_col / class) | Sentinel-6 (R_col / class) | Starlink V2 (R_col / class)\n');
fprintf('------+-------------------------+----------------------------+-------------------------------\n');

env_vals = zeros(1,numel(alphas));
sen_vals = zeros(1,numel(alphas));
sta_vals = zeros(1,numel(alphas));

for k = 1:numel(alphas)
    a = alphas(k);
    % Scaling factor: ratio (f_frag at new alpha) / (f_frag at 0.5)
    scale_env = (A_env/A_ref)^(a-0.5);
    scale_sen = (A_sen/A_ref)^(a-0.5);
    scale_sta = (A_sta/A_ref)^(a-0.5);
    env_vals(k) = Rcol_env_05 * scale_env;
    sen_vals(k) = Rcol_sen_05 * scale_sen;
    sta_vals(k) = Rcol_sta_05 * scale_sta;
    fprintf(' %0.2f | %8.1f  / %3s         | %7.4f  / %3s            | %7.2f  / %3s\n', ...
        a, env_vals(k), classOf(env_vals(k)), ...
           sen_vals(k), classOf(sen_vals(k)), ...
           sta_vals(k), classOf(sta_vals(k)));
end

fprintf('\nRanges for paper text:\n');
fprintf('  ENVISAT    R_col in [%.0f, %.0f]  (always VH)\n',  min(env_vals), max(env_vals));
fprintf('  Sentinel-6 R_col in [%.3f, %.3f]  (always L)\n',   min(sen_vals), max(sen_vals));
fprintf('  Starlink   R_col in [%.1f, %.1f]  (always VH)\n',  min(sta_vals), max(sta_vals));

% ====================================================================== %
%  FIGURE 5 - alpha-sensitivity (R_col vs alpha, with class thresholds)
% ====================================================================== %
SCRIPT_DIR = fileparts(mfilename('fullpath'));
REPO_ROOT  = fileparts(SCRIPT_DIR);
FIG_DIR    = fullfile(REPO_ROOT, 'figures');
if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

% Dense sweep for smooth lines
alphas_dense = linspace(0.30, 0.75, 60);
env_d = Rcol_env_05 .* (A_env/A_ref) .^ (alphas_dense - 0.5);
sen_d = Rcol_sen_05 .* (A_sen/A_ref) .^ (alphas_dense - 0.5);
sta_d = Rcol_sta_05 .* (A_sta/A_ref) .^ (alphas_dense - 0.5);

col_blue = [0.00 0.45 0.74];
col_oran = [0.85 0.33 0.10];
col_grn  = [0.20 0.55 0.25];
col_grey = [0.45 0.45 0.45];

f5 = figure('Visible','off','Units','centimeters','Position',[2 2 14 9]);
semilogy(alphas_dense, env_d, '-',  'Color',col_oran,'LineWidth',1.6, ...
         'DisplayName','ENVISAT'); hold on;
semilogy(alphas_dense, sta_d, '--', 'Color',col_blue,'LineWidth',1.6, ...
         'DisplayName','Starlink V2 Mini');
semilogy(alphas_dense, sen_d, '-.', 'Color',col_grn, 'LineWidth',1.6, ...
         'DisplayName','Sentinel-6');

% Markers at the discrete alphas used in Limitations text
% (HandleVisibility off so they do not appear in the legend)
for k = 1:numel(alphas)
    plot(alphas(k), env_vals(k), 'o', 'MarkerFaceColor',col_oran,'MarkerEdgeColor','none','MarkerSize',5,'HandleVisibility','off');
    plot(alphas(k), sen_vals(k), 'o', 'MarkerFaceColor',col_grn, 'MarkerEdgeColor','none','MarkerSize',5,'HandleVisibility','off');
    plot(alphas(k), sta_vals(k), 'o', 'MarkerFaceColor',col_blue,'MarkerEdgeColor','none','MarkerSize',5,'HandleVisibility','off');
end

% class threshold lines
thr  = [1e-2, 1e-1, 1e0, 1e1];
clab = {'L','M','H','VH'};
for t = 1:numel(thr)
    yline(thr(t),'LineStyle',':','Color',col_grey,'LineWidth',0.6,'HandleVisibility','off');
    text(0.76, thr(t), clab{t}, 'FontSize',8,'Color',col_grey, ...
         'VerticalAlignment','middle');
end

xlim([0.30 0.78]);
ylim([1e-2 1e4]);
xlabel('Fragmentation-weight exponent $\alpha$','Interpreter','latex');
ylabel('Normalised collective index $\tilde{R}^{col}$','Interpreter','latex');
legend('show','Location','east','FontSize',8,'Box','off');
grid on; set(gca,'FontSize',10,'Box','on','GridAlpha',0.15,'Layer','top');

exportgraphics(f5, fullfile(FIG_DIR,'fig5_alpha_sensitivity.pdf'),'ContentType','vector');
fprintf('\nFigure 5 (alpha sensitivity) written to %s\n', ...
        fullfile(FIG_DIR,'fig5_alpha_sensitivity.pdf'));

function c = classOf(R)
    if     R < 1e-2, c = 'VL';
    elseif R < 1e-1, c = 'L';
    elseif R < 1e0,  c = 'M';
    elseif R < 1e1,  c = 'H';
    else,            c = 'VH';
    end
end
