clc; clear; close all;
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'outputs');
if ~exist(outDir, 'dir'), mkdir(outDir); end

cfg.order = 3;
cfg.Nx = 513;
cfg.Ny = 513;
cfg.Ng = cfg.Nx * cfg.Ny;
cfg.N_list = unique(round(logspace(log10(8), log10(2000), 140)));
cfg.actual_file = fullfile(rootDir, 'outputfiles');

% Relative weights for the theoretical model. These are not wall-clock seconds.
cfg.w_coeff_direct = 1.0;
cfg.w_coeff_spectral = 0.25;
cfg.w_recon_direct = 1.0;
cfg.w_recon_spectral_terms = 8.0;
cfg.w_fft = 5.0;

Nc = cfg.N_list(:);
Ng = cfg.Ng;

% Term-count model
C2 = Nc .* (Nc - 1) / 2;
C3 = Nc .* (Nc - 1) .* (Nc - 2) / 6;
L2 = Nc .* (Nc - 1);

Nterms_direct = Nc + Nc + L2 + Nc + L2 + 4 * C3;
Nterms_spectral = Nc + Nc + L2 + Nc + L2 + C3;

Mcoeff_direct = (2/3) * Nc.^3 + 1 * Nc.^2 - (2/3) * Nc;
Mcoeff_spectral = (1/6) * Nc.^3 + 2 * Nc.^2 - (7/6) * Nc;
Mrecon_direct = Ng .* Nterms_direct;
Mrecon_spectral = cfg.w_recon_spectral_terms * Nterms_spectral + cfg.w_fft * Ng * log2(Ng);

Mtotal_direct = cfg.w_coeff_direct * Mcoeff_direct + cfg.w_recon_direct * Mrecon_direct;
Mtotal_spectral = cfg.w_coeff_spectral * Mcoeff_spectral + Mrecon_spectral;

assert(all(diff(Nc) > 0), 'N_list must be strictly increasing.');
assert(all(Mtotal_direct > 0) && all(diff(Mtotal_direct) >= 0), 'Direct complexity must be positive and monotone.');
assert(all(Mtotal_spectral > 0) && all(diff(Mtotal_spectral) >= 0), 'Spectral complexity must be positive and monotone.');

actual = read_actual_benchmark(cfg.actual_file);

plotData = struct();
plotData.cfg = cfg;
plotData.theory = table(Nc, Nterms_direct, Nterms_spectral, Mcoeff_direct, Mcoeff_spectral, ...
    Mrecon_direct, Mrecon_spectral, Mtotal_direct, Mtotal_spectral, ...
    Mtotal_direct ./ Mtotal_spectral, ...
    'VariableNames', {'Nc', 'Nterms_direct', 'Nterms_spectral', 'Mcoeff_direct', 'Mcoeff_spectral', ...
    'Mrecon_direct', 'Mrecon_spectral', 'Mtotal_direct', 'Mtotal_spectral', 'speedup_direct_over_spectral'});
plotData.actual = actual;

c_direct = [0.06, 0.35, 0.67];
c_spec = [0.86, 0.34, 0.12];
c_direct_soft = [0.48, 0.66, 0.86];
c_spec_soft = [0.96, 0.66, 0.48];
font_axis = 12.5;
font_label = 14.5;
font_legend = 11.5;
font_panel = 14;
lw_main = 2.8;
lw_aux = 1.7;
ms = 6.0;

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [1.5 1.5 24.0 10.8]);
tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'tight');

% -------------------- Left: theoretical --------------------
ax1 = nexttile;
ax1.Color = 'w';
h1 = semilogy(ax1, Nc, Mtotal_direct, '-', 'Color', c_direct, 'LineWidth', lw_main);
hold(ax1, 'on');
h2 = semilogy(ax1, Nc, Mtotal_spectral, '-', 'Color', c_spec, 'LineWidth', lw_main);
semilogy(ax1, Nc, Mcoeff_direct, '--', 'Color', c_direct_soft, 'LineWidth', lw_aux);
semilogy(ax1, Nc, Mrecon_direct, ':', 'Color', c_direct_soft, 'LineWidth', lw_aux);
semilogy(ax1, Nc, Mcoeff_spectral, '--', 'Color', c_spec_soft, 'LineWidth', lw_aux);
semilogy(ax1, Nc, Mrecon_spectral, ':', 'Color', c_spec_soft, 'LineWidth', lw_aux);
marker_idx = 1:8:numel(Nc);
semilogy(ax1, Nc(marker_idx), Mtotal_direct(marker_idx), 'o', ...
    'Color', c_direct, 'MarkerFaceColor', c_direct, 'MarkerSize', ms, 'LineStyle', 'none');
semilogy(ax1, Nc(marker_idx), Mtotal_spectral(marker_idx), 's', ...
    'Color', c_spec, 'MarkerFaceColor', c_spec, 'MarkerSize', ms, 'LineStyle', 'none');
style_axes(ax1, font_axis);
xlabel(ax1, 'Number of retained components N_c', 'FontName', 'Times New Roman', 'FontSize', font_label);
ylabel(ax1, 'Computational complexity', 'FontName', 'Times New Roman', 'FontSize', font_label);
hCoeff = plot(ax1, NaN, NaN, 'k--', 'LineWidth', lw_aux);
hRecon = plot(ax1, NaN, NaN, 'k:', 'LineWidth', lw_aux);
text(ax1, 0.03, 0.96, '(a)', 'Units', 'normalized', ...
    'FontName', 'Times New Roman', 'FontSize', font_panel, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'top');

% -------------------- Right: actual --------------------
ax2 = nexttile;
ax2.Color = 'w';
h3 = semilogy(ax2, actual.Nc, actual.direct_total_s, '-o', ...
    'Color', c_direct, 'LineWidth', lw_main, 'MarkerSize', ms, 'MarkerFaceColor', c_direct);
hold(ax2, 'on');
h4 = semilogy(ax2, actual.Nc, actual.spectral_total_s, '-s', ...
    'Color', c_spec, 'LineWidth', lw_main, 'MarkerSize', ms, 'MarkerFaceColor', c_spec);
semilogy(ax2, actual.Nc, actual.direct_coeff_s, '--', 'Color', c_direct_soft, 'LineWidth', lw_aux);
semilogy(ax2, actual.Nc, actual.direct_recon_s, ':', 'Color', c_direct_soft, 'LineWidth', lw_aux);
semilogy(ax2, actual.Nc, actual.spectral_coeff_s, '--', 'Color', c_spec_soft, 'LineWidth', lw_aux);
semilogy(ax2, actual.Nc, actual.spectral_recon_s, ':', 'Color', c_spec_soft, 'LineWidth', lw_aux);
style_axes(ax2, font_axis);
xlabel(ax2, 'Number of retained components N_c', 'FontName', 'Times New Roman', 'FontSize', font_label);
ylabel(ax2, 'Measured runtime (s)', 'FontName', 'Times New Roman', 'FontSize', font_label);
text(ax2, 0.03, 0.96, '(b)', 'Units', 'normalized', ...
    'FontName', 'Times New Roman', 'FontSize', font_panel, 'FontWeight', 'bold', ...
    'VerticalAlignment', 'top');

lgd = legend(ax1, [h1, h2, hCoeff, hRecon], ...
    {'Direct reconstruction', 'Spectral reconstruction', ...
    'Coefficient term', 'Reconstruction term'}, ...
    'Orientation', 'vertical', 'NumColumns', 1, ...
    'Location', 'southeast', 'Box', 'on', ...
    'Color', 'w', 'EdgeColor', [0.82 0.82 0.82], ...
    'FontName', 'Times New Roman', 'FontSize', font_legend);

out_png = fullfile(outDir, 'mf12_theory_vs_actual_time_scaling.png');
out_mat = fullfile(outDir, 'mf12_theory_vs_actual_time_scaling.mat');
out_csv = fullfile(outDir, 'mf12_theory_vs_actual_time_scaling_data.csv');

csvData = outerjoin(plotData.theory, plotData.actual, ...
    'Keys', 'Nc', ...
    'MergeKeys', true, ...
    'Type', 'full');
csvData = sortrows(csvData, 'Nc');

save(out_mat, 'plotData');
writetable(csvData, out_csv);
exportgraphics(fig, out_png, 'Resolution', 600);

fprintf('Saved figure: %s\n', out_png);
fprintf('Saved data  : %s\n', out_mat);
fprintf('Saved csv   : %s\n', out_csv);
fprintf('Trend check: theory and measurement should be qualitatively similar, not identical.\n');
fprintf('Reason: asymptotics capture growth rate, while real timings also include constants, MATLAB overhead, memory traffic, and implementation details.\n');

function style_axes(ax, font_axis)
grid(ax, 'on');
set(ax, ...
    'FontName', 'Times New Roman', ...
    'FontSize', font_axis, ...
    'LineWidth', 1.0, ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'GridColor', [0.82 0.82 0.82], ...
    'GridAlpha', 0.36, ...
    'TickDir', 'out', ...
    'YScale', 'log');
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
end

function T = read_actual_benchmark(filename)
if ~isfile(filename)
    error('Actual benchmark file not found: %s', filename);
end

txt = fileread(filename);
lines = splitlines(string(txt));
vals = [];
pat = ['N=\s*(\d+)\s*\|\s*direct total\s*([0-9eE+\-.]+)s\s*' ...
    '\(coeff\s*([0-9eE+\-.]+)s\s*\+\s*recon\s*([0-9eE+\-.]+)s\)\s*\|\s*' ...
    'spectral total\s*([0-9eE+\-.]+)s\s*' ...
    '\(coeff\s*([0-9eE+\-.]+)s\s*\+\s*recon\s*([0-9eE+\-.]+)s\)\s*\|\s*' ...
    'speedup\s*([0-9eE+\-.]+)x\s*\|\s*max\|deta\|\s*([0-9eE+\-.]+)'];
for i = 1:numel(lines)
    tok = regexp(strtrim(lines(i)), pat, 'tokens', 'once');
    if ~isempty(tok)
        vals(end+1, :) = str2double(tok); %#ok<AGROW>
    end
end

if isempty(vals)
    error('No benchmark rows could be parsed from: %s', filename);
end
T = table( ...
    vals(:,1), vals(:,2), vals(:,3), vals(:,4), ...
    vals(:,5), vals(:,6), vals(:,7), vals(:,8), vals(:,9), ...
    'VariableNames', {'Nc', 'direct_total_s', 'direct_coeff_s', 'direct_recon_s', ...
    'spectral_total_s', 'spectral_coeff_s', 'spectral_recon_s', ...
    'speedup_direct_over_spectral', 'eta_max_abs_err'});
end

