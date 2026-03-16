% BENCHMARK_DIRECT_VS_SPECTRAL_VSN
% Compare end-to-end runtime vs N:
%   Direct sum:      mf12_direct_coefficients + mf12_direct_surface
%   Spectral (new):  mf12_spectral_coefficients + mf12_spectral_surface

clc;
clear;
close all;

scriptDir = fileparts(mfilename('fullpath'));
matlabDir = fileparts(scriptDir);
repoDir = fileparts(matlabDir);
run(fullfile(matlabDir, 'setup_paths.m'));
outDir = fullfile(repoDir, 'outputs');
if ~exist(outDir, 'dir'), mkdir(outDir); end
setenv('MF12_PROGRESS', '0');

% -------------------- Config --------------------
cfg.g = 9.81;
cfg.h = 150;
cfg.Ux = 0;
cfg.Uy = 0;
cfg.order = 3;
cfg.t = 0.0;

cfg.Lx = 3000;
cfg.Ly = 3000;
cfg.Nx = 64;
cfg.Ny = 64;

cfg.N_list = [40 60 80];
cfg.repeats =1;
cfg.seed = 20260306;
cfg.plot_style_name = 'mf12_benchmark_v1';

% Directional wave-group-like setup
cfg.Tp = 12;
cfg.theta1 = deg2rad(25);
cfg.theta2 = deg2rad(-35);
cfg.sig_k_factor = 0.12;
cfg.sig_t1 = deg2rad(10);
cfg.sig_t2 = deg2rad(12);
cfg.w1 = 0.55;
cfg.w2 = 0.45;
cfg.Hs_target = 3.5;

rng(cfg.seed);

% -------------------- Build candidate k-space once --------------------
dkx = 2*pi/cfg.Lx;
dky = 2*pi/cfg.Ly;
kx_idx_all = (-floor(cfg.Nx/2)):(ceil(cfg.Nx/2)-1);
ky_idx_all = (-floor(cfg.Ny/2)):(ceil(cfg.Ny/2)-1);
[KXI, KYI] = meshgrid(kx_idx_all, ky_idx_all);
kx_all = KXI(:) * dkx;
ky_all = KYI(:) * dky;
kmag_all = hypot(kx_all, ky_all);
theta_all = atan2(ky_all, kx_all);

keep = (kx_all > 0) | (kx_all == 0 & ky_all > 0);
kx_all = kx_all(keep);
ky_all = ky_all(keep);
kmag_all = kmag_all(keep);
theta_all = theta_all(keep);

kp = (2*pi/cfg.Tp)^2 / cfg.g;
sig_k = cfg.sig_k_factor * kp;
Sk = exp(-0.5*((kmag_all - kp)/sig_k).^2);
St1 = exp(-0.5*(angle(exp(1i*(theta_all - cfg.theta1)))/cfg.sig_t1).^2);
St2 = exp(-0.5*(angle(exp(1i*(theta_all - cfg.theta2)))/cfg.sig_t2).^2);
W = Sk .* (cfg.w1*St1 + cfg.w2*St2);

valid = W > 1e-8 * max(W);
kx_all = kx_all(valid);
ky_all = ky_all(valid);
W = W(valid);
[~, idx_sort] = sort(W, 'descend');

if max(cfg.N_list) > numel(idx_sort)
    error('Max N=%d exceeds available energetic candidates=%d.', max(cfg.N_list), numel(idx_sort));
end

% Fixed grid for direct reconstruction
x = (0:cfg.Nx-1) * (cfg.Lx/cfg.Nx);
y = (0:cfg.Ny-1) * (cfg.Ly/cfg.Ny);
[X, Y] = meshgrid(x, y);

rows = [];
k = 0;

fprintf('=== Direct vs Spectral end-to-end runtime vs N ===\n');
fprintf('Grid=%dx%d, order=%d, repeats=%d\n\n', cfg.Nx, cfg.Ny, cfg.order, cfg.repeats);

for N = cfg.N_list
    idx = idx_sort(1:N);
    kx = kx_all(idx); kx = kx(:).';
    ky = ky_all(idx); ky = ky(:).';
    Wsel = W(idx); Wsel = Wsel(:).';

    amp_raw = Wsel / max(Wsel);
    m0_raw = 0.5 * sum(amp_raw.^2);
    Hs_raw = 4 * sqrt(max(m0_raw, eps));
    amp = amp_raw * (cfg.Hs_target / Hs_raw);

    xf = 0.5 * cfg.Lx;
    yf = 0.5 * cfg.Ly;
    omega_lin = sqrt(cfg.g*hypot(kx,ky).*tanh(cfg.h*hypot(kx,ky)));
    phase = -(kx*xf + ky*yf) + omega_lin*cfg.t;
    a = amp .* cos(phase); a = a(:).';
    b = amp .* sin(phase); b = b(:).';

    t_direct_coeff = zeros(cfg.repeats, 1);
    t_direct_recon = zeros(cfg.repeats, 1);
    t_direct_total = zeros(cfg.repeats, 1);
    t_spectral_coeff = zeros(cfg.repeats, 1);
    t_spectral_recon = zeros(cfg.repeats, 1);
    t_spectral_total = zeros(cfg.repeats, 1);
    eta_ref = [];
    eta_spec = [];

    for r = 1:cfg.repeats
        tic;
        coeffs_d = mf12_direct_coefficients(cfg.order, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy);
        t_direct_coeff(r) = toc;

        tic;
        [eta_d, ~] = mf12_direct_surface(cfg.order, coeffs_d, X, Y, cfg.t);
        t_direct_recon(r) = toc;
        t_direct_total(r) = t_direct_coeff(r) + t_direct_recon(r);

        tic;
        coeffs_s = mf12_spectral_coefficients(cfg.order, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy, struct('enable_subharmonic', false));
        t_spectral_coeff(r) = toc;
        coeffs_s.third_order_subharmonic_mode = 'skip';

        tic;
        [eta_s, ~] = mf12_spectral_surface(coeffs_s, cfg.Lx, cfg.Ly, cfg.Nx, cfg.Ny, cfg.t);
        t_spectral_recon(r) = toc;
        t_spectral_total(r) = t_spectral_coeff(r) + t_spectral_recon(r);

        if r == 1
            eta_ref = eta_d;
            eta_spec = eta_s;
        end
    end

    k = k + 1;
    rows(k).N = N; %#ok<SAGROW>
    rows(k).direct_coeff_mean_s = mean(t_direct_coeff);
    rows(k).direct_recon_mean_s = mean(t_direct_recon);
    rows(k).direct_total_mean_s = mean(t_direct_total);
    rows(k).spectral_coeff_mean_s = mean(t_spectral_coeff);
    rows(k).spectral_recon_mean_s = mean(t_spectral_recon);
    rows(k).spectral_total_mean_s = mean(t_spectral_total);
    rows(k).speedup_direct_over_spectral = rows(k).direct_total_mean_s / max(rows(k).spectral_total_mean_s, eps);
    rows(k).eta_max_abs_err = max(abs(eta_ref(:) - eta_spec(:)));

    fprintf('N=%3d | direct total %.3fs (coeff %.3fs + recon %.3fs) | spectral total %.3fs (coeff %.3fs + recon %.3fs) | speedup %.2fx | max|deta| %.3e\n', ...
        N, rows(k).direct_total_mean_s, rows(k).direct_coeff_mean_s, rows(k).direct_recon_mean_s, ...
        rows(k).spectral_total_mean_s, rows(k).spectral_coeff_mean_s, rows(k).spectral_recon_mean_s, ...
        rows(k).speedup_direct_over_spectral, rows(k).eta_max_abs_err);
end

T = struct2table(rows);
disp(T);

ts = datestr(now, 'yyyymmdd_HHMMSS');
csv_name = fullfile(outDir, ['benchmark_direct_vs_spectral_vsN_' ts '.csv']);
mat_name = fullfile(outDir, ['benchmark_direct_vs_spectral_vsN_' ts '.mat']);
png_name = fullfile(outDir, ['benchmark_direct_vs_spectral_vsN_' ts '.png']);
writetable(T, csv_name);
plotData = build_plot_data(T, cfg, ts, csv_name, mat_name, png_name);
save(mat_name, 'plotData');

fig = plot_benchmark_figure(plotData);
exportgraphics(fig, png_name, 'Resolution', 180);

fprintf('\nSaved: %s\n', csv_name);
fprintf('Saved: %s\n', mat_name);
fprintf('Saved: %s\n', png_name);

function plotData = build_plot_data(T, cfg, ts, csv_name, mat_name, png_name)
plotData = struct();
plotData.timestamp = ts;
plotData.cfg = cfg;
plotData.table = T;
plotData.files = struct('csv', csv_name, 'mat', mat_name, 'png', png_name);
plotData.series = struct( ...
    'N', T.N, ...
    'direct_total', T.direct_total_mean_s, ...
    'spectral_total', T.spectral_total_mean_s, ...
    'direct_coeff', T.direct_coeff_mean_s, ...
    'direct_recon', T.direct_recon_mean_s, ...
    'spectral_coeff', T.spectral_coeff_mean_s, ...
    'spectral_recon', T.spectral_recon_mean_s, ...
    'speedup', T.speedup_direct_over_spectral, ...
    'eta_err', T.eta_max_abs_err);
end

function fig = plot_benchmark_figure(plotData)
T = plotData.table;
style = benchmark_plot_style();

fig = figure('Color', style.figure_color, 'Position', [80 80 1280 520]);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile;
hold(ax1, 'on');
plot(ax1, T.N, T.direct_total_mean_s, '-o', ...
    'Color', style.direct_color, 'LineWidth', 2.2, ...
    'MarkerSize', 7, 'MarkerFaceColor', style.direct_color);
plot(ax1, T.N, T.spectral_total_mean_s, '-s', ...
    'Color', style.spectral_color, 'LineWidth', 2.2, ...
    'MarkerSize', 7, 'MarkerFaceColor', style.spectral_color);
hold(ax1, 'off');
format_axes(ax1, style);
set(ax1, 'YScale', 'log');
xlabel(ax1, 'N components');
ylabel(ax1, 'Mean total runtime (s)');
title(ax1, 'End-to-End Runtime');
legend(ax1, {'Direct sum', 'Spectral'}, 'Location', 'northwest', 'Box', 'off');

ax2 = nexttile;
hold(ax2, 'on');
plot(ax2, T.N, T.direct_coeff_mean_s, '--o', ...
    'Color', style.direct_light, 'LineWidth', 1.8, ...
    'MarkerSize', 6, 'MarkerFaceColor', style.direct_light);
plot(ax2, T.N, T.direct_recon_mean_s, '-o', ...
    'Color', style.direct_color, 'LineWidth', 2.0, ...
    'MarkerSize', 6, 'MarkerFaceColor', style.direct_color);
plot(ax2, T.N, T.spectral_coeff_mean_s, '--s', ...
    'Color', style.spectral_light, 'LineWidth', 1.8, ...
    'MarkerSize', 6, 'MarkerFaceColor', style.spectral_light);
plot(ax2, T.N, T.spectral_recon_mean_s, '-s', ...
    'Color', style.spectral_color, 'LineWidth', 2.0, ...
    'MarkerSize', 6, 'MarkerFaceColor', style.spectral_color);
hold(ax2, 'off');
format_axes(ax2, style);
set(ax2, 'YScale', 'log');
xlabel(ax2, 'N components');
ylabel(ax2, 'Mean runtime (s)');
title(ax2, 'Coefficient vs Reconstruction');
legend(ax2, {'Direct coeff', 'Direct recon', 'Spectral coeff', 'Spectral recon'}, ...
    'Location', 'northwest', 'Box', 'off');

ax3 = nexttile;
yyaxis(ax3, 'left');
plot(ax3, T.N, T.speedup_direct_over_spectral, '-d', ...
    'Color', style.accent_color, 'LineWidth', 2.2, ...
    'MarkerSize', 7, 'MarkerFaceColor', style.accent_color);
ylabel(ax3, 'Speedup (direct / spectral)');

yyaxis(ax3, 'right');
plot(ax3, T.N, T.eta_max_abs_err, ':^', ...
    'Color', style.error_color, 'LineWidth', 1.8, ...
    'MarkerSize', 6, 'MarkerFaceColor', style.error_color);
ylabel(ax3, 'Max |eta_{direct} - eta_{spectral}|');

format_axes(ax3, style);
set(ax3, 'YScale', 'log');
xlabel(ax3, 'N components');
title(ax3, 'Speedup and Consistency');
legend(ax3, {'Speedup', 'Max |eta error|'}, 'Location', 'northwest', 'Box', 'off');

sgtitle(fig, sprintf('MF12 Direct vs Spectral Benchmark  |  Grid %dx%d  |  t = %.2f s', ...
    plotData.cfg.Nx, plotData.cfg.Ny, plotData.cfg.t), ...
    'FontName', style.font_name, 'FontSize', 15, 'FontWeight', 'bold');
end

function style = benchmark_plot_style()
style = struct();
style.figure_color = [1.0 1.0 1.0];
style.axes_color = [0.98 0.985 0.99];
style.grid_color = [0.83 0.86 0.90];
style.direct_color = [0.07 0.37 0.70];
style.direct_light = [0.42 0.63 0.86];
style.spectral_color = [0.82 0.33 0.18];
style.spectral_light = [0.92 0.62 0.48];
style.accent_color = [0.12 0.56 0.38];
style.error_color = [0.58 0.22 0.67];
style.font_name = 'Helvetica';
end

function format_axes(ax, style)
set(ax, ...
    'Color', style.axes_color, ...
    'FontName', style.font_name, ...
    'FontSize', 11, ...
    'LineWidth', 0.9, ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'GridColor', style.grid_color, ...
    'GridAlpha', 0.55, ...
    'MinorGridColor', style.grid_color, ...
    'MinorGridAlpha', 0.25);
grid(ax, 'on');
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'on';
end
