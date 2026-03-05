% BENCHMARK_DIRECT_VS_SPECTRAL_VSN
% Compare reconstruction runtime vs N:
%   Direct sum:    surfaceMF12_new
%   Spectral (new): surfaceMF12_spectral

clc;
clear;
close all;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(genpath(fullfile(rootDir, 'irregularWavesMF12')));
outDir = fullfile(rootDir, 'outputs');
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

cfg.N_list = [40 60 80 100 150];
cfg.repeats = 1;
cfg.seed = 20260306;

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

fprintf('=== Direct vs Spectral runtime vs N ===\n');
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

    coeffs = coeffsMF12(cfg.order, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy);

    % Warmup
    surfaceMF12_new(cfg.order, coeffs, X, Y, cfg.t);
    coeffs_w = coeffs;
    coeffs_w.third_order_subharmonic_mode = 'skip';
    surfaceMF12_spectral(coeffs_w, cfg.Lx, cfg.Ly, cfg.Nx, cfg.Ny, cfg.t);

    t_direct = zeros(cfg.repeats, 1);
    t_spectral = zeros(cfg.repeats, 1);
    eta_ref = [];
    eta_spec = [];

    for r = 1:cfg.repeats
        tic;
        [eta_d, ~] = surfaceMF12_new(cfg.order, coeffs, X, Y, cfg.t);
        t_direct(r) = toc;

        coeffs_s = coeffs;
        coeffs_s.third_order_subharmonic_mode = 'skip';
        tic;
        [eta_s, ~] = surfaceMF12_spectral(coeffs_s, cfg.Lx, cfg.Ly, cfg.Nx, cfg.Ny, cfg.t);
        t_spectral(r) = toc;

        if r == 1
            eta_ref = eta_d;
            eta_spec = eta_s;
        end
    end

    k = k + 1;
    rows(k).N = N; %#ok<SAGROW>
    rows(k).direct_mean_s = mean(t_direct);
    rows(k).spectral_mean_s = mean(t_spectral);
    rows(k).speedup_direct_over_spectral = rows(k).direct_mean_s / max(rows(k).spectral_mean_s, eps);
    rows(k).eta_max_abs_err = max(abs(eta_ref(:) - eta_spec(:)));

    fprintf('N=%3d | direct %.3fs | spectral %.3fs | speedup %.2fx | max|deta| %.3e\n', ...
        N, rows(k).direct_mean_s, rows(k).spectral_mean_s, ...
        rows(k).speedup_direct_over_spectral, rows(k).eta_max_abs_err);
end

T = struct2table(rows);
disp(T);

ts = datestr(now, 'yyyymmdd_HHMMSS');
csv_name = fullfile(outDir, ['benchmark_direct_vs_spectral_vsN_' ts '.csv']);
png_name = fullfile(outDir, ['benchmark_direct_vs_spectral_vsN_' ts '.png']);
writetable(T, csv_name);

fig = figure('Color', 'w', 'Position', [100 100 980 420]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile;
plot(T.N, T.direct_mean_s, 'o-', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
plot(T.N, T.spectral_mean_s, 's-', 'LineWidth', 1.8, 'MarkerSize', 7); hold off;
grid on;
xlabel('N components');
ylabel('Mean runtime (s)');
title('Runtime vs N');
set(gca, 'YScale', 'log');
legend({'Direct sum', 'Spectral'}, 'Location', 'northwest');

nexttile;
plot(T.N, T.speedup_direct_over_spectral, 'd-', 'LineWidth', 1.8, 'MarkerSize', 7);
grid on;
xlabel('N components');
ylabel('Speedup (direct / spectral)');
title('Speedup vs N');

exportgraphics(fig, png_name, 'Resolution', 160);

fprintf('\nSaved: %s\n', csv_name);
fprintf('Saved: %s\n', png_name);
