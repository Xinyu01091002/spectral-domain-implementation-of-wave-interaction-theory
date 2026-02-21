clc; clear; close all;
addpath(genpath(fullfile(pwd, 'irregularWavesMF12')));

cfg.g = 9.81;
cfg.h = 100;
cfg.Ux = 0;
cfg.Uy = 0;
cfg.Lx = 3000;
cfg.Ly = 3000;
cfg.t = 0.0;
cfg.order = 3;
cfg.repeats = 1;
cfg.seed = 20260221;

% [N_components, Nx, Ny]
case_matrix = [
    20, 64, 64;
    40, 64, 64;
    60, 64, 64;
    100, 64,64
];

rng(cfg.seed);
rows = [];
k = 0;

fprintf('=== Direct MF12 vs Spectral MF12 ===\n');
fprintf('Order=%d, repeats=%d\n\n', cfg.order, cfg.repeats);

for c = 1:size(case_matrix, 1)
    N = case_matrix(c, 1);
    Nx = case_matrix(c, 2);
    Ny = case_matrix(c, 3);
    [a, b, kx, ky] = build_case_components(N, cfg.Lx, cfg.Ly, Nx, Ny);

    dx = cfg.Lx / Nx;
    dy = cfg.Ly / Ny;
    x = (0:Nx-1) * dx;
    y = (0:Ny-1) * dy;
    [X, Y] = meshgrid(x, y);

    fprintf('Case %d: N=%d, grid=%dx%d\n', c, N, Nx, Ny);

    % Shared coefficients to isolate reconstruction speed.
    % Use full coeffsMF12 so surfaceMF12_new indexing is fully compatible.
    coeffs = coeffsMF12(cfg.order, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy);

    t_direct = zeros(cfg.repeats, 1);
    t_spectral = zeros(cfg.repeats, 1);
    mem_direct = zeros(cfg.repeats, 1);
    mem_spectral = zeros(cfg.repeats, 1);

    eta_direct_ref = [];
    phi_direct_ref = [];
    eta_spec_ref = [];
    phi_spec_ref = [];

    % Warmup
    [~, ~] = surfaceMF12_new(cfg.order, coeffs, X, Y, cfg.t);
    [~, ~] = surfaceMF12_spectral(coeffs, cfg.Lx, cfg.Ly, Nx, Ny, cfg.t);

    for r = 1:cfg.repeats
        m0 = safe_memory_snapshot();
        tic;
        [eta_d, phi_d] = surfaceMF12_new(cfg.order, coeffs, X, Y, cfg.t);
        t_direct(r) = toc;
        m1 = safe_memory_snapshot();
        if ~isempty(m0) && ~isempty(m1), mem_direct(r) = (m1 - m0)/(1024^2); else, mem_direct(r) = NaN; end

        m0 = safe_memory_snapshot();
        tic;
        [eta_s, phi_s] = surfaceMF12_spectral(coeffs, cfg.Lx, cfg.Ly, Nx, Ny, cfg.t);
        t_spectral(r) = toc;
        m1 = safe_memory_snapshot();
        if ~isempty(m0) && ~isempty(m1), mem_spectral(r) = (m1 - m0)/(1024^2); else, mem_spectral(r) = NaN; end

        if r == 1
            eta_direct_ref = eta_d;
            phi_direct_ref = phi_d;
            eta_spec_ref = eta_s;
            phi_spec_ref = phi_s;
        end
    end

    [eta_max_err, eta_nan_count] = max_abs_err_finite(eta_direct_ref, eta_spec_ref);
    [phi_max_err, phi_nan_count] = max_abs_err_finite(phi_direct_ref, phi_spec_ref);
    speedup = mean(t_direct) / mean(t_spectral);

    k = k + 1;
    rows(k).case_id = c; %#ok<SAGROW>
    rows(k).N = N;
    rows(k).Nx = Nx;
    rows(k).Ny = Ny;
    rows(k).direct_mean_s = mean(t_direct);
    rows(k).direct_min_s = min(t_direct);
    rows(k).spectral_mean_s = mean(t_spectral);
    rows(k).spectral_min_s = min(t_spectral);
    rows(k).speedup_direct_over_spectral = speedup;
    rows(k).direct_mem_mean_mb = mean(mem_direct, 'omitnan');
    rows(k).spectral_mem_mean_mb = mean(mem_spectral, 'omitnan');
    rows(k).eta_max_abs_err = eta_max_err;
    rows(k).phi_max_abs_err = phi_max_err;
    rows(k).eta_nan_count = eta_nan_count;
    rows(k).phi_nan_count = phi_nan_count;

    fprintf('  direct mean %.3fs | spectral mean %.3fs | speedup %.2fx | eta err %.3e | eta NaN %d\n', ...
        rows(k).direct_mean_s, rows(k).spectral_mean_s, rows(k).speedup_direct_over_spectral, rows(k).eta_max_abs_err, rows(k).eta_nan_count);
end

T = struct2table(rows);
disp(T);

ts = datestr(now, 'yyyymmdd_HHMMSS');
csv_name = ['benchmark_direct_vs_spectral_' ts '.csv'];
mat_name = ['benchmark_direct_vs_spectral_' ts '.mat'];

writetable(T, csv_name);
save(mat_name, 'T', 'rows', 'cfg', 'case_matrix');

% Quick plot
fig = figure('Color', 'w', 'Position', [100 100 1200 420]);
tiledlayout(1,3);

nexttile;
bar([T.direct_mean_s, T.spectral_mean_s], 'grouped');
grid on;
title('Mean Runtime');
ylabel('s');
set(gca, 'YScale', 'log');
xticklabels("N="+string(T.N)+", "+string(T.Nx)+"x"+string(T.Ny));
xtickangle(15);
legend({'direct', 'spectral'}, 'Location', 'northwest');

nexttile;
bar(T.speedup_direct_over_spectral);
grid on;
title('Speedup (direct/spectral)');
ylabel('x');
xticklabels("N="+string(T.N));

nexttile;
err_plot = T.eta_max_abs_err;
err_plot(err_plot <= 0) = 1e-18;
bar(err_plot);
set(gca, 'YScale', 'log');
grid on;
title('Max |eta_{direct}-eta_{spectral}|');
ylabel('log scale');
xticklabels("N="+string(T.N));

png_name = ['benchmark_direct_vs_spectral_' ts '.png'];
saveas(fig, png_name);

fprintf('\nSaved: %s\n', csv_name);
fprintf('Saved: %s\n', mat_name);
fprintf('Saved: %s\n', png_name);
fprintf('Note: memory deltas are process-level approximations, not true peak memory.\n');

function [a, b, kx, ky] = build_case_components(N, Lx, Ly, Nx, Ny)
    dkx = 2*pi/Lx;
    dky = 2*pi/Ly;

    % Build unique spectral candidates from one half-plane to avoid duplicates.
    kx_idx_all = (-floor(Nx/2)):(ceil(Nx/2)-1);
    ky_idx_all = (-floor(Ny/2)):(ceil(Ny/2)-1);
    [KXI, KYI] = meshgrid(kx_idx_all, ky_idx_all);
    KXI = KXI(:);
    KYI = KYI(:);
    keep = (KXI > 0) | (KXI == 0 & KYI > 0); % half-plane, exclude (0,0)
    KXI = KXI(keep);
    KYI = KYI(keep);

    nCand = numel(KXI);
    if N > nCand
        error('N=%d exceeds available unique half-plane modes=%d for grid %dx%d.', N, nCand, Nx, Ny);
    end

    pick = randperm(nCand, N).';
    kx = KXI(pick) * dkx;
    ky = KYI(pick) * dky;
    kmag = hypot(kx, ky);
    amp = 0.03 * exp(-0.2 * (kmag / max(kmag + eps)));
    phase = 2*pi*rand(N, 1);
    a = amp .* cos(phase);
    b = amp .* sin(phase);
end

function mem_used = safe_memory_snapshot()
    mem_used = [];
    try
        m = memory;
        mem_used = m.MemUsedMATLAB;
    catch
        mem_used = [];
    end
end

function [err, nan_count] = max_abs_err_finite(A, B)
    d = A - B;
    finite_mask = isfinite(d);
    nan_count = numel(d) - nnz(finite_mask);
    if any(finite_mask(:))
        v = abs(d(finite_mask));
        err = max(v(:));
    else
        err = NaN;
    end
end
