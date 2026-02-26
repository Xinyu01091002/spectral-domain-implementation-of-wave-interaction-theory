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
cfg.seed = 20260225;

% [N_components, Nx, Ny]
% Keep moderate sizes because direct full-order pipeline is expensive.
case_matrix = [
    20, 64, 64;
    40, 64, 64;
    60, 64, 64;
    80, 64, 64
];

rng(cfg.seed);
rows = [];
k = 0;

fprintf('=== MF12 Speed/Memory Benchmark ===\n');
fprintf('Pipelines:\n');
fprintf('  Direct   = coeffsMF12 + surfaceMF12_new\n');
fprintf('  Spectral = coeffsMF12_superharmonic + surfaceMF12_spectral\n');
fprintf('  Streaming= surfaceMF12_integrated_opt(useStreaming=true)\n');
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

    t_direct = zeros(cfg.repeats, 1);
    t_spec = zeros(cfg.repeats, 1);
    t_stream = zeros(cfg.repeats, 1);

    m_direct = zeros(cfg.repeats, 1);
    m_spec = zeros(cfg.repeats, 1);
    m_stream = zeros(cfg.repeats, 1);

    eta_direct_ref = [];
    eta_spec_ref = [];
    eta_stream_ref = [];

    % Warmup (once per case)
    coeffs_d = coeffsMF12(cfg.order, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy);
    [~, ~] = surfaceMF12_new(cfg.order, coeffs_d, X, Y, cfg.t);
    coeffs_s = coeffsMF12_superharmonic(cfg.order, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy);
    [~, ~] = surfaceMF12_spectral(coeffs_s, cfg.Lx, cfg.Ly, Nx, Ny, cfg.t);
    opts = struct('useStreaming', true, 'legacyMode', false);
    [~, ~] = surfaceMF12_integrated_opt(cfg.order, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy, cfg.Lx, cfg.Ly, Nx, Ny, cfg.t, opts);

    for r = 1:cfg.repeats
        % Direct (full)
        m0 = safe_memory_snapshot();
        tic;
        coeffs_d = coeffsMF12(cfg.order, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy);
        [eta_d, ~] = surfaceMF12_new(cfg.order, coeffs_d, X, Y, cfg.t);
        t_direct(r) = toc;
        m1 = safe_memory_snapshot();
        m_direct(r) = mem_delta_mb(m0, m1);

        % Spectral (superharmonic coeffs + spectral reconstruction)
        m0 = safe_memory_snapshot();
        tic;
        coeffs_s = coeffsMF12_superharmonic(cfg.order, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy);
        [eta_s, ~] = surfaceMF12_spectral(coeffs_s, cfg.Lx, cfg.Ly, Nx, Ny, cfg.t);
        t_spec(r) = toc;
        m1 = safe_memory_snapshot();
        m_spec(r) = mem_delta_mb(m0, m1);

        % Streaming spectral (integrated)
        m0 = safe_memory_snapshot();
        tic;
        [eta_st, ~] = surfaceMF12_integrated_opt(cfg.order, cfg.g, cfg.h, a, b, kx, ky, cfg.Ux, cfg.Uy, cfg.Lx, cfg.Ly, Nx, Ny, cfg.t, opts);
        t_stream(r) = toc;
        m1 = safe_memory_snapshot();
        m_stream(r) = mem_delta_mb(m0, m1);

        if r == 1
            eta_direct_ref = eta_d;
            eta_spec_ref = eta_s;
            eta_stream_ref = eta_st;
        end
    end

    [eta_spec_vs_direct_max_err, eta_spec_vs_direct_nan] = max_abs_err_finite(eta_spec_ref, eta_direct_ref);
    [eta_stream_vs_direct_max_err, eta_stream_vs_direct_nan] = max_abs_err_finite(eta_stream_ref, eta_direct_ref);
    [eta_stream_vs_spec_max_err, eta_stream_vs_spec_nan] = max_abs_err_finite(eta_stream_ref, eta_spec_ref);

    k = k + 1;
    rows(k).case_id = c; %#ok<SAGROW>
    rows(k).N = N;
    rows(k).Nx = Nx;
    rows(k).Ny = Ny;

    rows(k).direct_time_s = mean(t_direct);
    rows(k).spectral_time_s = mean(t_spec);
    rows(k).streaming_time_s = mean(t_stream);

    rows(k).direct_mem_mb = mean(m_direct, 'omitnan');
    rows(k).spectral_mem_mb = mean(m_spec, 'omitnan');
    rows(k).streaming_mem_mb = mean(m_stream, 'omitnan');

    rows(k).speedup_direct_over_spectral = rows(k).direct_time_s / rows(k).spectral_time_s;
    rows(k).speedup_spectral_over_streaming = rows(k).spectral_time_s / rows(k).streaming_time_s;
    rows(k).eta_spec_vs_direct_max_err = eta_spec_vs_direct_max_err;
    rows(k).eta_spec_vs_direct_nan = eta_spec_vs_direct_nan;
    rows(k).eta_stream_vs_direct_max_err = eta_stream_vs_direct_max_err;
    rows(k).eta_stream_vs_direct_nan = eta_stream_vs_direct_nan;
    rows(k).eta_stream_vs_spec_max_err = eta_stream_vs_spec_max_err;
    rows(k).eta_stream_vs_spec_nan = eta_stream_vs_spec_nan;

    fprintf('  time(s): direct %.3f | spectral %.3f | streaming %.3f\n', ...
        rows(k).direct_time_s, rows(k).spectral_time_s, rows(k).streaming_time_s);
    fprintf('  mem(MB): direct %.2f | spectral %.2f | streaming %.2f\n', ...
        rows(k).direct_mem_mb, rows(k).spectral_mem_mb, rows(k).streaming_mem_mb);
    fprintf('  eta err max: spec-direct %.3e | stream-direct %.3e | stream-spec %.3e\n', ...
        rows(k).eta_spec_vs_direct_max_err, rows(k).eta_stream_vs_direct_max_err, rows(k).eta_stream_vs_spec_max_err);
end

T = struct2table(rows);
disp(T);

ts = datestr(now, 'yyyymmdd_HHMMSS');
csv_name = ['benchmark_mf12_speed_memory_' ts '.csv'];
mat_name = ['benchmark_mf12_speed_memory_' ts '.mat'];
png_name = ['benchmark_mf12_speed_memory_' ts '.png'];

writetable(T, csv_name);
save(mat_name, 'T', 'rows', 'cfg', 'case_matrix');

% Plot
c_direct = [0.00, 0.33, 0.62];
c_spec = [0.84, 0.37, 0.00];
c_stream = [0.15, 0.60, 0.25];

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 7.5]);
tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile;
semilogy(T.N, T.direct_time_s, '-o', 'Color', c_direct, 'LineWidth', 2, 'MarkerSize', 5); hold on;
semilogy(T.N, T.spectral_time_s, '-s', 'Color', c_spec, 'LineWidth', 2, 'MarkerSize', 5);
semilogy(T.N, T.streaming_time_s, '-^', 'Color', c_stream, 'LineWidth', 2, 'MarkerSize', 5);
grid on; grid minor;
xlabel('Number of wavenumber components', 'FontName', 'Times New Roman', 'FontSize', 10);
ylabel('Runtime (s)', 'FontName', 'Times New Roman', 'FontSize', 10);
set(ax1, 'FontName', 'Times New Roman', 'FontSize', 9, 'LineWidth', 0.8, 'Box', 'on');
legend({'Direct', 'Spectral', 'Streaming'}, 'Location', 'northwest', ...
    'FontName', 'Times New Roman', 'FontSize', 8, 'Box', 'on', 'Color', 'w');
text(0.97, 0.95, '(a)', 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'FontName', 'Times New Roman', 'FontSize', 11, ...
    'FontWeight', 'bold', 'BackgroundColor', 'w', 'Margin', 1);

ax2 = nexttile;
semilogy(T.N, max(T.direct_mem_mb, eps), '-o', 'Color', c_direct, 'LineWidth', 2, 'MarkerSize', 5); hold on;
semilogy(T.N, max(T.spectral_mem_mb, eps), '-s', 'Color', c_spec, 'LineWidth', 2, 'MarkerSize', 5);
semilogy(T.N, max(T.streaming_mem_mb, eps), '-^', 'Color', c_stream, 'LineWidth', 2, 'MarkerSize', 5);
grid on; grid minor;
xlabel('Number of wavenumber components', 'FontName', 'Times New Roman', 'FontSize', 10);
ylabel('Memory delta (MB)', 'FontName', 'Times New Roman', 'FontSize', 10);
set(ax2, 'FontName', 'Times New Roman', 'FontSize', 9, 'LineWidth', 0.8, 'Box', 'on');
legend({'Direct', 'Spectral', 'Streaming'}, 'Location', 'northwest', ...
    'FontName', 'Times New Roman', 'FontSize', 8, 'Box', 'on', 'Color', 'w');
text(0.97, 0.95, '(b)', 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'FontName', 'Times New Roman', 'FontSize', 11, ...
    'FontWeight', 'bold', 'BackgroundColor', 'w', 'Margin', 1);

exportgraphics(fig, png_name, 'Resolution', 600);

fprintf('\nSaved: %s\n', csv_name);
fprintf('Saved: %s\n', mat_name);
fprintf('Saved: %s\n', png_name);
fprintf('Note: memory deltas are process-level approximations, not true peak memory.\n');

function [a, b, kx, ky] = build_case_components(N, Lx, Ly, Nx, Ny)
    dkx = 2*pi/Lx;
    dky = 2*pi/Ly;

    kx_idx_all = (-floor(Nx/2)):(ceil(Nx/2)-1);
    ky_idx_all = (-floor(Ny/2)):(ceil(Ny/2)-1);
    [KXI, KYI] = meshgrid(kx_idx_all, ky_idx_all);
    KXI = KXI(:);
    KYI = KYI(:);
    keep = (KXI > 0) | (KXI == 0 & KYI > 0);
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

function dmb = mem_delta_mb(m0, m1)
    if ~isempty(m0) && ~isempty(m1)
        dmb = (m1 - m0) / (1024^2);
    else
        dmb = NaN;
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
