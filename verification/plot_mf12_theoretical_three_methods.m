clc; clear; close all;
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
outDir = fullfile(rootDir, 'outputs');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Theoretical end-to-end comparison (no hardware timing):
%   Direct   = coeffsMF12 + surfaceMF12_new
%   Spectral = coeffsMF12_superharmonic + surfaceMF12_spectral
%   Streaming= surfaceMF12_streaming_super (integrated on-the-fly)
%
% Notes:
% 1) This is a structural complexity model from loop/data shapes in code.
% 2) Curves are "relative units", not measured seconds.
% 3) Full polynomial terms are kept (not only leading order).

cfg.order = 3;
cfg.Nx = 512;
cfg.Ny = 512;
cfg.N_list = unique(round(logspace(log10(8), log10(2000), 140)));

G = cfg.Nx * cfg.Ny;
N = cfg.N_list(:);
C2 = N .* (N - 1) / 2;
C3 = N .* (N - 1) .* (N - 2) / 6;
L2 = N .* (N - 1);  % = 2*C2

% ------------------------
% Term counts by pipeline
% ------------------------
% Precompute (full):
M_pre_direct = N + 3*L2 + 4*C3;                % coeffsMF12

% Precompute (superharmonic):
M_pre_spec = N + (2*L2 + C2) + C3;             % coeffsMF12_superharmonic

% Direct reconstruction loop count:
M_recon_direct = 3*N + 2*L2 + 4*C3;            % surfaceMF12_new

% Spectral reconstruction term count:
M_recon_spec = 3*N + 2*L2 + C3;                % surfaceMF12_spectral(super-only)

% Streaming on-the-fly coefficient/transfer calculations:
% Includes extra pair-correction loops in streaming source.
M_stream_calc = 5*N + 3*L2 + C3;               % surfaceMF12_streaming_super

% ------------------------
% Time model (relative)
% ------------------------
% Constants are model weights to combine unlike operations.
k_fft = 5;               % ifft2 weight
k_dep = 8;               % non-stream spectral bilinear deposition
k_stream_grid = 0.20;    % streaming add_to_spec accumarray([Ny*Nx,1]) per call

% Direct: every reconstruction term acts on full grid.
time_direct = M_pre_direct + G .* M_recon_direct;

% Spectral: precompute + term deposition + one FFT.
time_spectral = M_pre_spec + k_dep .* M_recon_spec + k_fft .* G .* log2(G);

% Streaming: no large coeff arrays; on-the-fly transfer + per-call full-grid accumarray + FFT.
time_stream = M_stream_calc + k_stream_grid .* G .* M_recon_spec + k_fft .* G .* log2(G);

% ------------------------
% Memory model (bytes)
% ------------------------
bytes_real = 8;
bytes_complex = 16;

% Direct precompute coefficient storage (rough field-count model).
coeff_fields_full = 31*N + 30*L2 + 36*C3 + N.^2;
mem_pre_direct = coeff_fields_full .* bytes_real;

% Spectral precompute coefficient storage (superharmonic reduced).
coeff_fields_super = 31*N + 18*L2 + 9*C3 + 2*N.^2;
mem_pre_spec = coeff_fields_super .* bytes_real;

% Direct reconstruction peak memory (grid dominated).
mem_recon_direct = (8 * bytes_real * G) .* ones(size(N));

% Non-stream spectral reconstruction peak:
% base grids + large idx/value accumulation vectors.
mem_spec_base = (2*bytes_complex + 5*bytes_real) * G;
mem_spec_idxvals = (4*(4*M_recon_spec)*bytes_real) + (2*(4*M_recon_spec)*bytes_complex);
mem_recon_spec = mem_spec_base + mem_spec_idxvals;

% Streaming reconstruction peak:
% mainly fixed spectral grids and transient full-grid addv buffer, plus O(N) arrays.
mem_stream_base = (2*bytes_complex + 6*bytes_real) * G;
mem_stream_temp = (1*bytes_complex + 1*bytes_real) * G;
mem_stream_linearN = (40 * bytes_real) .* N;
mem_recon_stream = mem_stream_base + mem_stream_temp + mem_stream_linearN;

mem_direct_mb = (mem_pre_direct + mem_recon_direct) / 1024^2;
mem_spec_mb = (mem_pre_spec + mem_recon_spec) / 1024^2;
mem_stream_mb = (mem_recon_stream) / 1024^2;

% ------------------------
% Plot
% ------------------------
c_direct = [0.00, 0.33, 0.62];
c_spec = [0.84, 0.37, 0.00];
c_stream = [0.15, 0.60, 0.25];

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 7.8]);
tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile;
loglog(ax1, N, time_direct, '-', 'Color', c_direct, 'LineWidth', 2.2); hold(ax1, 'on');
loglog(ax1, N, time_spectral, '-', 'Color', c_spec, 'LineWidth', 2.2);
loglog(ax1, N, time_stream, '-', 'Color', c_stream, 'LineWidth', 2.2);
grid(ax1, 'on'); grid(ax1, 'minor');
xlabel(ax1, 'Number of wavenumber components', 'FontName', 'Times New Roman', 'FontSize', 10);
ylabel(ax1, 'Relative operation count', 'FontName', 'Times New Roman', 'FontSize', 10);
set(ax1, 'FontName', 'Times New Roman', 'FontSize', 9, 'LineWidth', 0.8, 'Box', 'on');
legend(ax1, {'Direct', 'Spectral', 'Streaming'}, ...
    'Location', 'northwest', 'FontName', 'Times New Roman', 'FontSize', 8, ...
    'Box', 'on', 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
text(ax1, 0.97, 0.95, '(a)', 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'FontName', 'Times New Roman', 'FontSize', 11, ...
    'FontWeight', 'bold', 'BackgroundColor', 'w', 'Margin', 1);

ax2 = nexttile;
semilogy(ax2, N, mem_direct_mb, '-', 'Color', c_direct, 'LineWidth', 2.2); hold(ax2, 'on');
semilogy(ax2, N, mem_spec_mb, '-', 'Color', c_spec, 'LineWidth', 2.2);
semilogy(ax2, N, mem_stream_mb, '-', 'Color', c_stream, 'LineWidth', 2.2);
grid(ax2, 'on'); grid(ax2, 'minor');
xlabel(ax2, 'Number of wavenumber components', 'FontName', 'Times New Roman', 'FontSize', 10);
ylabel(ax2, 'Estimated peak memory (MB)', 'FontName', 'Times New Roman', 'FontSize', 10);
set(ax2, 'FontName', 'Times New Roman', 'FontSize', 9, 'LineWidth', 0.8, 'Box', 'on');
legend(ax2, {'Direct', 'Spectral', 'Streaming'}, ...
    'Location', 'northwest', 'FontName', 'Times New Roman', 'FontSize', 8, ...
    'Box', 'on', 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
text(ax2, 0.97, 0.95, '(b)', 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'FontName', 'Times New Roman', 'FontSize', 11, ...
    'FontWeight', 'bold', 'BackgroundColor', 'w', 'Margin', 1);

ts = datestr(now, 'yyyymmdd_HHMMSS');
out_png = fullfile(outDir, ['mf12_theory_three_methods_' ts '.png']);
exportgraphics(fig, out_png, 'Resolution', 600);

fprintf('Saved figure: %s\n', out_png);
fprintf('Model: full polynomial term counts included, order=%d, grid=%dx%d.\n', cfg.order, cfg.Nx, cfg.Ny);
fprintf('At N=%d (MB): Direct=%.1f, Spectral=%.1f, Streaming=%.1f\n', ...
    N(end), mem_direct_mb(end), mem_spec_mb(end), mem_stream_mb(end));
fprintf('At N=%d (time units): Direct=%.3e, Spectral=%.3e, Streaming=%.3e\n', ...
    N(end), time_direct(end), time_spectral(end), time_stream(end));
