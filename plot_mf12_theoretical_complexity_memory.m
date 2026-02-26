clc; clear; close all;

% End-to-end theoretical scaling model:
%   Direct   : coeffsMF12 (full) + surfaceMF12_new
%   Spectral : coeffsMF12_superharmonic + surfaceMF12_spectral

cfg.order = 3;
cfg.Nx = 512;
cfg.Ny = 512;
cfg.N_list = unique(round(logspace(log10(8), log10(2000), 120)));

G = cfg.Nx * cfg.Ny;
N = cfg.N_list(:);

% Combinatorics
C2 = N .* (N - 1) / 2;
C3 = N .* (N - 1) .* (N - 2) / 6;
L2 = N .* (N - 1);      % = 2*C2

% ---------- Reconstruction term counts ----------
% Direct surfaceMF12_new loop structure (includes loops even when some terms are skipped by if):
M_recon_direct = ...
    N + ...             % 1st order
    N + ...             % 2nd self
    L2 + ...            % 2nd pair (+/-)
    N + ...             % 3rd self
    L2 + ...            % 3rd double loop over +/- m
    4*C3;               % 3rd triple loop over +/-m and +/-p

% Spectral surfaceMF12_spectral with superharmonic coeffs:
M_recon_spec = ...
    N + ...             % 1st order
    N + ...             % 2nd self
    L2 + ...            % 2nd pair (+/-)
    N + ...             % 3rd self
    L2 + ...            % 3rd double: (n+2m) and (2n+m)
    C3;                 % 3rd triple: (n+m+p) only

% ---------- Precompute loop counts ----------
% coeffsMF12 (full)
M_pre_direct = ...
    N + ...             % O(N) single loops
    3*L2 + ...          % three O(N^2) loop groups
    4*C3;               % full +/- triple combinations

% coeffsMF12_superharmonic
M_pre_spec = ...
    N + ...             % O(N) single loops
    (2*L2 + C2) + ...   % second pair + third correction + super double
    C3;                 % super triple only

% ---------- Time model (relative units) ----------
time_direct_units = M_pre_direct + G .* M_recon_direct;
deposition_units = 8 .* M_recon_spec;          % eta+phi bilinear deposition
fft_units = 5 .* G .* log2(G);                 % visualization constant factor
time_spectral_units = M_pre_spec + deposition_units + fft_units;

% ---------- Memory model ----------
bytes_real = 8;
bytes_complex = 16;

% Reconstruction memory
direct_grid_fields = 8;                        % rough peak full-size real arrays
mem_recon_direct_bytes = direct_grid_fields * bytes_real * G * ones(size(N));

spec_complex_grids = 2;                    % spec_eta, spec_phi
spec_real_grids = 5;                       % eta, phi_wave, phiS, X, Y
mem_recon_spectral_base = spec_complex_grids * bytes_complex * G + ...
                          spec_real_grids * bytes_real * G;
idx_arrays_bytes = 4 .* (4 .* M_recon_spec) .* bytes_real;
val_arrays_bytes = 2 .* (4 .* M_recon_spec) .* bytes_complex;
mem_recon_spectral_bytes = mem_recon_spectral_base + idx_arrays_bytes + val_arrays_bytes;

% Precompute memory (rough coefficient-array storage + key temporary maps)
% Full coeffsMF12: many O(N^2) arrays and full 4*C3 triple arrays.
coeff_fields_full = ...
    31*N + ...            % linear + self terms + 3rd singles (rough)
    30*L2 + ...           % pair and double families (rough)
    36*C3 + ...           % triple families (A/B/F/G/mu/kappa/omega/kx/ky with 4*C3)
    N.^2;                 % M_nm-like indexing map
mem_pre_direct_bytes = coeff_fields_full .* bytes_real;

% Superharmonic coeffs: reduced double/triple storage, but includes two N^2 maps.
coeff_fields_super = ...
    31*N + ...            % linear + self terms + 3rd singles (rough)
    18*L2 + ...           % reduced pair/double families
    9*C3 + ...            % super triple only
    2*N.^2;               % M_nm_row_odd + M_nm_col_odd
mem_pre_spectral_bytes = coeff_fields_super .* bytes_real;

mem_direct_mb = (mem_pre_direct_bytes + mem_recon_direct_bytes) / 1024^2;
mem_spectral_mb = (mem_pre_spectral_bytes + mem_recon_spectral_bytes) / 1024^2;
mem_pre_direct_mb = mem_pre_direct_bytes / 1024^2;
mem_recon_direct_mb = mem_recon_direct_bytes / 1024^2;
mem_pre_spectral_mb = mem_pre_spectral_bytes / 1024^2;
mem_recon_spectral_mb = mem_recon_spectral_bytes / 1024^2;

% ---------- Sanity checks ----------
assert(all(diff(N) > 0), 'N_list must be strictly increasing.');
assert(all(diff(time_direct_units) >= 0) && all(time_direct_units > 0), 'time_direct must be positive and nondecreasing.');
assert(all(diff(time_spectral_units) >= 0) && all(time_spectral_units > 0), 'time_spectral must be positive and nondecreasing.');
assert(all(diff(mem_spectral_mb) >= 0) && all(mem_spectral_mb > 0), 'mem_spectral must be positive and nondecreasing.');
assert(all(mem_direct_mb > 0), 'mem_direct must be positive.');

% ---------- Publication-style plotting ----------
c_direct = [0.00, 0.33, 0.62];
c_spec = [0.84, 0.37, 0.00];
lw = 2.2;
ms = 5;

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 7.5]);
tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile;
loglog(ax1, N, time_direct_units, '-o', 'Color', c_direct, 'LineWidth', lw, ...
    'MarkerSize', ms, 'MarkerIndices', 1:6:numel(N)); hold(ax1, 'on');
loglog(ax1, N, time_spectral_units, '-s', 'Color', c_spec, 'LineWidth', lw, ...
    'MarkerSize', ms, 'MarkerIndices', 1:6:numel(N));
grid(ax1, 'on'); grid(ax1, 'minor');
xlabel(ax1, 'Number of wavenumber components', 'FontName', 'Times New Roman', 'FontSize', 10);
ylabel(ax1, 'Relative operation count', 'FontName', 'Times New Roman', 'FontSize', 10);
set(ax1, 'FontName', 'Times New Roman', 'FontSize', 9, 'LineWidth', 0.8, 'Box', 'on');
legend(ax1, {'Direct', 'Spectral'}, ...
    'Location', 'northwest', 'FontName', 'Times New Roman', 'FontSize', 8, ...
    'Box', 'on', 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
text(ax1, 0.97, 0.95, '(a)', 'Units', 'normalized', 'FontName', 'Times New Roman', ...
    'FontSize', 11, 'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'right', ...
    'BackgroundColor', 'w', 'Margin', 1);

ax2 = nexttile;
hold(ax2, 'on');
hDT = semilogy(ax2, N, mem_direct_mb, '-', 'Color', c_direct, 'LineWidth', 2.4);
hST = semilogy(ax2, N, mem_spectral_mb, '-', 'Color', c_spec, 'LineWidth', 2.4);

% Decomposition curves (same color family, lighter line styles)
semilogy(ax2, N, mem_pre_direct_mb, '--', 'Color', c_direct, 'LineWidth', 1.4);
semilogy(ax2, N, mem_recon_direct_mb, ':', 'Color', c_direct, 'LineWidth', 1.4);
semilogy(ax2, N, mem_pre_spectral_mb, '--', 'Color', c_spec, 'LineWidth', 1.4);
semilogy(ax2, N, mem_recon_spectral_mb, ':', 'Color', c_spec, 'LineWidth', 1.4);

grid(ax2, 'on'); grid(ax2, 'minor');
xlabel(ax2, 'Number of wavenumber components', 'FontName', 'Times New Roman', 'FontSize', 10);
ylabel(ax2, 'Estimated peak memory (MB)', 'FontName', 'Times New Roman', 'FontSize', 10);
set(ax2, 'FontName', 'Times New Roman', 'FontSize', 9, 'LineWidth', 0.8, 'Box', 'on');
% Compact in-axis legend: totals + style meaning
hPre = plot(ax2, NaN, NaN, 'k--', 'LineWidth', 1.4);
hRec = plot(ax2, NaN, NaN, 'k:', 'LineWidth', 1.4);
legend(ax2, [hDT, hST, hPre, hRec], ...
    {'Direct total', 'Spectral total', 'Precompute (dashed)', 'Reconstruction (dotted)'}, ...
    'Location', 'northwest', 'NumColumns', 1, ...
    'FontName', 'Times New Roman', 'FontSize', 8, 'Box', 'on', ...
    'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
text(ax2, 0.97, 0.95, '(b)', 'Units', 'normalized', 'FontName', 'Times New Roman', ...
    'FontSize', 11, 'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'right', ...
    'BackgroundColor', 'w', 'Margin', 1);

ts = datestr(now, 'yyyymmdd_HHMMSS');
out_png = ['mf12_theoretical_complexity_memory_pub_' ts '.png'];

exportgraphics(fig, out_png, 'Resolution', 600);

fprintf('Saved figure (PNG, 600 dpi): %s\n', out_png);
fprintf('Assumptions: end-to-end model, fixed grid %dx%d, order=%d.\n', cfg.Nx, cfg.Ny, cfg.order);
fprintf('Pipelines: Direct=coeffsMF12+surfaceMF12_new, Spectral=coeffsMF12_superharmonic+surfaceMF12_spectral.\n');

% Full polynomial expressions from loop counts (not only leading terms)
% Direct precompute: N + 3N(N-1) + 4*C3 = (2/3)N^3 + 1*N^2 - (2/3)N
% Spectral precompute: N + (2N(N-1)+C2) + C3 = (1/6)N^3 + 2N^2 - (7/6)N
% Direct reconstruction terms: 3N + 2N(N-1) + 4*C3 = (2/3)N^3 + (7/3)N
% Spectral reconstruction terms: 3N + 2N(N-1) + C3 = (1/6)N^3 + (3/2)N^2 + (4/3)N
fprintf('Full count model:\n');
fprintf('  M_pre_direct(N)    = (2/3)N^3 + 1*N^2 - (2/3)N\n');
fprintf('  M_pre_spectral(N)  = (1/6)N^3 + 2*N^2 - (7/6)N\n');
fprintf('  M_recon_direct(N)  = (2/3)N^3 + (7/3)N\n');
fprintf('  M_recon_spectral(N)= (1/6)N^3 + (3/2)N^2 + (4/3)N\n');

fprintf('Sanity @ N_max=%d: time ratio direct/spectral = %.3g\n', N(end), ...
    time_direct_units(end) / time_spectral_units(end));
fprintf('Sanity @ N_max=%d: mem direct = %.3f MB, mem spectral = %.3f MB\n', ...
    N(end), mem_direct_mb(end), mem_spectral_mb(end));

% Memory breakdown at N_max
iN = numel(N);
md_pre = mem_pre_direct_bytes(iN) / 1024^2;
md_rec = mem_recon_direct_bytes(iN) / 1024^2;
ms_pre = mem_pre_spectral_bytes(iN) / 1024^2;
ms_rec = mem_recon_spectral_bytes(iN) / 1024^2;
fprintf('Memory breakdown @ N=%d (MB):\n', N(end));
fprintf('  Direct   = precompute %.3f + reconstruction %.3f = %.3f\n', md_pre, md_rec, md_pre + md_rec);
fprintf('  Spectral = precompute %.3f + reconstruction %.3f = %.3f\n', ms_pre, ms_rec, ms_pre + ms_rec);
