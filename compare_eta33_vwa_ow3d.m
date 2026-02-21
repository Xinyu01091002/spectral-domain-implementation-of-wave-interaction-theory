
fid_prog = fopen('progress.log', 'w');
fprintf(fid_prog, 'Start Script\n');
fclose(fid_prog);

% compare_eta33_vwa_ow3d.m
% Compare OW3D and VWA third-order eta slices (centerline/off-centerline).
clc; clear; close all;
addpath(genpath('irregularWavesMF12'));

% Enable MF12 progress/timing output for this run
setenv('MF12_PROGRESS', '1');


fid_prog = fopen('progress.log', 'a'); fprintf(fid_prog, 'Init Done. Loaded Data.\n'); fclose(fid_prog);
config.dataset_label = 'test6';
config.kd = 5.0;
config.spread_deg = 15;
config.Akp = 0.12;
config.time_step = 1200;
lambda = 225;
export_folder = fullfile(pwd, 'figures', 'comparison');
if ~isfolder(export_folder)
    mkdir(export_folder);
end

proc_dir = fullfile(pwd, 'processed_eta33');
if ~isfolder(proc_dir)
    error('processed_eta33 folder missing; run extract_eta33_from_OW3D first.');
end
snapshot_tag = sprintf('OW3D_eta33_%s_kd%.1f_spread%d_Akp%.2f_t%05d.mat', ...
    config.dataset_label, config.kd, config.spread_deg, config.Akp, config.time_step);
snapshot_path = fullfile(proc_dir, snapshot_tag);
if ~isfile(snapshot_path)
    error('Missing snapshot %s', snapshot_tag);
end
S = load(snapshot_path);

X = S.X;
Y = S.Y;
eta11 = S.eta11;
eta33_ow3d = S.eta33; % assume field stored
x_vec = X(:,1);
y_vec = Y(1,:);
x_norm = x_vec / lambda;

g = 9.81;
kp_ref = 2*pi/lambda;

[center_idx, off_idx, center_y, off_y] = locate_center_offline(eta11, y_vec, lambda, 'second');
eta33_center_ow3d = eta33_ow3d(:, center_idx);
eta33_off_ow3d = eta33_ow3d(:, off_idx);
eta11_center_ow3d = eta11(:, center_idx);
eta11_off_ow3d = eta11(:, off_idx);


eta33_vwa = compute_eta33_vwa(eta11, x_vec, y_vec, config, lambda);
eta33_center_vwa = eta33_vwa(:, center_idx);
eta33_off_vwa = eta33_vwa(:, off_idx);

A = config.Akp / kp_ref;
eta_norm_factor = (A^3) * (kp_ref^2);
if ~isfinite(eta_norm_factor) || eta_norm_factor <= 0
    eta_norm_factor = 1;
end

eta33_center_ow3d = eta33_center_ow3d / eta_norm_factor;
eta33_off_ow3d = eta33_off_ow3d / eta_norm_factor;
eta33_center_vwa = eta33_center_vwa / eta_norm_factor;
eta33_off_vwa = eta33_off_vwa / eta_norm_factor;

% =========================================================================
% New Block: Spectral Implementation (MF12) Comparison
% =========================================================================
fprintf('Computing MF12 Spectral Implementation...\n');
t_mf12_total = tic;
[Nx_p, Ny_p] = size(eta11); % Dim 1 is X (Nx), Dim 2 is Y (Ny)
dx = x_vec(2) - x_vec(1);
dy = y_vec(2) - y_vec(1);
Lx = Nx_p * dx; 
Ly = Ny_p * dy;
h = config.kd / kp_ref;
Ux = 0; Uy = 0;

% 1. FFT
% coeff_grid(i,j) corresponds to freq i in dim1 (x) and freq j in dim2 (y)
coeff_grid = fft2(eta11) / (Nx_p * Ny_p);

% 2. Wavenumber Grid (NDGRID to match Matrix Indexing)
dkx = 2*pi/Lx;
dky = 2*pi/Ly;

kx_ind = [0:ceil(Nx_p/2)-1, -floor(Nx_p/2):-1]'; 
ky_ind = [0:ceil(Ny_p/2)-1, -floor(Ny_p/2):-1];
kx_vec = kx_ind * dkx;
ky_vec = ky_ind * dky;

[KX, KY] = ndgrid(kx_vec, ky_vec); 
% NOW KX(i,j) is kx_vec(i)

% Energy filter (keep 99% energy)
energy = abs(coeff_grid).^2;
[sorted_E, sort_idx] = sort(energy(:), 'descend');
cum_E = cumsum(sorted_E);
total_E = sum(energy(:));
cutoff_idx = find(cum_E >= 0.9999 * total_E, 1, 'first'); 

% For linear reconstruction check, we can use many components.
% For nonlinear MF12 calculation, we must be strict.
max_components_nonlinear = 800; % Reduced to 600 for speed test (approx 40-60s)
max_components_linear = 5000;

current_cutoff = min([cutoff_idx, max_components_linear, length(sorted_E)]);
keep_indices = sort_idx(1:current_cutoff);

kx_list = []; ky_list = []; a_list = []; b_list = [];

for k = 1:length(keep_indices)
    idx = keep_indices(k);
    
    kkx = KX(idx);
    kky = KY(idx);
    
    % Half-plane filter
    if kkx < 0 || (kkx == 0 && kky < 0)
        continue;
    end
    
    val = coeff_grid(idx);
    
    if (kkx == 0 && kky == 0)
        amp_a = real(val); amp_b = 0;
    else
        % Reverting to Factor 2 based on test_fft.m and Factor=4 giving 2.0x recovery.
        % If Factor 2 gives 0.5, then there is a mystery.
        amp_a = 2 * real(val);
        amp_b = -2 * imag(val);
    end
    
    kx_list(end+1, 1) = kkx;
    ky_list(end+1, 1) = kky;
    a_list(end+1, 1) = amp_a;
    b_list(end+1, 1) = amp_b;
end

% DEBUG VALUES
fid_dbg = fopen('debug_values.txt', 'w');
fprintf(fid_dbg, 'Mean Abs Val: %.6f\n', mean(abs(coeff_grid(:))));
fprintf(fid_dbg, 'Mean Abs AmpA (Linear): %.6f\n', mean(abs(a_list)));
fprintf(fid_dbg, 'Max Abs AmpA: %.6f\n', max(abs(a_list)));
fprintf(fid_dbg, 'Num Components: %d\n', length(a_list));
fclose(fid_dbg);


fprintf('MF12: Extracted %d components for linear check.\n', length(a_list));
fid_prog = fopen('progress.log', 'a'); fprintf(fid_prog, 'Extract Linear Done.\n'); fclose(fid_prog);

% --- Linear Check (Full Spectrum subset) ---
% Using all extracted components
coeffs1_full = coeffsMF12(1, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy);
[eta_mf12_lin_check, ~] = surfaceMF12_spectral(coeffs1_full, Lx, Ly, Nx_p, Ny_p, 0);

fprintf('LINEAR CHECK (Top %d comp):\n', length(a_list));
max_input = max(abs(eta11(:)));
max_recon = max(abs(eta_mf12_lin_check(:)));
recovery_ratio = max_recon / max_input;
fprintf('  Max Input: %.4f, Max Recon: %.4f, Ratio: %.4f\n', max_input, max_recon, recovery_ratio);

% RMS Error Check (Phase Check)
% Note: surfaceMF12 outputs (Ny, Nx), but eta11 seems to be (Nx, Ny) based on later code requiring transpose.
diff_field = eta11 - eta_mf12_lin_check.'; 
rms_error = sqrt(mean(diff_field(:).^2));
rel_rms_error = rms_error / sqrt(mean(eta11(:).^2));
fprintf('  RMS Error: %.4f\n', rms_error);
fprintf('  Rel RMS Error: %.4f (Should be near 0 if phases are correct)\n', rel_rms_error);

fid = fopen('linear_check_result_sliced.txt', 'w');
fprintf(fid, 'Max Input eta11: %.6f\n', max_input);
fprintf(fid, 'Max Recon MF12: %.6f\n', max_recon);
fprintf(fid, 'Recovery Ratio: %.6f\n', recovery_ratio);
fprintf(fid, 'Rel RMS Error: %.6f\n', rel_rms_error);
fclose(fid);


% --- Reduce for Nonlinear Calculation ---
% Keep only top N for nonlinear (assuming a_list is roughly sorted by energy, but sort_idx was sorted by E)
% Since keep_indices was sorted by E, a_list is also sorted by E order (iterated 1:length).
if length(a_list) > max_components_nonlinear
    fprintf('Reducing to %d components for Nonlinear MF12...\n', max_components_nonlinear);
    a_list = a_list(1:max_components_nonlinear);
    b_list = b_list(1:max_components_nonlinear);
    kx_list = kx_list(1:max_components_nonlinear);
    ky_list = ky_list(1:max_components_nonlinear);
end

fid_prog = fopen('progress.log', 'a'); fprintf(fid_prog, 'Reducing Components Done. Starting 3rd Order Calc.\n'); fclose(fid_prog);

fprintf('MF12: Using %d components for 3rd Order calc.\n', length(a_list));

try
    % Compute 2nd and 3rd order
    % Note: t=0 because we extracted phases from snapshot directly
    % USE SUPERHARMONIC OPTIMIZED VERSION
    fprintf('Using coeffsMF12_superharmonic for speed...\n');
    coeffs3 = coeffsMF12_superharmonic(3, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy);
    fid_prog = fopen('progress.log', 'a'); fprintf(fid_prog, 'Coeffs3 Done.\n'); fclose(fid_prog);
    coeffs2 = coeffsMF12_superharmonic(2, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy);
    fid_prog = fopen('progress.log', 'a'); fprintf(fid_prog, 'Coeffs2 Done.\n'); fclose(fid_prog);

    % Reconstruct Surfaces
    [eta_mf12_3, ~] = surfaceMF12_spectral(coeffs3, Lx, Ly, Nx_p, Ny_p, 0);
    fid_prog = fopen('progress.log', 'a'); fprintf(fid_prog, 'Surface3 Done.\n'); fclose(fid_prog);
    [eta_mf12_2, ~] = surfaceMF12_spectral(coeffs2, Lx, Ly, Nx_p, Ny_p, 0);

    % --- Linear Check (SLICED) ---
    coeffs1 = coeffsMF12(1, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy);
    [eta_mf12_lin_check, ~] = surfaceMF12_spectral(coeffs1, Lx, Ly, Nx_p, Ny_p, 0);

    fprintf('LINEAR CHECK (SLICED - %d comp):\n', length(a_list));
    max_input = max(abs(eta11(:)));
    max_recon = max(abs(eta_mf12_lin_check(:)));
    recovery_ratio = max_recon / max_input;

    fprintf('  Max Input eta11:    %.4f\n', max_input);
    fprintf('  Max Recon MF12:     %.4f\n', max_recon);
    fprintf('  Recovery Ratio:     %.4f\n', recovery_ratio);

    % DEBUG: Check ranges
    eta_mf12_3_min = min(eta_mf12_3(:));
    eta_mf12_3_max = max(eta_mf12_3(:));
    fprintf('DEBUG: eta_mf12_3 range: %.4e to %.4e\n', eta_mf12_3_min, eta_mf12_3_max);
    fprintf('DEBUG: eta_mf12_2 range: %.4e to %.4e\n', min(eta_mf12_2(:)), max(eta_mf12_2(:)));

    fid = fopen('nonlinear_check_result.txt', 'w');
    fprintf(fid, '3rd Order Min: %.6e\n', eta_mf12_3_min);
    fprintf(fid, '3rd Order Max: %.6e\n', eta_mf12_3_max);
    fclose(fid);

    % Transpose to get (Nx_p, Ny_p) match eta11
    eta33_mf12 = (eta_mf12_3 - eta_mf12_2).';
    mf12_elapsed = toc(t_mf12_total);
    fprintf('MF12 Spectral total time: %.2f s\n', mf12_elapsed);
    fid_prog = fopen('progress.log', 'a');
    fprintf(fid_prog, 'MF12 Spectral total time: %.2f s\n', mf12_elapsed);
    fclose(fid_prog);

    % =========================================================================
    % 4. NEW INTEGRATED OPTIMIZED METHOD TEST
    % =========================================================================
    fprintf('\n--- Testing Integrated Optimized Function (Low Memory) ---\n');
    t_opt = tic;
    
    % Call the new function (Order 3)
    % Note: It handles 3rd order loops internally.
    [eta_opt_3, phi_opt_3] = surfaceMF12_integrated_opt(3, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy, Lx, Ly, Nx_p, Ny_p, 0);
    
    % For consistent comparison we need 2nd order too (Wait, the function computes total surface up to Order X)
    % surfaceMF12_integrated_opt returns the TOTAL accumulated surface (Linear + 2nd + 3rd...) depending on 'order' input?
    % Let's check the function logic.
    % Yes: "Accumulate Linear", then "if order>=2 Accumulate 2nd", "if order>=3 Accumulate 3rd".
    % So eta_opt_3 is the full eta_3.
    % To get eta33, we need eta_3 - eta_2.
    
    [eta_opt_2, ~] = surfaceMF12_integrated_opt(2, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy, Lx, Ly, Nx_p, Ny_p, 0);
    
    t_opt_elapsed = toc(t_opt);
    fprintf('Integrated Opt Total Time: %.2f s\n', t_opt_elapsed);
    
    % Validation
    eta33_opt = (eta_opt_3 - eta_opt_2).';
    diff_opt = eta33_opt - eta33_mf12;
    err_opt = max(abs(diff_opt(:)));
    fprintf('Difference between Standard and Integrated: %.6e\n', err_opt);
    
    if err_opt < 1e-4
        fprintf('>>> SUCCESS: Integrated method matches standard method.\n');
        % Use the optimized result for plotting if desired
        % eta33_mf12 = eta33_opt; 
    else
        fprintf('>>> WARNING: Match failed. Check implementation.\n');
    end
    % =========================================================================

catch e
    fprintf('ERROR in Nonlinear Calc: %s\n', e.message);
    fprintf('Stack: %s line %d\n', e.stack(1).name, e.stack(1).line);
    mf12_elapsed = toc(t_mf12_total);
    fprintf('MF12 Spectral time (failed): %.2f s\n', mf12_elapsed);
    fid_prog = fopen('progress.log', 'a'); 
    fprintf(fid_prog, 'ERROR: %s\n', e.message); 
    fprintf(fid_prog, 'Line: %d\n', e.stack(1).line);
    fprintf(fid_prog, 'MF12 Spectral time (failed): %.2f s\n', mf12_elapsed);
    fclose(fid_prog);
end


fid = fopen('nonlinear_check_result.txt', 'w');
fprintf(fid, '3rd Order Min: %.6e\n', eta_mf12_3_min);
fprintf(fid, '3rd Order Max: %.6e\n', eta_mf12_3_max);
fclose(fid);


% Transpose to get (Nx_p, Ny_p) match eta11
eta33_mf12 = (eta_mf12_3 - eta_mf12_2).';
% eta33_mf12=eta33_mf12';
% Slicing
eta33_center_mf12 = eta33_mf12(:, center_idx);
eta33_off_mf12 = eta33_mf12(:, off_idx);

% Normalization & Filtering
eta33_center_mf12 = eta33_center_mf12 / eta_norm_factor;
eta33_off_mf12 = eta33_off_mf12 / eta_norm_factor;

fprintf('DEBUG: eta33_center_mf12 range: %.4e to %.4e\n', min(eta33_center_mf12), max(eta33_center_mf12));
fprintf('\n=== MAGNITUDE CHECK ===\n');
fprintf('Centerline Ranges (Normalized):\n');
fprintf('  OW3D (Blue): %.4e to %.4e\n', min(eta33_center_ow3d), max(eta33_center_ow3d));
fprintf('  VWA  (Red):  %.4e to %.4e\n', min(eta33_center_vwa), max(eta33_center_vwa));
fprintf('  MF12 (Grn):  %.4e to %.4e\n', min(eta33_center_mf12), max(eta33_center_mf12));
fprintf('=======================\n\n');

% eta33_center_mf12 = frequency_filtering(eta33_center_mf12.', x_vec, kp_ref, 3).';
% eta33_off_mf12 = frequency_filtering(eta33_off_mf12.', x_vec, kp_ref, 3).';

mf12Color = [0.4660 0.6740 0.1880]; % Greenish

crest_width_lambda = 3;
center_env = abs(hilbert(eta11_center_ow3d));
off_env = abs(hilbert(eta11_off_ow3d));
[~, center_peak_idx] = max(center_env);
[~, off_peak_idx] = max(off_env);
center_xlim = [x_norm(1), x_norm(end)];
off_xlim = center_xlim;
if ~isempty(center_peak_idx)
    px = x_norm(center_peak_idx);
    center_xlim = [max(x_norm(1), px - crest_width_lambda), min(x_norm(end), px + crest_width_lambda)];
end
if ~isempty(off_peak_idx)
    px = x_norm(off_peak_idx);
    off_xlim = [max(x_norm(1), px - crest_width_lambda), min(x_norm(end), px + crest_width_lambda)];
end

pub_line_width = 3.2;
pub_font_size = 26;
axis_line_width = 1.8;
label_font_size = pub_font_size + 6;
annotation_font_size = pub_font_size - 2;
exactColor = [0.1176 0.2980 0.5765];
vwaColor = [0.7961 0.1882 0.2157];
eta33_center_ow3d=frequency_filtering(eta33_center_ow3d,x_vec, kp_ref,3);
eta33_off_ow3d=frequency_filtering(eta33_off_ow3d,x_vec, kp_ref,3);
fig_center = figure('Color', 'w', 'Position', [120, 120, 900, 900]);
set(fig_center, 'PaperUnits', 'inches', 'PaperSize', [9 9], 'PaperPositionMode', 'auto');
ax_center = axes(fig_center);
plot(ax_center, x_norm, eta33_center_ow3d, 'Color', exactColor, 'LineWidth', pub_line_width, 'DisplayName', 'FNPF');
hold(ax_center, 'on');
plot(ax_center, x_norm, eta33_center_vwa, '--', 'Color', vwaColor, 'LineWidth', pub_line_width, 'DisplayName', 'VWA');
plot(ax_center, x_norm, eta33_center_mf12, '-', 'Color', mf12Color, 'LineWidth', pub_line_width, 'DisplayName', 'MF12 (Spectral)');
grid(ax_center, 'on');
grid(ax_center, 'minor');
set(ax_center, 'GridAlpha', 0.3, 'GridLineStyle', '-', 'MinorGridLineStyle', ':', 'MinorGridAlpha', 0.2);
xlabel(ax_center, '$x/\lambda$', 'Interpreter', 'latex', 'FontSize', label_font_size);
ylabel(ax_center, '$\eta^{(33)}/(A^3 k_p^2)$', 'Interpreter', 'latex', 'FontSize', label_font_size);
legend(ax_center, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', pub_font_size+4, 'Box', 'off');
set(ax_center, 'FontSize', pub_font_size, 'LineWidth', axis_line_width, 'TickLabelInterpreter', 'latex');
ytickformat(ax_center, '%.2f');
ax_center.YAxis.Exponent = computeAxisExponent(ax_center, 'Y');
title(ax_center, sprintf('Centerline ($y/\\lambda=%.3f$)', center_y / lambda-20), 'Interpreter', 'latex', 'FontSize', annotation_font_size);
text(ax_center, -0.05, 0.97, '(a)', 'Units', 'normalized', 'FontSize', label_font_size, 'FontWeight', 'bold', 'Interpreter', 'latex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
box(ax_center, 'on');
set(ax_center, 'LooseInset', [0 0 0 0]);
if center_xlim(2) > center_xlim(1)
    xlim(ax_center, center_xlim);
end
center_filename = sprintf('eta33_centerline_comparison_spread%d_kd%.1f.png', config.spread_deg, config.kd);
exportgraphics(fig_center, fullfile(export_folder, center_filename), 'Resolution', 300);

fig_off = figure('Color', 'w', 'Position', [120, 120, 900, 900]);
set(fig_off, 'PaperUnits', 'inches', 'PaperSize', [9 9], 'PaperPositionMode', 'auto');
ax_off = axes(fig_off);
plot(ax_off, x_norm, eta33_off_ow3d, 'Color', exactColor, 'LineWidth', pub_line_width, 'DisplayName', 'FNPF');
hold(ax_off, 'on');
plot(ax_off, x_norm, eta33_off_vwa, '--', 'Color', vwaColor, 'LineWidth', pub_line_width, 'DisplayName', 'VWA');
plot(ax_off, x_norm, eta33_off_mf12, '-', 'Color', mf12Color, 'LineWidth', pub_line_width, 'DisplayName', 'MF12 (Spectral)');
grid(ax_off, 'on');
grid(ax_off, 'minor');
set(ax_off, 'GridAlpha', 0.3, 'GridLineStyle', '-', 'MinorGridLineStyle', ':', 'MinorGridAlpha', 0.2);
xlabel(ax_off, '$x/\lambda$', 'Interpreter', 'latex', 'FontSize', label_font_size);
ylabel(ax_off, '$\eta^{(33)}/(A^3 k_p^2)$', 'Interpreter', 'latex', 'FontSize', label_font_size);
legend(ax_off, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', pub_font_size+4, 'Box', 'off');
set(ax_off, 'FontSize', pub_font_size, 'LineWidth', axis_line_width, 'TickLabelInterpreter', 'latex');
ytickformat(ax_off, '%.2f');
ax_off.YAxis.Exponent = computeAxisExponent(ax_off, 'Y');
title(ax_off, sprintf('Off-centerline ($y/\\lambda=%.3f$)', off_y / lambda-20), 'Interpreter', 'latex', 'FontSize', annotation_font_size);
text(ax_off, -0.05, 0.97, '(b)', 'Units', 'normalized', 'FontSize', label_font_size, 'FontWeight', 'bold', 'Interpreter', 'latex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
box(ax_off, 'on');
set(ax_off, 'LooseInset', [0 0 0 0]);
if off_xlim(2) > off_xlim(1)
    xlim(ax_off, off_xlim);
end
offline_filename = sprintf('eta33_offline_comparison_spread%d_kd%.1f.png', config.spread_deg, config.kd);
exportgraphics(fig_off, fullfile(export_folder, offline_filename), 'Resolution', 300);

function eta33_vwa = compute_eta33_vwa(eta11, x_vec, y_vec, config, lambda)
    eta11_field = eta11.';
    [~, ~, KX, ~] = compute_kxky(x_vec, y_vec);

    g = 9.81;
    kp_ref = 2 * pi / lambda;
    h = config.kd / kp_ref;

    sigma = tanh(abs(KX) * h);
    sigma = max(sigma, 1e-6);
    B33 = (27 - 9 * sigma.^2 + 9 * sigma.^4 - 3 * sigma.^6) ./ (64 * sigma.^6);
    B33(~isfinite(B33)) = 0;
    B33(abs(KX).*h<0.3) = 0;

    Kmag = abs(KX);
    Omega = sqrt(g .* Kmag .* tanh(Kmag * h));
    Omega = max(Omega, 1e-6);
    cp = g ./ Omega;
    cp(~isfinite(cp)) = 0;

    Eta_fft = fft2(eta11_field);
    shiftKX = fftshift(KX);
    shiftB33 = fftshift(B33);
    shiftCp = fftshift(cp);

    kappa_vwa_field = ifft2(Eta_fft .* shiftKX .* shiftB33.* shiftKX);
    k_vwa_field = ifft2(Eta_fft .* shiftKX );
    % eta33_vwa_field = real(analytic2D(eta11_field) .* analytic2D(kappa_vwa_field) .* analytic2D(k_vwa_field));
    eta33_vwa_field = real(analytic2D(eta11_field) .* analytic2D(kappa_vwa_field) .* analytic2D(eta11_field));
    eta33_vwa = eta33_vwa_field.';
end

function [center_idx, off_idx, center_y, off_y] = locate_center_offline(eta11, y_values, lambda, orientation)
    if nargin < 4
        orientation = 'second';
    end
    switch lower(orientation)
        case 'second'
            Ny = size(eta11, 2);
            envelope = abs(hilbert(eta11));
            max_env = max(envelope, [], 1);
            [~, center_idx] = max(max_env);
            center_idx = max(1, min(Ny, center_idx));
            center_y = y_values(center_idx);
        otherwise
            error('Unsupported orientation: %s', orientation);
    end
    target_amp = 0.5 * max_env(center_idx);
    absdiff = abs(max_env - target_amp);
    minDiff = min(absdiff);
    tol = max(1e-12, minDiff);
    candidates = find(absdiff <= minDiff + tol);
    candidates(candidates == center_idx) = [];
    if isempty(candidates)
        off_idx = center_idx;
    else
        [~, best] = max(abs(candidates - center_idx));
        off_idx = candidates(best);
    end
    off_idx = max(1, min(off_idx, Ny));
    off_y = y_values(off_idx);
end

function exponent = computeAxisExponent(ax, axisName)
    if nargin < 2
        axisName = 'Y';
    end
    limField = sprintf('%sLim', axisName);
    if ~isprop(ax, limField)
        exponent = 0;
        return;
    end
    axisLimits = ax.(limField);
    maxVal = max(abs(axisLimits));
    if maxVal == 0
        exponent = 0;
        return;
    end
    exponent = floor(log10(maxVal));
end
