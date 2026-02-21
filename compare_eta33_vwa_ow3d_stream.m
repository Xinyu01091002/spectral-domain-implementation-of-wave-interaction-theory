% compare_eta33_vwa_ow3d_stream.m
% Compare OW3D, VWA, and MF12 (Integrated Stream Version) third-order eta.
% Optimized for high component count without memory overflow.

clc; clear; close all;
addpath(genpath('irregularWavesMF12'));

% --- Configuration ---
config.dataset_label = 'test6';
config.kd = 5.0;
config.spread_deg = 15;
config.Akp = 0.12;
config.time_step = 800;
lambda = 225;

export_folder = fullfile(pwd, 'figures', 'comparison_stream');
if ~isfolder(export_folder)
    mkdir(export_folder);
end

% Check dependencies
if ~exist('compute_eta33_vwa', 'file') && ~exist('compare_eta33_vwa_ow3d.m', 'file')
    warning('Reference scripts not found suitable for VWA function copying. Ensuring function exists below.');
end

% --- Load OW3D Data ---
proc_dir = fullfile(pwd, 'processed_eta33');
snapshot_tag = sprintf('OW3D_eta33_%s_kd%.1f_spread%d_Akp%.2f_t%05d.mat', ...
    config.dataset_label, config.kd, config.spread_deg, config.Akp, config.time_step);
snapshot_path = fullfile(proc_dir, snapshot_tag);

if ~isfile(snapshot_path)
    error('Missing snapshot: %s', snapshot_path);
end
fprintf('Loading snapshot: %s\n', snapshot_path);
S = load(snapshot_path);

X = S.X; Y = S.Y;
eta11 = S.eta11;
eta33_ow3d = S.eta33;
x_vec = X(:,1);
y_vec = Y(1,:);
x_norm = x_vec / lambda;

[Nx_p, Ny_p] = size(eta11);
Lx = Nx_p * (x_vec(2)-x_vec(1));
Ly = Ny_p * (y_vec(2)-y_vec(1));

% Physics Params
g = 9.81;
kp_ref = 2*pi/lambda;
h = config.kd / kp_ref;
Ux = 0; Uy = 0;

% Normalization Factor
A = config.Akp / kp_ref;
eta_norm_factor = (A^3) * (kp_ref^2);
if ~isfinite(eta_norm_factor) || eta_norm_factor <= 0, eta_norm_factor = 1; end

% Locate Centers
[center_idx, off_idx, center_y, off_y] = locate_center_offline(eta11, y_vec, lambda, 'second');

% --- 1. Compute VWA ---
fprintf('Computing VWA...\n');
eta33_vwa = compute_eta33_vwa(eta11, x_vec, y_vec, config, lambda);

% --- 2. MF12 Initialization (Extract Components) ---
fprintf('Extracting Wave Components (FFT)...\n');
coeff_grid = fft2(eta11) / (Nx_p * Ny_p);
dkx = 2*pi/Lx; dky = 2*pi/Ly;
kx_ind = [0:ceil(Nx_p/2)-1, -floor(Nx_p/2):-1]'; 
ky_ind = [0:ceil(Ny_p/2)-1, -floor(Ny_p/2):-1];
kx_vec_grid = kx_ind * dkx;
ky_vec_grid = ky_ind * dky;
[KX, KY] = ndgrid(kx_vec_grid, ky_vec_grid); 

energy = abs(coeff_grid).^2;
[sorted_E, sort_idx] = sort(energy(:), 'descend');
cum_E = cumsum(sorted_E);
total_E = sum(energy(:));
cutoff_idx = find(cum_E >= 0.999 * total_E, 1, 'first'); 

% Increase component count for precision since we use the optimized method
% Standard memory limit is ~200-800. Stream mode can handle 3000-5000 easily.
max_components = 1800; 
current_cutoff = min([cutoff_idx, max_components, length(sorted_E)]);
keep_indices = sort_idx(1:current_cutoff);

kx_list = []; ky_list = []; a_list = []; b_list = [];
for k = 1:length(keep_indices)
    idx = keep_indices(k);
    kkx = KX(idx); kky = KY(idx);
    if kkx < 0 || (kkx == 0 && kky < 0), continue; end
    
    val = coeff_grid(idx);
    if (kkx == 0 && kky == 0)
        amp_a = real(val); amp_b = 0;
    else
        amp_a = 2 * real(val); amp_b = -2 * imag(val);
    end
    
    kx_list(end+1, 1) = kkx;
    ky_list(end+1, 1) = kky;
    a_list(end+1, 1) = amp_a;
    b_list(end+1, 1) = amp_b;
end
fprintf('Using %d components for MF12 Stream Calculation.\n', length(a_list));

% --- 3. MF12 Integrated Stream Calculation ---
fprintf('Computing MF12 (Integrated Stream Version)...\n');
try
    t_opt = tic;
    % Calculate Full Surface up to Order 3
    fprintf('  >> Pass 1: Order 3 (Full)...\n');
    [eta_opt_3, ~] = surfaceMF12_integrated_opt(3, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy, Lx, Ly, Nx_p, Ny_p, 0);
    
    % Calculate Full Surface up to Order 2
    fprintf('  >> Pass 2: Order 2 (Subtraction)...\n');
    [eta_opt_2, ~] = surfaceMF12_integrated_opt(2, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy, Lx, Ly, Nx_p, Ny_p, 0);
    
    % Difference is purely 3rd Order (and technically higher order terms if they were incl, but logic filters)
    % Note: surfaceMF12_integrated_opt accumulates L + Q (+ C). 
    eta33_mf12_field = (eta_opt_3 - eta_opt_2);
    eta33_mf12 = eta33_mf12_field.'; % Transpose to match (Nx, Ny) of script conventions
    
    t_total = toc(t_opt);
    fprintf('MF12 Stream Calc Complete in %.2f s\n', t_total);
    
catch ME
    fprintf('Error in MF12 Stream Calc: %s\n', ME.message);
    rethrow(ME);
end

% --- 4. Slicing & Comparison ---
eta33_center_ow3d = eta33_ow3d(:, center_idx) / eta_norm_factor;
eta33_off_ow3d = eta33_ow3d(:, off_idx) / eta_norm_factor;

eta33_center_vwa = eta33_vwa(:, center_idx) / eta_norm_factor;
eta33_off_vwa = eta33_vwa(:, off_idx) / eta_norm_factor;

eta33_center_mf12 = eta33_mf12(:, center_idx) / eta_norm_factor;
eta33_off_mf12 = eta33_mf12(:, off_idx) / eta_norm_factor;

% Filtering (Optional, based on original script, can be enabled/disabled)
% eta33_center_ow3d = frequency_filtering(eta33_center_ow3d, x_vec, kp_ref, 3);
% eta33_center_mf12 = frequency_filtering(eta33_center_mf12, x_vec, kp_ref, 3);

fprintf('\n=== RESULTS SUMMARY ===\n');
fprintf('Max Amplitudes (Normalized):\n');
fprintf('  OW3D: %.4f\n', max(abs(eta33_center_ow3d)));
fprintf('  VWA:  %.4f\n', max(abs(eta33_center_vwa)));
fprintf('  MF12: %.4f\n', max(abs(eta33_center_mf12)));
fprintf('=======================\n');

% --- 5. Plotting ---
mf12Color = [0.4660 0.6740 0.1880]; % Green
exactColor = [0.1176 0.2980 0.5765]; % Blue
vwaColor = [0.7961 0.1882 0.2157]; % Red
pub_line_width = 2.0;
pub_font_size = 14;

% Auto-focus on Linear Peak (to handle t=1200 or shifting waves)
eta11_center = eta11(:, center_idx);
eta11_off = eta11(:, off_idx);
peak_width_view = 3; % View +/- 3 wavelengths around peak

% Centerline Envelope Peak
env_cen = abs(hilbert(eta11_center));
[~, p_idx_cen] = max(env_cen);
px_cen = x_norm(p_idx_cen);
xlim_cen = [max(x_norm(1), px_cen - peak_width_view), min(x_norm(end), px_cen + peak_width_view)];

% Off-line Envelope Peak
env_off = abs(hilbert(eta11_off));
[~, p_idx_off] = max(env_off);
px_off = x_norm(p_idx_off);
xlim_off = [max(x_norm(1), px_off - peak_width_view), min(x_norm(end), px_off + peak_width_view)];

% Figure 1: Centerline
f1 = figure('Color','w', 'Position', [100 100 800 500]);
ax1 = axes(f1);
plot(x_norm, eta33_center_ow3d, 'Color', exactColor, 'LineWidth', pub_line_width, 'DisplayName', 'OW3D (FNPF)'); hold on;
plot(x_norm, eta33_center_vwa, '--', 'Color', vwaColor, 'LineWidth', pub_line_width, 'DisplayName', 'VWA');
plot(x_norm, eta33_center_mf12, '-.', 'Color', mf12Color, 'LineWidth', pub_line_width, 'DisplayName', 'MF12 (Stream)');
grid on; xlabel('$x/\lambda$', 'Interpreter','latex'); ylabel('$\eta^{(33)}/(A^3 k_p^2)$', 'Interpreter','latex');
title(sprintf('Centerline Comparison ($y/\\lambda=%.3f$)', center_y/lambda), 'Interpreter','latex');
legend('Location', 'best');
xlim(xlim_cen); 

% Figure 2: Off-Centerline
f2 = figure('Color','w', 'Position', [150 150 800 500]);
ax2 = axes(f2);
plot(x_norm, eta33_off_ow3d, 'Color', exactColor, 'LineWidth', pub_line_width, 'DisplayName', 'OW3D (FNPF)'); hold on;
plot(x_norm, eta33_off_vwa, '--', 'Color', vwaColor, 'LineWidth', pub_line_width, 'DisplayName', 'VWA');
plot(x_norm, eta33_off_mf12, '-.', 'Color', mf12Color, 'LineWidth', pub_line_width, 'DisplayName', 'MF12 (Stream)');
grid on; xlabel('$x/\lambda$', 'Interpreter','latex'); ylabel('$\eta^{(33)}/(A^3 k_p^2)$', 'Interpreter','latex');
title(sprintf('Off-Centerline Comparison ($y/\\lambda=%.3f$)', off_y/lambda), 'Interpreter','latex');
legend('Location', 'best');
xlim(xlim_off);

% Export
exportgraphics(f1, fullfile(export_folder, 'eta33_stream_centerline.png'), 'Resolution', 300);
exportgraphics(f2, fullfile(export_folder, 'eta33_stream_offline.png'), 'Resolution', 300);
fprintf('Figures saved to: %s\n', export_folder);


% --- HELPER FUNCTIONS ---

function eta33_vwa = compute_eta33_vwa(eta11, x_vec, y_vec, config, lambda)
    % VWA Implementation
    eta11_field = eta11.';
    [~, ~, KX, ~] = compute_kxky(x_vec, y_vec);

    g = 9.81;
    kp_ref = 2 * pi / lambda;
    h = config.kd / kp_ref;

    sigma = tanh(abs(KX) * h);
    sigma = max(sigma, 1e-6);
    % VWA B33 Kernel
    B33 = (27 - 9 * sigma.^2 + 9 * sigma.^4 - 3 * sigma.^6) ./ (64 * sigma.^6);
    B33(~isfinite(B33)) = 0;
    B33(abs(KX).*h<0.3) = 0; % Shallow cut

    Kmag = abs(KX);
    Omega = sqrt(g .* Kmag .* tanh(Kmag * h));
    Omega = max(Omega, 1e-6);
    cp = g ./ Omega;
    cp(~isfinite(cp)) = 0;

    Eta_fft = fft2(eta11_field);
    shiftKX = fftshift(KX);
    shiftB33 = fftshift(B33);

    kappa_vwa_field = ifft2(Eta_fft .* shiftKX .* shiftB33.* shiftKX);
    % Note: Simplified VWA term approximation used in previous script
    % explicit call to file on path is expected
    eta33_vwa_field = real(analytic2D(eta11_field) .* analytic2D(kappa_vwa_field) .* analytic2D(eta11_field));
    eta33_vwa = eta33_vwa_field.';
end

function [center_idx, off_idx, center_y, off_y] = locate_center_offline(eta11, y_values, lambda, orientation)
    % Locate peak y
    Ny = size(eta11, 2);
    envelope = abs(hilbert(eta11));
    max_env = max(envelope, [], 1);
    [~, center_idx] = max(max_env);
    center_idx = max(1, min(Ny, center_idx));
    center_y = y_values(center_idx);
    
    % Locate ~50% amplitude offline
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

