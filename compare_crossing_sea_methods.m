clc;clear ;close all;
addpath(genpath(fullfile(pwd, 'irregularWavesMF12')));
rehash toolboxcache;


% compare_crossing_sea_methods.m
% Compare Original MF12 3rd Order (Spectral FFT) vs Flow+Spectral Integrated Opt 3rd Order
% for a Crossing Sea state.

try
    fprintf('--- Starting Comparison Script ---\n');
    diary('compare_log.txt');
    diary on;
    % 1. Basic Physics & Grid Setup
g = 9.81;
h = 1500; 
Lx = 3000; 
Ly = 3000; 
% Reduced grid for feasible 3rd order full coefficient calculation
Nx = 128;  
Ny = 32;   

dx = Lx/Nx;
dy = Ly/Ny;
x_vec = (0:Nx-1)*dx;
y_vec = (0:Ny-1)*dy;
[X_mesh, Y_mesh] = meshgrid(x_vec, y_vec);

dkx = 2*pi/Lx;
dky = 2*pi/Ly;

% 2. Wave Parameters
Tp = 12; 
kp = (2*pi/Tp)^2/g; 
lambda = 2*pi/kp; % Wavelength
Akp = 2; 
alpha = 8; 
kw_left = 0.004606; 
kw_right = kw_left * 1;

A_focus_total = Akp / kp; 

% Crossing Sea
theta_1_deg = 0; spread_1_deg = 5; weight_1 = 0.5;
theta_2_deg = -45; spread_2_deg = 5; weight_2 = 0.5;

tf = 0.0; 
xf = Lx/2; 
yf = Ly/2;

fprintf('Generating Wave Components...\n');
% 3. Spectral Generation
[Idx_x, Idx_y] = meshgrid(-Nx/2 : Nx/2-1, -Ny/2 : Ny/2-1);
KX_grid = Idx_x * dkx;
KY_grid = Idx_y * dky;

K_mag = sqrt(KX_grid.^2 + KY_grid.^2);
valid_mask = K_mag > 0;
KX = KX_grid(valid_mask);
KY = KY_grid(valid_mask);
K_mag = K_mag(valid_mask);
Theta = atan2(KY, KX); 

kw_vec = ones(size(K_mag)) * kw_right;
kw_vec(K_mag <= kp) = kw_left;
S_vec = exp(-(K_mag - kp).^2 ./ (2*kw_vec.^2));

theta_1 = deg2rad(theta_1_deg);
theta_2 = deg2rad(theta_2_deg);
w_1 = deg2rad(spread_1_deg);
w_2 = deg2rad(spread_2_deg);

ang_diff = @(t, t_ref) angle(exp(1i*(t - t_ref)));
D_vec_1 = exp(-(ang_diff(Theta, theta_1)).^2 / (2*w_1^2));
D_vec_2 = exp(-(ang_diff(Theta, theta_2)).^2 / (2*w_2^2));

Amp_raw = S_vec .* (weight_1 * D_vec_1 + weight_2 * D_vec_2);

% Energy Cutoff & Limit
Energy_raw = Amp_raw.^2;
[sorted_E, sort_idx] = sort(Energy_raw, 'descend');
cum_E = cumsum(sorted_E);
total_E = sum(Energy_raw);

cutoff_ratio = 0.95; 
last_idx = find(cum_E >= cutoff_ratio * total_E, 1, 'first');

% LIMIT COMPONENT COUNT to prevent O(N^3) explosion in memory/time for the naive method
max_waves = 100; 
if last_idx > max_waves
    fprintf('  Warning: Limiting components from %d to %d for performance.\n', last_idx, max_waves);
    last_idx = max_waves;
end

keep_indices = sort_idx(1:last_idx);
temp_kx = KX(keep_indices);
temp_ky = KY(keep_indices);
temp_amp = Amp_raw(keep_indices);

% Scaling
if isempty(temp_amp)
    error('No components generated.');
end
scale_factor = A_focus_total / sum(temp_amp);
temp_amp = temp_amp * scale_factor;

% Phase (Focusing)
kx_list = []; ky_list = []; a_list = []; b_list = [];
for i = 1:length(temp_kx)
    kk_x = temp_kx(i);
    kk_y = temp_ky(i);
    kk = sqrt(kk_x^2 + kk_y^2);
    ww = sqrt(g * kk * tanh(kk * h)); 
    
    phase = -kk_x * xf - kk_y * yf + ww * tf; 
    
    kx_list(end+1,1) = kk_x;
    ky_list(end+1,1) = kk_y;
    a_list(end+1,1) = temp_amp(i) * cos(phase);
    b_list(end+1,1) = temp_amp(i) * sin(phase);
end

Ux = 0; Uy = 0;
N_comp = length(a_list);
fprintf('  Final N = %d components.\n', N_comp);

% =========================================================================
% METHOD 1: Original MF12 (Spectral / Full Coeff Matrix)
% =========================================================================
fprintf('\n--- Method 1: Original MF12 (Full Coefficients) ---\n');
tic;
% 1st Order
fprintf('  Computing 1st Order Coeffs...\n');
coeffs1 = coeffsMF12(1, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy);
[eta_mf12_1,~] = surfaceMF12_new(1, coeffs1, X_mesh, Y_mesh, tf);

% 2nd Order
fprintf('  Computing 2nd Order Coeffs...\n');
coeffs2 = coeffsMF12(2, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy);
fprintf('  Reconstructing 2nd Order Surface (surfaceMF12_new)...\n');
[eta_mf12_2,~] = surfaceMF12_new(2, coeffs2, X_mesh, Y_mesh, tf);

% 3rd Order
fprintf('  Computing 3rd Order Coeffs (O(N^3))...\n');
coeffs3 = coeffsMF12(3, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy);
fprintf('  Reconstructing 3rd Order Surface (surfaceMF12_new)...\n');
[eta_mf12_3,~] = surfaceMF12_new(3, coeffs3, X_mesh, Y_mesh, tf);

% Extract Pure Terms
eta_mf12_3rd_term = eta_mf12_3 - eta_mf12_2;
eta_mf12_2nd_term = eta_mf12_2 - eta_mf12_1;

t_mf12 = toc;
fprintf('  >> Time: %.2f s\n', t_mf12);

% =========================================================================
% METHOD 2: Flow + Spectral (Integrated Optimization)
% =========================================================================
fprintf('\n--- Method 2: Flow+Spectral (Integrated Opt) ---\n');
tic;
% Note: surfaceMF12_integrated_opt returns magnitude accumulated, so 3rd call includes 2nd.
fprintf('  Computing Order 3 Integrated...\n');
[eta_stream_3, ~] = surfaceMF12_integrated_opt(3, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy, Lx, Ly, Nx, Ny, tf);

fprintf('  Computing Order 2 Integrated...\n');
[eta_stream_2, ~] = surfaceMF12_integrated_opt(2, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy, Lx, Ly, Nx, Ny, tf);

fprintf('  Computing Order 1 Integrated (Linear)...\n');
[eta_stream_1, ~] = surfaceMF12_integrated_opt(1, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy, Lx, Ly, Nx, Ny, tf);

% Extract Pure 3rd Order & 2nd Order
eta_stream_3rd_term = eta_stream_3 - eta_stream_2;
eta_stream_2nd_term = eta_stream_2 - eta_stream_1;
t_stream = toc;
fprintf('  >> Time: %.2f s\n', t_stream);

% =========================================================================
% METHOD 3: Modified Superharmonic With Subharmonics (Verification)
% =========================================================================
fprintf('\n--- Method 3: Modified Superharmonic (With Subharmonics) ---\n');
% Uses coeffsMF12_superharmonic (which we modified to include 2nd order subharmonics)
% + surfaceMF12_spectral

% 1st Order
coeffs_sup_1 = coeffsMF12_superharmonic(1, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy);
[eta_sup_1, ~] = surfaceMF12_spectral(coeffs_sup_1, Lx, Ly, Nx, Ny, tf);

% 2nd Order
coeffs_sup_2 = coeffsMF12_superharmonic(2, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy);
[eta_sup_2, ~] = surfaceMF12_spectral(coeffs_sup_2, Lx, Ly, Nx, Ny, tf);
eta_sup_2nd_term = eta_sup_2 - eta_sup_1; % Remove linear

% 3rd Order (Modified: Linear + 2nd(Sum+Diff) + 3rd(Sum Only))
coeffs_sup_3 = coeffsMF12_superharmonic(3, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy);
[eta_sup_3, ~] = surfaceMF12_spectral(coeffs_sup_3, Lx, Ly, Nx, Ny, tf);

% Compare 1st Order (Linear): Full (Method 2) vs Modified (Method 3)
% Note: eta_stream_1 is Method 2's linear result.
diff_1st = eta_stream_1 - eta_sup_1;
max_diff_1st = max(abs(diff_1st(:)));
fprintf('  Max Diff (Full vs Modified 1st Order): %.6e\n', max_diff_1st);

% Compare 2nd Order: Full (Method 2 Relative to Method 1) & Method 2 vs Method 3
% Note: eta_stream_2nd_term is Method 2's pure 2nd order.
diff_2nd_sup = eta_stream_2nd_term - eta_sup_2nd_term;
max_diff_2nd_sup = max(abs(diff_2nd_sup(:)));
fprintf('  Max Diff (Full vs Modified 2nd Order): %.6e\n', max_diff_2nd_sup);


% Compare Total Nonlinear Surface: Full (Method 2) vs Modified (Method 3)
diff_total_sup = eta_stream_3 - eta_sup_3;
max_diff_total_sup = max(abs(diff_total_sup(:)));
fprintf('  Max Diff (Full vs Modified Total Surface): %.6e\n', max_diff_total_sup);
fprintf('  (Note: Diff is expected due to missing 3rd order subharmonics in Method 3)\n');



% =========================================================================
% COMPARISON
% =========================================================================
diff = eta_mf12_3rd_term - eta_stream_3rd_term;
max_diff = max(abs(diff(:)));
max_val = max(abs(eta_mf12_3rd_term(:)));
rel_error = max_diff / max_val * 100;

fprintf('\n=== Comparison Result ===\n');
fprintf('  Max Approx Value (MF12): %.6f m\n', max_val);
fprintf('  Max Diff: %.6e m\n', max_diff);
fprintf('  Rel Error: %.4f%%\n', rel_error);


% Compare Total 3rd Order Surface as well
max_val_total = max(abs(eta_mf12_3(:)));
diff_total = eta_mf12_3 - eta_stream_3;
max_diff_total = max(abs(diff_total(:)));
rel_error_total = max_diff_total / max_val_total * 100;

fprintf('\n=== Nonlinear Surface Comparison ===\n');
fprintf('  Max Val: %.6f\n', max_val_total);
fprintf('  Max Diff: %.6e\n', max_diff_total);
fprintf('  Rel Error: %.4f%%\n', rel_error_total);

% Plot
try
    figure('Position',[50 100 1500 800], 'Visible', 'off');

    subplot(2,3,1);
    surf(X_mesh, Y_mesh, eta_mf12_3, 'EdgeColor', 'none');
    view(2); colormap jet; colorbar; shading interp;
    title('MF12 Original (Nonlinear Surface)');
    xlabel('x'); ylabel('y'); axis equal; xlim([0 Lx]); ylim([0 Ly]);
    c_lim_tot = clim;

    subplot(2,3,2);
    surf(X_mesh, Y_mesh, eta_stream_3, 'EdgeColor', 'none');
    view(2); colormap jet; colorbar; shading interp;
    title('Flow+Spectral (Nonlinear Surface)');
    xlabel('x'); ylabel('y'); axis equal; xlim([0 Lx]); ylim([0 Ly]);
    clim(c_lim_tot);

    subplot(2,3,3);
    surf(X_mesh, Y_mesh, abs(diff_total), 'EdgeColor', 'none');
    view(2); colormap jet; colorbar; shading interp;
    title(sprintf('Diff Total (Max: %.2e)', max_diff_total));
    xlabel('x'); ylabel('y'); axis equal; xlim([0 Lx]); ylim([0 Ly]);

    subplot(2,3,4);
    surf(X_mesh, Y_mesh, eta_mf12_3rd_term, 'EdgeColor', 'none');
    view(2); colormap jet; colorbar; shading interp;
    title('MF12 Original (3rd Order Term)');
    xlabel('x'); ylabel('y'); axis equal; xlim([0 Lx]); ylim([0 Ly]);
    c_lim = clim;

    subplot(2,3,5);
    surf(X_mesh, Y_mesh, eta_stream_3rd_term, 'EdgeColor', 'none');
    view(2); colormap jet; colorbar; shading interp;
    title('Flow+Spectral (3rd Order Term)');
    xlabel('x'); ylabel('y'); axis equal; xlim([0 Lx]); ylim([0 Ly]);
    clim(c_lim);

    subplot(2,3,6);
    surf(X_mesh, Y_mesh, abs(diff), 'EdgeColor', 'none');
    view(2); colormap jet; colorbar; shading interp;
    title(sprintf('Diff 3rd Term (Max: %.2e)', max_diff));
    xlabel('x'); ylabel('y'); axis equal; xlim([0 Lx]); ylim([0 Ly]);

    saveas(gcf, 'compare_3rd_order_crossing_full.png');
    fprintf('Saved figure compare_3rd_order_crossing_full.png\n');

    % --- New Plot: Centerline & Diagonal Comparison ---
    figure('Position',[100 100 1000 600]);
    
    % 1. Centerline (y = Ly/2)
    mid_y_idx = floor(Ny/2) + 1;
    y_val_center = y_vec(mid_y_idx);
    
    eta_mf12_center = eta_mf12_3rd_term(mid_y_idx, :);
    eta_stream_center = eta_stream_3rd_term(mid_y_idx, :);
    
    subplot(2,1,1);
    plot(x_vec/lambda, eta_mf12_center, 'b-', 'LineWidth', 1.5, 'DisplayName', 'MF12 Original');
    hold on;
    plot(x_vec/lambda, eta_stream_center, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Flow+Spectral');
    plot(x_vec/lambda, abs(eta_mf12_center - eta_stream_center), 'k:', 'LineWidth', 1.0, 'DisplayName', 'Abs Diff');
    
    grid on; legend('Location', 'best');
    title(sprintf('Centerline Profile (3rd Order Term) at y = %.1f m', y_val_center));
    xlabel('x / \lambda'); ylabel('\eta^{(33)} [m]');
    
    % 2. Diagonal (0,0) to (Lx, Ly)
    % Grid might be non-square in index count (Nx=128, Ny=32)
    num_diag = 200;
    x_diag = linspace(0, Lx, num_diag);
    y_diag = linspace(0, Ly, num_diag);
    diag_dist = sqrt(x_diag.^2 + y_diag.^2);
    
    eta_mf12_diag = interp2(X_mesh, Y_mesh, eta_mf12_3rd_term, x_diag, y_diag, 'linear');
    eta_stream_diag = interp2(X_mesh, Y_mesh, eta_stream_3rd_term, x_diag, y_diag, 'linear');
    
    subplot(2,1,2);
    plot(diag_dist/lambda, eta_mf12_diag, 'b-', 'LineWidth', 1.5, 'DisplayName', 'MF12 Original');
    hold on;
    plot(diag_dist/lambda, eta_stream_diag, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Flow+Spectral');
    plot(diag_dist/lambda, abs(eta_mf12_diag - eta_stream_diag), 'k:', 'LineWidth', 1.0, 'DisplayName', 'Abs Diff');
    
    grid on; legend('Location', 'best');
    title('Diagonal Profile (3rd Order Term) (0,0) to (Lx, Ly)');
    xlabel('Distance along diagonal / \lambda'); ylabel('\eta^{(33)} [m]');
    
    saveas(gcf, 'compare_3rd_order_crossing_lines.png');
    fprintf('Saved figure compare_3rd_order_crossing_lines.png\n');
    
    % --- New Plot: 2nd Order Comparison ---
    figure('Position',[150 150 1500 400], 'Visible', 'off');
    subplot(1,3,1);
    surf(X_mesh, Y_mesh, eta_stream_2nd_term, 'EdgeColor', 'none');
    view(2); colormap jet; colorbar; shading interp;
    title('Full MF12 (2nd Order Term)'); xlabel('x'); ylabel('y'); axis equal; xlim([0 Lx]); ylim([0 Ly]);
    c_lim_2 = clim;

    subplot(1,3,2);
    surf(X_mesh, Y_mesh, eta_sup_2nd_term, 'EdgeColor', 'none');
    view(2); colormap jet; colorbar; shading interp;
    title('Modified Super+Sub (2nd Order Term)'); xlabel('x'); ylabel('y'); axis equal; xlim([0 Lx]); ylim([0 Ly]);
    clim(c_lim_2);

    subplot(1,3,3);
    surf(X_mesh, Y_mesh, abs(diff_2nd_sup), 'EdgeColor', 'none');
    view(2); colormap jet; colorbar; shading interp;
    title(sprintf('Diff 2nd (Max: %.2e)', max_diff_2nd_sup));
    xlabel('x'); ylabel('y'); axis equal; xlim([0 Lx]); ylim([0 Ly]);
    
    saveas(gcf, 'compare_2nd_order_verification.png');
    fprintf('Saved figure compare_2nd_order_verification.png\n');
    
    % --- New Plot: Centerline & Diagonal Comparison (2nd Order) ---
    figure('Position',[150 150 1000 600], 'Visible', 'off');
    
    % 1. Centerline (y = Ly/2)
    % Re-use mid_y_idx, y_val_center from previous block
    
    eta_stream_2nd_center = eta_stream_2nd_term(mid_y_idx, :);
    eta_sup_2nd_center = eta_sup_2nd_term(mid_y_idx, :);
    
    subplot(2,1,1);
    plot(x_vec/lambda, eta_stream_2nd_center, 'b-', 'LineWidth', 1.5, 'DisplayName', 'MF12 Integrated (Ref)');
    hold on;
    plot(x_vec/lambda, eta_sup_2nd_center, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Modified Super+Sub');
    plot(x_vec/lambda, abs(eta_stream_2nd_center - eta_sup_2nd_center), 'k:', 'LineWidth', 1.0, 'DisplayName', 'Abs Diff');
    
    grid on; legend('Location', 'best');
    title(sprintf('Centerline Profile (2nd Order Term) at y = %.1f m', y_val_center));
    xlabel('x / \lambda'); ylabel('\eta^{(2)} [m]');
    
    % 2. Diagonal (0,0) to (Lx, Ly)
    % Re-use x_diag, y_diag, diag_dist
    
    eta_stream_2nd_diag = interp2(X_mesh, Y_mesh, eta_stream_2nd_term, x_diag, y_diag, 'linear');
    eta_sup_2nd_diag = interp2(X_mesh, Y_mesh, eta_sup_2nd_term, x_diag, y_diag, 'linear');
    
    subplot(2,1,2);
    plot(diag_dist/lambda, eta_stream_2nd_diag, 'b-', 'LineWidth', 1.5, 'DisplayName', 'MF12 Integrated (Ref)');
    hold on;
    plot(diag_dist/lambda, eta_sup_2nd_diag, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Modified Super+Sub');
    plot(diag_dist/lambda, abs(eta_stream_2nd_diag - eta_sup_2nd_diag), 'k:', 'LineWidth', 1.0, 'DisplayName', 'Abs Diff');
    
    grid on; legend('Location', 'best');
    title('Diagonal Profile (2nd Order Term) (0,0) to (Lx, Ly)');
    xlabel('Distance along diagonal / \lambda'); ylabel('\eta^{(2)} [m]');
    
    saveas(gcf, 'compare_2nd_order_crossing_lines.png');
    fprintf('Saved figure compare_2nd_order_crossing_lines.png\n');
    
    % --- New Plot: Total Nonlinear Surface Comparison (Full vs Modified) ---
    figure('Position',[150 150 1000 600], 'Visible', 'off');
    
    % 1. Centerline (y = Ly/2)
    eta_full_center = eta_stream_3(mid_y_idx, :);
    eta_mod_center = eta_sup_3(mid_y_idx, :);
    
    subplot(2,1,1);
    plot(x_vec/lambda, eta_full_center, 'b-', 'LineWidth', 1.5, 'DisplayName', 'MF12 Full (Ref)');
    hold on;
    plot(x_vec/lambda, eta_mod_center, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Modified (Partial 3rd)');
    plot(x_vec/lambda, abs(eta_full_center - eta_mod_center), 'k:', 'LineWidth', 1.0, 'DisplayName', 'Abs Diff');
    
    grid on; legend('Location', 'best');
    title(sprintf('Centerline Profile (Total Surface) at y = %.1f m', y_val_center));
    xlabel('x / \lambda'); ylabel('\eta_{total} [m]');
    
    % 2. Diagonal (0,0) to (Lx, Ly)
    eta_full_diag = interp2(X_mesh, Y_mesh, eta_stream_3, x_diag, y_diag, 'linear');
    eta_mod_diag = interp2(X_mesh, Y_mesh, eta_sup_3, x_diag, y_diag, 'linear');
    
    subplot(2,1,2);
    plot(diag_dist/lambda, eta_full_diag, 'b-', 'LineWidth', 1.5, 'DisplayName', 'MF12 Full (Ref)');
    hold on;
    plot(diag_dist/lambda, eta_mod_diag, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Modified (Partial 3rd)');
    plot(diag_dist/lambda, abs(eta_full_diag - eta_mod_diag), 'k:', 'LineWidth', 1.0, 'DisplayName', 'Abs Diff');
    
    grid on; legend('Location', 'best');
    title('Diagonal Profile (Total Surface) (0,0) to (Lx, Ly)');
    xlabel('Distance along diagonal / \lambda'); ylabel('\eta_{total} [m]');
    
    saveas(gcf, 'compare_total_surface_lines.png');
    fprintf('Saved figure compare_total_surface_lines.png\n');
    
    % --- New Plot: 1st Order (Linear) Comparison ---
    figure('Position',[150 150 1000 600], 'Visible', 'off');
    
    eta_lin_center_mf12 = eta_mf12_1(mid_y_idx, :);
    eta_lin_center_sup = eta_sup_1(mid_y_idx, :);
    
    subplot(2,1,1);
    plot(x_vec/lambda, eta_lin_center_mf12, 'b-', 'LineWidth', 1.5, 'DisplayName', 'MF12 Full (Reference)');
    hold on;
    plot(x_vec/lambda, eta_lin_center_sup, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Modified Super');
    plot(x_vec/lambda, abs(eta_lin_center_mf12 - eta_lin_center_sup), 'k:', 'LineWidth', 1.0, 'DisplayName', 'Abs Diff');
    
    grid on; legend('Location', 'best');
    title(sprintf('Centerline Profile (1st Order Linear) at y = %.1f m', y_val_center));
    xlabel('x / \lambda'); ylabel('\eta^{(1)} [m]');
    
    eta_lin_diag_mf12 = interp2(X_mesh, Y_mesh, eta_mf12_1, x_diag, y_diag, 'linear');
    eta_lin_diag_sup = interp2(X_mesh, Y_mesh, eta_sup_1, x_diag, y_diag, 'linear');
    
    subplot(2,1,2);
    plot(diag_dist/lambda, eta_lin_diag_mf12, 'b-', 'LineWidth', 1.5, 'DisplayName', 'MF12 Full (Reference)');
    hold on;
    plot(diag_dist/lambda, eta_lin_diag_sup, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Modified Super');
    plot(diag_dist/lambda, abs(eta_lin_diag_mf12 - eta_lin_diag_sup), 'k:', 'LineWidth', 1.0, 'DisplayName', 'Abs Diff');
    
    grid on; legend('Location', 'best');
    title('Diagonal Profile (1st Order Linear) (0,0) to (Lx, Ly)');
    xlabel('Distance along diagonal / \lambda'); ylabel('\eta^{(1)} [m]');
    
    saveas(gcf, 'compare_1st_order_lines.png');
    fprintf('Saved figure compare_1st_order_lines.png\n');
    
catch ME
    printf('Plotting error: %s\n', ME.message);
end


catch ME_outer; fid = fopen('compare_fatal_error.txt','w'); fprintf(fid, 'Fatal Error: %s\n', ME_outer.message); fclose(fid); exit(1); end
