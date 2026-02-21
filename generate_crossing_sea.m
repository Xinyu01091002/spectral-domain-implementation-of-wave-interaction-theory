clc;clear ;close all;
addpath(genpath('irregularWavesMF12'));

% generate_crossing_sea.m
% Generates a crossing sea state with two directional wave groups.
% This script creates two wave systems with different directions and superimposes them.

% 1. Basic Physics & Grid Setup
g = 9.81;
h = 150; 
Lx = 3000; 
Ly = 3000; 
Nx = 256;  
Ny = 64;   % Sufficient resolution for directional spreading

dx = Lx/Nx;
dy = Ly/Ny;
x_vec = (0:Nx-1)*dx;
y_vec = (0:Ny-1)*dy;
[X_mesh, Y_mesh] = meshgrid(x_vec, y_vec);

dkx = 2*pi/Lx;
dky = 2*pi/Ly;

% 2. Wave Parameters

% Frequency Spectrum Parameters (JONSWAP-like or Gaussian)
Tp = 12; 
kp = (2*pi/Tp)^2/g; 
Akp = 0.12;        % Overall Steepness parameter
alpha =1;         % Asymmetry for frequency spectrum
kw_left = 0.004606; 
kw_right = kw_left * alpha;

% Target Amplitude
A_focus_total = Akp / kp; 

% Crossing Sea Parameters
% System 1
theta_1_deg = 30;           % Angle for 1st system [deg]
spread_1_deg = 5;          % Spreading width [deg]
weight_1 = 0.5;             % Relative weight

% System 2
theta_2_deg = -30;          % Angle for 2nd system [deg]
spread_2_deg = 5;          % Spreading width [deg]
weight_2 = 0.5;             % Relative weight

% Focus Time & Location
tf = .0; 
xf = Lx/2; 
yf = Ly/2;

fprintf('Crossing Sea Generation:\n');
fprintf('  System 1: Angle = %.1f deg, Spread = %.1f deg, Weight = %.2f\n', theta_1_deg, spread_1_deg, weight_1);
fprintf('  System 2: Angle = %.1f deg, Spread = %.1f deg, Weight = %.2f\n', theta_2_deg, spread_2_deg, weight_2);
fprintf('  Total A_focus = %.3f m (Akp=%.2f)\n', A_focus_total, Akp);

% 3. Spectral Generation on Grid
[Idx_x, Idx_y] = meshgrid(-Nx/2 : Nx/2-1, -Ny/2 : Ny/2-1);
KX_grid = Idx_x * dkx;
KY_grid = Idx_y * dky;

% Filter: Exclude only DC component. 
% We allow negative KX to support waves from all directions.
K_mag = sqrt(KX_grid.^2 + KY_grid.^2);
% valid_mask = K_mag > 0 & K_mag < (Nx/2)*dkx; % Remove Nyquist or high freq if needed, mainly remove 0
valid_mask = K_mag > 0;
KX = KX_grid(valid_mask);
KY = KY_grid(valid_mask);
K_mag = K_mag(valid_mask);
Theta = atan2(KY, KX); % Result in (-pi, pi]

% Frequency Spectrum S(k)
kw_vec = ones(size(K_mag)) * kw_right;
kw_vec(K_mag <= kp) = kw_left;
S_vec = exp(-(K_mag - kp).^2 ./ (2*kw_vec.^2));

% Directional Spreading D(theta) - Bimodal
theta_1 = deg2rad(theta_1_deg);
theta_2 = deg2rad(theta_2_deg);
w_1 = deg2rad(spread_1_deg);
w_2 = deg2rad(spread_2_deg);

% Angular difference function for correct wrapping [-pi, pi]
ang_diff = @(t, t_ref) angle(exp(1i*(t - t_ref)));

D_vec_1 = exp(-(ang_diff(Theta, theta_1)).^2 / (2*w_1^2));
D_vec_2 = exp(-(ang_diff(Theta, theta_2)).^2 / (2*w_2^2));

% Combined Raw Amplitude Weights
% We combine them linearly.
Amp_raw = S_vec .* (weight_1 * D_vec_1 + weight_2 * D_vec_2);

% 4. Energy Based Cutoff
Energy_raw = Amp_raw.^2;
[sorted_E, sort_idx] = sort(Energy_raw, 'descend');
cum_E = cumsum(sorted_E);
total_E = sum(Energy_raw);
cutoff_energy_ratio = 0.99; 

last_idx = find(cum_E >= cutoff_energy_ratio * total_E, 1, 'first');
if isempty(last_idx), last_idx = length(cum_E); end

keep_indices = sort_idx(1:last_idx);
temp_kx = KX(keep_indices);
temp_ky = KY(keep_indices);
temp_amp = Amp_raw(keep_indices);

fprintf('  -> Selected %d wave components (%.1f%% Energy)\n', length(temp_amp), cutoff_energy_ratio*100);

% 5. Scaling to Target Amplitude (Linear Focusing)
if isempty(temp_amp)
    warning('No components generated.');
else
    scale_factor = A_focus_total / sum(temp_amp);
    temp_amp = temp_amp * scale_factor;
end

% 6. Phase Generation for Focusing
kx_list = [];
ky_list = [];
a_list = [];
b_list = [];

for i = 1:length(temp_kx)
    kk_x = temp_kx(i);
    kk_y = temp_ky(i);
    kk = sqrt(kk_x^2 + kk_y^2);
    
    % Linear dispersion relation
    ww = sqrt(g * kk * tanh(kk * h)); 
    
    % Focusing Phase: phi = -kx*xf - ky*yf + w*tf
    phase = -kk_x * xf - kk_y * yf + ww * tf; 
    
    kx_list(end+1,1) = kk_x;
    ky_list(end+1,1) = kk_y;
    a_list(end+1,1) = temp_amp(i) * cos(phase);
    b_list(end+1,1) = temp_amp(i) * sin(phase);
end

Ux = 0; Uy = 0;

% 7. Compute Surface (Linear + 2nd Order) using Spectral Method
fprintf('Computing 2nd Order Coefficients...\n');
coeffs2 = coeffsMF12(2, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy);

fprintf('Computing Surface Elevation (Linear + 2nd Order)...\n');
t_eval = tf; 
[eta_spec_2, ~] = surfaceMF12_spectral(coeffs2, Lx, Ly, Nx, Ny, t_eval);

% 8. Plotting
fig = figure('Position',[100 100 1200 500]);

subplot(1,2,1);
scatter(kx_list, ky_list, 20, temp_amp, 'filled');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');
title('Wave Number Spectrum (Amplitudes)');
axis equal; grid on;
xlim([-0.1 0.1]); ylim([-0.1 0.1]); % Zoom in near origin

subplot(1,2,2);
surf(X_mesh, Y_mesh, eta_spec_2, 'EdgeColor', 'none');
view(2); colormap jet; colorbar;
shading interp;
clim([-A_focus_total, A_focus_total]*1);
title(sprintf('Surface Elevation at t=%.1fs (Crossing Sea)', t_eval));
xlabel('x [m]'); ylabel('y [m]');
axis equal; axis([0 Lx 0 Ly]);

saveas(fig, 'check_crossing_sea_result.png');
fprintf('Done. Result check_crossing_sea_result.png saved.\n');
