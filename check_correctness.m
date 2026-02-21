clc;clear ;close all;
addpath(genpath('irregularWavesMF12'));

% check_correctness.m
% 验证 surfaceMF12_spectral 与 surfaceMF12_new 结果的一致性
% 并测试表面生成的正确性

% 1. 基本参数设置
g = 9.81;
h = 150; % 水深
Tp = 12; % 峰值周期
kp = (2*pi/Tp)^2/g; % 深水近似，或者用迭代码求
% 使用 Linear_focus_envelope 中的类似逻辑生成 Gaussian 谱
% 为了方便对比 FFT，我们需要确保物理域和波数域是离散匹配的

Lx = 3000; % 物理域长度 [m]
Ly = 3000; % 物理域宽度 [m] (假设正方形)
Nx = 256;  % Reduced for faster validation loop (Direct method is slow)
Ny = 16;

dx = Lx/Nx;
dy = Ly/Ny;
x_vec = (0:Nx-1)*dx;
y_vec = (0:Ny-1)*dy;
[X_mesh, Y_mesh] = meshgrid(x_vec, y_vec);

dkx = 2*pi/Lx;
dky = 2*pi/Ly;

% 2. 生成高斯波谱成分 (只在一维方向上分布，方便对比，或者简单的二维分布)
% 这里的目的是生成 dummy coeffs 数据
% 为了严格的 FFT 对比，波数必须正好落在 dkx 的整数倍上

k_modes_x = -Nx/2 : Nx/2-1; % 模式数
k_modes_y = -Ny/2 : Ny/2-1;

% 选取一部分非零能量的模态
% 这里的逻辑：只选取正方向的传播
N_wave = 0;
kx_list = [];
ky_list = [];
a_list = [];
b_list = [];

% 高斯谱参数 (Updated to match Custom Linear_focus_envelope logic)
kp = 0.0279;       % 峰值波数 [rad/m]
Akp = 0.12;        % 波陡 (Steepness parameter k*A)
alpha = 8;         % 谱宽不对称系数 (控制右侧谱宽)
tf = 60.0;         % 聚焦时间 [s]

% 定义非对称谱宽
kw_left = 0.004606; 
kw_right = kw_left * alpha; % 右侧谱宽随 alpha 变化

% 根据波陡计算目标聚焦波幅 (Linear Amplitude Sum)
% A_focus * kp = Akp  =>  A_focus = Akp / kp
A_focus = Akp / kp; 

fprintf('Setup Parameters: kp=%.4f rad/m, Akp=%.2f, alpha=%.1f\n', kp, Akp, alpha);
fprintf('  -> Grid Size: Nx=%d, Ny=%d\n', Nx, Ny);
fprintf('  -> Left Bandwidth kw_L=%.6f, Right Bandwidth kw_R=%.6f\n', kw_left, kw_right);
fprintf('  -> Target Focus Amplitude A_focus=%.3f m\n', A_focus);

% Directional Spreading Parameters (Roughly implemented as requested)
spreading_deg = 2;            % Spreading width [deg]
theta_w = deg2rad(spreading_deg);
theta_p = 0;                   % Mean direction 0 (x-axis)
fprintf('  -> Directional Spreading: %.1f deg\n', spreading_deg);

% 聚焦位置 (Center of the 2D domain)
xf = Lx/2; 
yf = Ly/2;

% Iterate over 2D wavenumber grid using Vectorization for sorting
[Idx_x, Idx_y] = meshgrid(-Nx/2 : Nx/2-1, -Ny/2 : Ny/2-1);
KX_grid = Idx_x * dkx;
KY_grid = Idx_y * dky;

% Filter: only forward (+x) waves
valid_mask = KX_grid > 0;
KX = KX_grid(valid_mask);
KY = KY_grid(valid_mask);

% Calculate Spectrum Weights
K_mag = sqrt(KX.^2 + KY.^2);
Theta = atan2(KY, KX);

% 1. Frequency Spectrum S(k) (Asymmetric Gaussian)
% Vectorized S calc
kw_vec = ones(size(K_mag)) * kw_right;
kw_vec(K_mag <= kp) = kw_left;
S_vec = exp(-(K_mag - kp).^2 ./ (2*kw_vec.^2));

% 2. Directional Spreading D(theta)
D_vec = exp(-(Theta - theta_p).^2 / (2*theta_w^2));

% Raw Amplitude Weights
Amp_raw = S_vec .* D_vec;

% 3. Energy Based Cutoff (Retain 99% Energy)
% Energy is proportional to amplitude squared
Energy_raw = Amp_raw.^2;
[sorted_E, sort_idx] = sort(Energy_raw, 'descend');
cum_E = cumsum(sorted_E);
total_E = sum(Energy_raw);
cutoff_energy_ratio = 0.99;

% Find index where cumulative energy passes threshold
last_idx = find(cum_E >= cutoff_energy_ratio * total_E, 1, 'first');
if isempty(last_idx), last_idx = length(cum_E); end

% Select components
keep_indices = sort_idx(1:last_idx);

temp_kx = KX(keep_indices);
temp_ky = KY(keep_indices);
temp_amp = Amp_raw(keep_indices);

fprintf('  -> Energy Cutoff: Retaining %.1f%% energy.\n', cutoff_energy_ratio*100);
fprintf('  -> Components reduced from %d to %d\n', length(Amp_raw), length(temp_amp));
N_wave = length(temp_amp);

% 归一化振幅以匹配目标聚焦波高
if isempty(temp_amp)
    warning('没有生成任何波浪成分，请检查 kp 和 dkx。');
else
    % Note: A_focus is the target linear amplitude sum
    scale_factor = A_focus / sum(temp_amp);
    temp_amp = temp_amp * scale_factor;
end

% 设置相位以在此处聚焦: eta ~ cos(kx*x + ky*y - w*t + phase)
% 聚焦条件: (kx*xf + ky*yf - w*tf + phase) = 0
for i = 1:length(temp_kx)
    kk_x = temp_kx(i);
    kk_y = temp_ky(i);
    kk = sqrt(kk_x^2 + kk_y^2);
    ww = sqrt(g * kk * tanh(kk * h)); % Linear dispersion
    
    % Fixed phase for focusing at (xf, yf, tf)
    phase = -kk_x * xf - kk_y * yf + ww * tf; 
    
    kx_list(end+1,1) = kk_x;
    ky_list(end+1,1) = kk_y;
    
    a_list(end+1,1) = temp_amp(i) * cos(phase);
    b_list(end+1,1) = temp_amp(i) * sin(phase);
end

Ux = 0; Uy = 0;

fprintf('Generated Wave Components N = %d\n', N_wave);

% =========================================================================
% 3. Generate Coefficients (2nd and 3rd order for separation)
% =========================================================================
fprintf('Computing 2nd Order Coefficients (coeffs2)...\n');
coeffs2 = coeffsMF12(2, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy);

fprintf('Computing 3rd Order Coefficients (coeffs3)...\n');
coeffs3 = coeffsMF12(3, g, h, a_list, b_list, kx_list, ky_list, Ux, Uy);

% =========================================================================
% 4. Separation of 3rd Order Contribution (Total 3rd - Total 2nd)
% =========================================================================
t = 0; % Focus Time (relative to simulation start t=0, wait, t param is simulation time)
       % Ideally if we want to see focus at t=0 of simulation, we set tf properly during phase gen.
       % Currently we are just viewing t=0 snapshot.

% --- A. Spectral Method (FFT) ---
fprintf('\n--- Spectral Method (FFT) ---\n');
tic;
fprintf('  Computing 2nd Order Result...\n');
[eta_spec_2, ~] = surfaceMF12_spectral(coeffs2, Lx, Ly, Nx, Ny, t);

fprintf('  Computing 3rd Order Result...\n');
[eta_spec_3, ~] = surfaceMF12_spectral(coeffs3, Lx, Ly, Nx, Ny, t);

% Extract Pure 3rd Order Contribution
eta_spec_3rd_term = eta_spec_3 - eta_spec_2;
time_spec = toc;
fprintf('  >> Time (FFT): %.4f s\n', time_spec);


% --- B. Direct Method (Old) ---
fprintf('\n--- Direct Method (Direct Summation) ---\n');
if exist('surfaceMF12_new', 'file')==2
    tic;
    fprintf('  Computing 2nd Order Result (Iterating summation)...\n');
    [eta_direct_2, ~] = surfaceMF12_new(2, coeffs2, X_mesh, Y_mesh, t);
    
    fprintf('  Computing 3rd Order Result (This may take a while)...\n');
    [eta_direct_3, ~] = surfaceMF12_new(3, coeffs3, X_mesh, Y_mesh, t);
    
    % Extract Pure 3rd Order Contribution
    eta_direct_3rd_term = eta_direct_3 - eta_direct_2;
    time_direct = toc;
    fprintf('  >> Time (Direct): %.4f s\n', time_direct);
else
    error('Cannot find surfaceMF12_new.m');
end

% =========================================================================
% 5. Comparison and Plotting (Focus on 3rd Order)
% =========================================================================

% Calculate Error
diff_term = abs(eta_spec_3rd_term - eta_direct_3rd_term);
max_diff = max(diff_term(:));
rel_err = max_diff / (max(abs(eta_direct_3rd_term(:))) + 1e-15); % Avoid division by zero

fprintf('\n--------------------------------------\n');
fprintf('Validation of 3rd Order Contribution:\n');
fprintf('  Max Absolute Error: %.4e (m)\n', max_diff);
fprintf('  Relative Error:     %.4e\n', rel_err);
fprintf('  Time Direct:        %.4f s\n', time_direct);
fprintf('  Time FFT:           %.4f s\n', time_spec);
fprintf('  Speedup Factor:     %.2f x\n', time_direct / time_spec);
fprintf('--------------------------------------\n');

if rel_err < 1e-8
    fprintf('>>> SUCCESS: Results match perfectly (within numerical precision).\n');
else
    fprintf('>>> WARNING: Large discrepancy found. Please check implementation.\n');
end

try
    figure('Color','w', 'Position', [100, 100, 1400, 500]);
    
    % Normalization Parameters
    lambda_p = 2*pi/kp;
    x_norm = x_vec / lambda_p;
    
    % Plot Caption Info
    kph = kp * h;
    caption_str = sprintf('Parameters: A k_p=%.2f, k_p h=%.2f, \\alpha=%.1f, t_{focus}=%.1f s, Grid: %dx%d', Akp, kph, alpha, tf, Nx, Ny);

    % Find Centerline Index
    mid_idx = floor(Ny/2) + 1;
    y_mid = y_vec(mid_idx);

    % Find Peak Location (Center of Wave Envelope)
    % Use the 3rd order signal magnitude to identify the packet center
    [~, peak_idx] = max(abs(eta_direct_3rd_term(mid_idx, :)));
    x_peak_norm = x_norm(peak_idx);
    
    % Define view range: Peak +/- 6 lambda (User requested) or 3? 
    % Let's use 6 as requested in prompt "envelope maximum around 6 lambda area" (meaning +/- 3? or +/- 6?) 
    % "envelope maximum周围的6 lambda区域" -> Likely means total width 6 lambda, i.e., +/- 3 lambda.
    xlim_range = [x_peak_norm - 3, x_peak_norm + 3];

    % 1. 2D Map (Direct)
    subplot(1,3,1); 
    imagesc(x_norm, y_vec, eta_direct_3rd_term); axis xy; 
    title({'Direct Method', '3rd Order Contribution Only'}); 
    colorbar; xlabel('x / \lambda_p'); ylabel('y (m)');
    xlim(xlim_range);
    
    % 2. 2D Map (Spectral)
    subplot(1,3,2); 
    imagesc(x_norm, y_vec, eta_spec_3rd_term); axis xy; 
    title({'Spectral FFT', '3rd Order Contribution Only'}); 
    colorbar; xlabel('x / \lambda_p'); ylabel('y (m)');
    xlim(xlim_range);

    % 3. Centerline Comparison
    subplot(1,3,3); 
    yyaxis left; % Explicitly set left axis
    hold on;
    % Plot Direct
    h1 = plot(x_norm, eta_direct_3rd_term(mid_idx, :), 'k-', 'LineWidth', 2, 'DisplayName', 'Direct (3rd Only)');
    % Plot Spectral
    h2 = plot(x_norm, eta_spec_3rd_term(mid_idx, :), 'r--', 'LineWidth', 2, 'DisplayName', 'FFT (3rd Only)');
    ylabel('Elevation \eta^{(3)} (m)');
    
    % Plot Difference on Right Axis
    diff_line = eta_spec_3rd_term(mid_idx, :) - eta_direct_3rd_term(mid_idx, :);
    
    yyaxis right; % Switch to right axis
    h3 = plot(x_norm, diff_line, 'b:', 'LineWidth', 1, 'DisplayName', 'Difference');
    ylabel('Difference (m)', 'Color', 'b');
    set(gca, 'ycolor', 'b');
    
    title(['Centerline Comparison (y = ' num2str(y_mid) ' m)']);
    xlabel('x / \lambda_p'); 
    xlim(xlim_range);
    
    % Combine legends if possible, or just let Matlab handle it (it usually does fine with yyaxis)
    legend([h1, h2, h3], 'Location', 'best');
    grid on; box on;
    
    sgtitle({'Validation of 3rd Order Wave Contribution', caption_str}, 'FontSize', 12);
    
catch ME
    warning('Plotting failed: %s', ME.message);
    fprintf('Error Stack:\n');
    disp(ME.stack(1));
end
