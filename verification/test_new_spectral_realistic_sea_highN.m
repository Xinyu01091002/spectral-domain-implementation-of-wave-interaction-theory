% TEST_NEW_SPECTRAL_REALISTIC_SEA_HIGHN
% New spectral-only test using many wave-number components
% (single in-memory coefficient precompute + spectral reconstruction).

clc;
clear;
close all;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(genpath(fullfile(rootDir, 'irregularWavesMF12')));

outDir = fullfile(rootDir, 'outputs');
if ~exist(outDir, 'dir'), mkdir(outDir); end
setenv('MF12_PROGRESS', '1');

% -------------------- Sea state / domain (based on plot_phi3_wavegroup_lines) --------------------
g = 9.81;
h = 150;
Ux = 0;
Uy = 0;

Lx = 3000;
Ly = 3000;
Nx = 64;
Ny = 64;

% Many components (raise to 240/300 if your machine allows longer runtimes).
N_components = 300;
order = 3;  % MF12 third-order superharmonic pipeline

% Peak period and directional crossing setup
Tp = 12;
kp = (2*pi/Tp)^2 / g;
lambda_p = 2*pi / kp;
theta1 = deg2rad(25);
theta2 = deg2rad(-35);
sig_k = 0.12 * kp;
sig_t1 = deg2rad(10);
sig_t2 = deg2rad(12);
w1 = 0.55;
w2 = 0.45;

% Realistic sea severity target (linear level): Hs ~ 3.5 m
Hs_target = 3.5;

% -------------------- Build candidate k-space and select energetic modes --------------------
dkx = 2*pi/Lx;
dky = 2*pi/Ly;
kx_idx_all = (-floor(Nx/2)):(ceil(Nx/2)-1);
ky_idx_all = (-floor(Ny/2)):(ceil(Ny/2)-1);
[KXI, KYI] = meshgrid(kx_idx_all, ky_idx_all);

kx_all = KXI(:) * dkx;
ky_all = KYI(:) * dky;
kmag_all = hypot(kx_all, ky_all);
theta_all = atan2(ky_all, kx_all);

% Keep half-plane to avoid duplicate conjugate modes
keep = (kx_all > 0) | (kx_all == 0 & ky_all > 0);
kx_all = kx_all(keep);
ky_all = ky_all(keep);
kmag_all = kmag_all(keep);
theta_all = theta_all(keep);

% Radial + directional weights (wave-group-like crossing sea)
Sk = exp(-0.5*((kmag_all - kp)/sig_k).^2);
St1 = exp(-0.5*(angle(exp(1i*(theta_all - theta1)))/sig_t1).^2);
St2 = exp(-0.5*(angle(exp(1i*(theta_all - theta2)))/sig_t2).^2);
W = Sk .* (w1*St1 + w2*St2);

% Remove very low-energy bins and select strongest modes
valid = W > 1e-8 * max(W);
kx_all = kx_all(valid);
ky_all = ky_all(valid);
W = W(valid);

[~, idx_sort] = sort(W, 'descend');
N_components = min(N_components, numel(idx_sort));
idx = idx_sort(1:N_components);
kx = kx_all(idx);
ky = ky_all(idx);
kx = kx(:).';
ky = ky(:).';
Wsel = W(idx);
Wsel = Wsel(:).';

% -------------------- Amplitudes and directional wave-group focusing phases --------------------
amp_raw = Wsel / max(Wsel);

% Scale amplitudes to match target linear Hs (m0 ~= sum(A_n^2)/2)
m0_raw = 0.5 * sum(amp_raw.^2);
Hs_raw = 4 * sqrt(max(m0_raw, eps));
amp = amp_raw * (Hs_target / Hs_raw);

% Deterministic wave-group phase (no random phase):
% focus at (xf, yf) and time t_focus so that components align there.
xf = 0.5 * Lx;
yf = 0.5 * Ly;
t_focus = 0.0;
omega_lin = sqrt(g*hypot(kx,ky).*tanh(h*hypot(kx,ky)));
phi_focus = -(kx*xf + ky*yf) + omega_lin*t_focus;
a = amp .* cos(phi_focus);
b = amp .* sin(phi_focus);
a = a(:).';
b = b(:).';

assert(isrow(a) && isrow(b) && isrow(kx) && isrow(ky), 'a,b,kx,ky must be row vectors');
assert(numel(a)==numel(kx) && numel(b)==numel(kx), 'a,b,kx,ky length mismatch');

% -------------------- New spectral method only (RAM workflow) --------------------
opts = struct();
opts.enable_subharmonic = false;  % explicit for clarity

fprintf('\n=== New Spectral-Only High-N Test ===\n');
fprintf('N components: %d\n', N_components);
fprintf('Target Hs (linear): %.2f m\n', Hs_target);
fprintf('Peak wavelength: %.1f m\n', lambda_p);
fprintf('Directional wave-group focus: (x,y,t) = (%.1f, %.1f, %.1f)\n', xf, yf, t_focus);

t_coeff = tic;
coeffs = coeffsMF12_superharmonic(order, g, h, a, b, kx, ky, Ux, Uy, opts);
fprintf('Coefficient stage done in %.2f s\n', toc(t_coeff));
fprintf('Storage mode: in-RAM\n');

% Reconstruct at multiple times
t_list = [0];
eta_all = cell(size(t_list));
phi_all = cell(size(t_list));
for it = 1:numel(t_list)
    tt = t_list(it);
    fprintf('\nReconstruct at t = %.2f s ...\n', tt);
    t_rec = tic;
    [eta, phi, X, Y] = surfaceMF12_spectral(coeffs, Lx, Ly, Nx, Ny, tt);
    fprintf('  done in %.2f s\n', toc(t_rec));
    eta_all{it} = eta;
    phi_all{it} = phi;
    fprintf('  eta range: [%.3f, %.3f] m\n', min(eta(:)), max(eta(:)));
    fprintf('  phi range: [%.3f, %.3f] m^2/s\n', min(phi(:)), max(phi(:)));
end

% -------------------- Quick visualization --------------------
x = X(1,:);
y = Y(:,1);
[~, iyc] = min(abs(y - Ly/2));

fig = figure('Color', 'w', 'Position', [120 120 1200 520]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
for it = 1:numel(t_list)
    nexttile;
    imagesc(x, y, eta_all{it}); axis image; set(gca, 'YDir', 'normal');
    title(sprintf('\\eta, t=%.2fs', t_list(it)));
    xlabel('x (m)'); ylabel('y (m)'); colorbar;

    nexttile;
    plot(x, eta_all{it}(iyc,:), 'b-', 'LineWidth', 1.2); grid on;
    title(sprintf('centerline \\eta, t=%.2fs', t_list(it)));
    xlabel('x (m)'); ylabel('\eta (m)');
end

out_png = fullfile(outDir, 'mf12_new_spectral_realistic_highN.png');
exportgraphics(fig, out_png, 'Resolution', 150);
fprintf('\nSaved figure: %s\n', out_png);
fprintf('Test finished.\n');
