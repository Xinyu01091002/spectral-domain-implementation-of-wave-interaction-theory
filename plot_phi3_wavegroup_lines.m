clc; clear; close all;
addpath(genpath(fullfile(pwd, 'irregularWavesMF12')));

% Wave-group style directional case
g = 9.81;
h = 100;
Ux = 0;
Uy = 0;
Lx = 3000;
Ly = 3000;
Nx = 32;
Ny = 32;
t = 0.0;

% Build a directional wave-group-like spectrum in k-space
rng(1234);
dkx = 2*pi/Lx;
dky = 2*pi/Ly;

kx_idx_all = (-floor(Nx/2)):(ceil(Nx/2)-1);
ky_idx_all = (-floor(Ny/2)):(ceil(Ny/2)-1);
[KXI, KYI] = meshgrid(kx_idx_all, ky_idx_all);
kx_all = KXI(:) * dkx;
ky_all = KYI(:) * dky;
kmag_all = hypot(kx_all, ky_all);
theta_all = atan2(ky_all, kx_all);

keep = (kx_all > 0) | (kx_all == 0 & ky_all > 0);
kx_all = kx_all(keep);
ky_all = ky_all(keep);
kmag_all = kmag_all(keep);
theta_all = theta_all(keep);

Tp = 10;
kp = (2*pi/Tp)^2 / g;
theta0 = pi/4;           % 45 deg mean direction
sig_k = 0.12*kp;         % narrow-band => wave-group behavior
sig_t = deg2rad(12);     % directional spread

Sk = exp(-0.5*((kmag_all - kp)/sig_k).^2);
St = exp(-0.5*(angle(exp(1i*(theta_all - theta0)))/sig_t).^2);
W = Sk .* St;

% Keep energetic modes
[~, idx_sort] = sort(W, 'descend');
N = 30;
idx = idx_sort(1:N);
kx = kx_all(idx);
ky = ky_all(idx);

% Amplitudes and focusing phase at domain center
A0 = 0.02;
amp = A0 * W(idx) / max(W(idx));
xf = Lx/2;
yf = Ly/2;
omega_lin = sqrt(g*hypot(kx,ky).*tanh(h*hypot(kx,ky)));
phase = -(kx*xf + ky*yf) + omega_lin*t;
a = amp .* cos(phase);
b = amp .* sin(phase);

% Grid
x = (0:Nx-1) * (Lx/Nx);
y = (0:Ny-1) * (Ly/Ny);
[X, Y] = meshgrid(x, y);

% Compute phi (order 2 and 3) for direct and spectral
c2 = coeffsMF12(2, g, h, a, b, kx, ky, Ux, Uy);
c3 = coeffsMF12(3, g, h, a, b, kx, ky, Ux, Uy);

[~, phi2_d] = surfaceMF12_new(2, c2, X, Y, t);
[~, phi3_d] = surfaceMF12_new(3, c3, X, Y, t);
phi3term_d = phi3_d - phi2_d;

[~, phi2_s] = surfaceMF12_spectral(c2, Lx, Ly, Nx, Ny, t);
[~, phi3_s] = surfaceMF12_spectral(c3, Lx, Ly, Nx, Ny, t);
phi3term_s = phi3_s - phi2_s;

% Centerline: y = Ly/2
[~, iyc] = min(abs(y - Ly/2));
phi_c_d = phi3term_d(iyc, :);
phi_c_s = phi3term_s(iyc, :);
phi_c_e = phi_c_d - phi_c_s;

% Diagonal: (0,0) -> (L,L)
Ldiag = min(Lx, Ly);
s = linspace(0, Ldiag, Nx);
phi_diag_d = interp2(X, Y, phi3term_d, s, s, 'linear');
phi_diag_s = interp2(X, Y, phi3term_s, s, s, 'linear');
phi_diag_e = phi_diag_d - phi_diag_s;

% Report
fprintf('Wave-group phi3-term compare:\n');
fprintf('  Centerline max |diff| = %.6e\n', max(abs(phi_c_e)));
fprintf('  Diagonal  max |diff| = %.6e\n', max(abs(phi_diag_e)));

% Plot lines only
fig = figure('Color','w','Position',[100 100 1300 780]);
tiledlayout(2,2);

nexttile;
plot(x, phi_c_d, 'k-', 'LineWidth', 1.2); hold on;
plot(x, phi_c_s, 'r--', 'LineWidth', 1.2);
grid on;
title('\phi_3-\phi_2 Centerline');
xlabel('x (m)'); ylabel('\phi (m^2/s)');
legend({'Direct','Spectral'}, 'Location','best');

nexttile;
plot(x, phi_c_e, 'b-', 'LineWidth', 1.2);
grid on;
title(sprintf('Centerline Diff (max %.2e)', max(abs(phi_c_e))));
xlabel('x (m)'); ylabel('\Delta\phi (m^2/s)');

nexttile;
plot(s, phi_diag_d, 'k-', 'LineWidth', 1.2); hold on;
plot(s, phi_diag_s, 'r--', 'LineWidth', 1.2);
grid on;
title('\phi_3-\phi_2 Diagonal');
xlabel('diag distance (m)'); ylabel('\phi (m^2/s)');
legend({'Direct','Spectral'}, 'Location','best');

nexttile;
plot(s, phi_diag_e, 'b-', 'LineWidth', 1.2);
grid on;
title(sprintf('Diagonal Diff (max %.2e)', max(abs(phi_diag_e))));
xlabel('diag distance (m)'); ylabel('\Delta\phi (m^2/s)');

out_png = ['phi3_wavegroup_lines_' datestr(now,'yyyymmdd_HHMMSS') '.png'];
saveas(fig, out_png);
fprintf('Saved: %s\n', out_png);

