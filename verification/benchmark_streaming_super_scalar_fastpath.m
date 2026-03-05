% BENCHMARK_STREAMING_SUPER_SCALAR_FASTPATH
% Compare surfaceMF12_streaming_super speed:
%   A) scalar fast path ON  (new)
%   B) scalar fast path OFF (old accumarray-every-call behavior)

clc;
clear;
close all;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(genpath(fullfile(rootDir, 'irregularWavesMF12')));
setenv('MF12_PROGRESS', '0');

% -------------------- Case setup --------------------
g = 9.81;
h = 150;
Ux = 0;
Uy = 0;
order = 3;

Lx = 1500;
Ly = 1500;
Nx = 256;
Ny = 256;
t = 0.0;

N_components = 300;  % adjust (e.g. 60/90/120)
Tp = 12;
kp = (2*pi/Tp)^2 / g;
theta1 = deg2rad(25);
theta2 = deg2rad(-35);
sig_k = 0.12 * kp;
sig_t1 = deg2rad(10);
sig_t2 = deg2rad(12);
w1 = 0.55;
w2 = 0.45;
Hs_target = 3.5;

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

Sk = exp(-0.5*((kmag_all - kp)/sig_k).^2);
St1 = exp(-0.5*(angle(exp(1i*(theta_all - theta1)))/sig_t1).^2);
St2 = exp(-0.5*(angle(exp(1i*(theta_all - theta2)))/sig_t2).^2);
W = Sk .* (w1*St1 + w2*St2);

valid = W > 1e-8 * max(W);
kx_all = kx_all(valid);
ky_all = ky_all(valid);
W = W(valid);
[~, idx_sort] = sort(W, 'descend');
N_components = min(N_components, numel(idx_sort));
idx = idx_sort(1:N_components);
kx = kx_all(idx); kx = kx(:).';
ky = ky_all(idx); ky = ky(:).';
Wsel = W(idx); Wsel = Wsel(:).';

amp_raw = Wsel / max(Wsel);
m0_raw = 0.5 * sum(amp_raw.^2);
Hs_raw = 4 * sqrt(max(m0_raw, eps));
amp = amp_raw * (Hs_target / Hs_raw);

xf = 0.5 * Lx;
yf = 0.5 * Ly;
omega_lin = sqrt(g*hypot(kx,ky).*tanh(h*hypot(kx,ky)));
phi_focus = -(kx*xf + ky*yf) + omega_lin*t;
a = amp .* cos(phi_focus); a = a(:).';
b = amp .* sin(phi_focus); b = b(:).';

fprintf('=== Benchmark: streaming_super scalar fast path ===\n');
fprintf('N=%d, grid=%dx%d, order=%d\n', N_components, Nx, Ny, order);

nWarm = 1;
nRun = 3;

% -------------------- OLD behavior (disable fast path) --------------------
setenv('MF12_DISABLE_SCALAR_FAST', '1');
for i = 1:nWarm
    surfaceMF12_streaming_super(order, g, h, a, b, kx, ky, Ux, Uy, Lx, Ly, Nx, Ny, t);
end
tic;
for i = 1:nRun
    [eta_old, phi_old] = surfaceMF12_streaming_super(order, g, h, a, b, kx, ky, Ux, Uy, Lx, Ly, Nx, Ny, t);
end
t_old = toc / nRun;

% -------------------- NEW behavior (fast path) --------------------
setenv('MF12_DISABLE_SCALAR_FAST', '0');
for i = 1:nWarm
    surfaceMF12_streaming_super(order, g, h, a, b, kx, ky, Ux, Uy, Lx, Ly, Nx, Ny, t);
end
tic;
for i = 1:nRun
    [eta_new, phi_new] = surfaceMF12_streaming_super(order, g, h, a, b, kx, ky, Ux, Uy, Lx, Ly, Nx, Ny, t);
end
t_new = toc / nRun;

% -------------------- Accuracy check --------------------
max_eta_diff = max(abs(eta_old(:) - eta_new(:)));
max_phi_diff = max(abs(phi_old(:) - phi_new(:)));

fprintf('OLD avg time: %.3f s\n', t_old);
fprintf('NEW avg time: %.3f s\n', t_new);
fprintf('Speedup (OLD/NEW): %.2f x\n', t_old / max(t_new, eps));
fprintf('Max |eta_old-eta_new|: %.3e\n', max_eta_diff);
fprintf('Max |phi_old-phi_new|: %.3e\n', max_phi_diff);

% reset environment flag
setenv('MF12_DISABLE_SCALAR_FAST', '');
