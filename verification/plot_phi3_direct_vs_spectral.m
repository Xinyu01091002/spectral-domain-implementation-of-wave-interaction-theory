clc; clear; close all;
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
run(fullfile(rootDir, 'setup_paths.m'));
outDir = fullfile(rootDir, 'outputs');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Config
g = 9.81;
h = 100;
Ux = 0;
Uy = 0;
Lx = 3000;
Ly = 3000;
Nx = 64;
Ny = 64;
t = 0.0;
N = 80;
rng(77);

% Build a stable directional test case from unique half-plane modes.
dkx = 2*pi/Lx;
dky = 2*pi/Ly;
kx_idx_all = (-floor(Nx/2)):(ceil(Nx/2)-1);
ky_idx_all = (-floor(Ny/2)):(ceil(Ny/2)-1);
[KXI, KYI] = meshgrid(kx_idx_all, ky_idx_all);
KXI = KXI(:);
KYI = KYI(:);
keep = (KXI > 0) | (KXI == 0 & KYI > 0);
KXI = KXI(keep);
KYI = KYI(keep);
pick = randperm(numel(KXI), N).';
kx = KXI(pick) * dkx;
ky = KYI(pick) * dky;
kmag = hypot(kx, ky);
amp = 0.03 * exp(-0.2 * (kmag / max(kmag)));
ph = 2*pi * rand(N, 1);
a = amp .* cos(ph);
b = amp .* sin(ph);

x = (0:Nx-1) * (Lx/Nx);
y = (0:Ny-1) * (Ly/Ny);

% Compute coefficients using the two public workflows.
c2_d = mf12_direct_coefficients(2, g, h, a, b, kx, ky, Ux, Uy);
c3_d = mf12_direct_coefficients(3, g, h, a, b, kx, ky, Ux, Uy);
c2_s = mf12_spectral_coefficients(2, g, h, a, b, kx, ky, Ux, Uy);
c3_s = mf12_spectral_coefficients(3, g, h, a, b, kx, ky, Ux, Uy);

% Direct method
[X, Y] = meshgrid(x, y);
[~, phi2_d] = mf12_direct_surface(2, c2_d, X, Y, t);
[~, phi3_d] = mf12_direct_surface(3, c3_d, X, Y, t);
phi3term_d = phi3_d - phi2_d;

% Spectral method
[~, phi2_s] = mf12_spectral_surface(c2_s, Lx, Ly, Nx, Ny, t);
[~, phi3_s] = mf12_spectral_surface(c3_s, Lx, Ly, Nx, Ny, t);
phi3term_s = phi3_s - phi2_s;

% Differences
d_total = phi3_d - phi3_s;
d_term = phi3term_d - phi3term_s;

max_total = max(abs(d_total(:)));
max_term = max(abs(d_term(:)));
rel_total = max_total / max(max(abs(phi3_d(:))), eps);
rel_term = max_term / max(max(abs(phi3term_d(:))), eps);

fprintf('phi order3 total: max diff = %.6e, rel = %.6e\n', max_total, rel_total);
fprintf('phi pure 3rd term: max diff = %.6e, rel = %.6e\n', max_term, rel_term);

% Plot
fig = figure('Color', 'w', 'Position', [80 80 1500 900]);
tiledlayout(2, 3);

nexttile;
imagesc(x, y, phi3_d); axis image; set(gca, 'YDir', 'normal');
title('Direct \phi_{order=3}');
xlabel('x'); ylabel('y'); colorbar;

nexttile;
imagesc(x, y, phi3_s); axis image; set(gca, 'YDir', 'normal');
title('Spectral \phi_{order=3}');
xlabel('x'); ylabel('y'); colorbar;

nexttile;
imagesc(x, y, d_total); axis image; set(gca, 'YDir', 'normal');
title(sprintf('Diff total (max %.2e)', max_total));
xlabel('x'); ylabel('y'); colorbar;

nexttile;
imagesc(x, y, phi3term_d); axis image; set(gca, 'YDir', 'normal');
title('Direct (\phi_3 - \phi_2)');
xlabel('x'); ylabel('y'); colorbar;

nexttile;
imagesc(x, y, phi3term_s); axis image; set(gca, 'YDir', 'normal');
title('Spectral (\phi_3 - \phi_2)');
xlabel('x'); ylabel('y'); colorbar;

nexttile;
imagesc(x, y, d_term); axis image; set(gca, 'YDir', 'normal');
title(sprintf('Diff pure-3rd (max %.2e)', max_term));
xlabel('x'); ylabel('y'); colorbar;

out_png = fullfile(outDir, ['phi3_direct_vs_spectral_' datestr(now, 'yyyymmdd_HHMMSS') '.png']);
saveas(fig, out_png);
fprintf('Saved: %s\n', out_png);
