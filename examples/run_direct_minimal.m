clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
run(fullfile(rootDir, 'setup_paths.m'));

g = 9.81;
h = 40;
Ux = 0;
Uy = 0;
order = 2;

Lx = 400;
Ly = 400;
Nx = 64;
Ny = 64;
t = 0;

kx = [0.080, 0.075, 0.070];
ky = [0.000, 0.012, -0.010];
a = [0.18, 0.10, 0.08];
b = [0.00, 0.04, -0.02];

x = (0:Nx-1) * (Lx/Nx);
y = (0:Ny-1) * (Ly/Ny);
[X, Y] = meshgrid(x, y);
coeffs = mf12_direct_coefficients(order, g, h, a, b, kx, ky, Ux, Uy);
[eta, phi] = mf12_direct_surface(order, coeffs, X, Y, t);

figure('Color', 'w');
subplot(1,2,1);
imagesc(x, y, eta);
axis image;
set(gca, 'YDir', 'normal');
title('Direct eta');
xlabel('x (m)');
ylabel('y (m)');
colorbar;

subplot(1,2,2);
imagesc(x, y, phi);
axis image;
set(gca, 'YDir', 'normal');
title('Direct phi');
xlabel('x (m)');
ylabel('y (m)');
colorbar;

fprintf('Direct example complete. max|eta| = %.3e, max|phi| = %.3e\n', max(abs(eta(:))), max(abs(phi(:))));
