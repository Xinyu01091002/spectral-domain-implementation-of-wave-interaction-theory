clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
matlabDir = fileparts(scriptDir);
run(fullfile(matlabDir, 'setup_paths.m'));

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

coeffs = mf12_spectral_coefficients(order, g, h, a, b, kx, ky, Ux, Uy, struct('enable_subharmonic', false));
[eta, phi, X, Y] = mf12_spectral_surface(coeffs, Lx, Ly, Nx, Ny, t);

figure('Color', 'w');
subplot(1,2,1);
imagesc(X(1,:), Y(:,1), eta);
axis image;
set(gca, 'YDir', 'normal');
title('Spectral eta');
xlabel('x (m)');
ylabel('y (m)');
colorbar;

subplot(1,2,2);
imagesc(X(1,:), Y(:,1), phi);
axis image;
set(gca, 'YDir', 'normal');
title('Spectral phi');
xlabel('x (m)');
ylabel('y (m)');
colorbar;

fprintf('Spectral example complete. max|eta| = %.3e, max|phi| = %.3e\n', max(abs(eta(:))), max(abs(phi(:))));
