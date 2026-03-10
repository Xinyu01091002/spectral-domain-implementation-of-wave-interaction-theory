clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(genpath(fullfile(rootDir, 'irregularWavesMF12')));

fprintf('=== Smoke Test: Directional Wave-Group MF12 ===\n');

g = 9.81;
Tp = 12;
t = 2 * Tp;
kp = (2*pi/Tp)^2 / g;
h = 1 / kp;
Ux = 0;
Uy = 0;
Lx = 2000;
Ly = 2000;
Nx = 32;
Ny = 32;

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

theta1 = deg2rad(25);
theta2 = deg2rad(-35);
sig_k = 0.10 * kp;
sig_t1 = deg2rad(16);
sig_t2 = deg2rad(18);
w1 = 0.55;
w2 = 0.45;

Sk = exp(-0.5*((kmag_all - kp)/sig_k).^2);
St1 = exp(-0.5*(angle(exp(1i*(theta_all - theta1)))/sig_t1).^2);
St2 = exp(-0.5*(angle(exp(1i*(theta_all - theta2)))/sig_t2).^2);
W = Sk .* (w1*St1 + w2*St2);

[~, idx_sort] = sort(W, 'descend');
N = 20;
idx = idx_sort(1:N);
kx = kx_all(idx);
ky = ky_all(idx);

A0 = 0.18;
amp = A0 * W(idx) / max(W(idx));
xf = Lx/2;
yf = Ly/2;
omega_lin = sqrt(g*hypot(kx,ky).*tanh(h*hypot(kx,ky)));
phase = -(kx*xf + ky*yf) + omega_lin*t;
a = amp .* cos(phase);
b = amp .* sin(phase);

x = (0:Nx-1) * (Lx/Nx);
y = (0:Ny-1) * (Ly/Ny);
[X, Y] = meshgrid(x, y);

c1_d = coeffsMF12(1, g, h, a, b, kx, ky, Ux, Uy);
c2_d = coeffsMF12(2, g, h, a, b, kx, ky, Ux, Uy);
c3_d = coeffsMF12(3, g, h, a, b, kx, ky, Ux, Uy);
c1_s = coeffsMF12_superharmonic(1, g, h, a, b, kx, ky, Ux, Uy);
c2_s = coeffsMF12_superharmonic(2, g, h, a, b, kx, ky, Ux, Uy);
c3_s = coeffsMF12_superharmonic(3, g, h, a, b, kx, ky, Ux, Uy);
c3_s.third_order_subharmonic_mode = 'skip';
c3_d_sup = select_third_order_superharmonic(c3_d);

[~, phi1_d] = surfaceMF12_new(1, c1_d, X, Y, t);
[~, phi2_d] = surfaceMF12_new(2, c2_d, X, Y, t);
[~, phi3_d] = surfaceMF12_new(3, c3_d_sup, X, Y, t);
[~, phi1_s] = surfaceMF12_spectral(c1_s, Lx, Ly, Nx, Ny, t);
[~, phi2_s] = surfaceMF12_spectral(c2_s, Lx, Ly, Nx, Ny, t);
[~, phi3_s] = surfaceMF12_spectral(c3_s, Lx, Ly, Nx, Ny, t);

phi1_err = max(abs(phi1_d(:) - phi1_s(:)));
phi3_err = max(abs((phi3_d(:) - phi2_d(:)) - (phi3_s(:) - phi2_s(:))));

fprintf('max|phi1_d - phi1_s| = %.3e\n', phi1_err);
fprintf('max|phi3sup_d - phi3sup_s| = %.3e\n', phi3_err);

tol = 1e-9;
if phi1_err > tol || phi3_err > tol
    error('Smoke test failed: phi1_err=%.3e, phi3_err=%.3e, tol=%.3e', phi1_err, phi3_err, tol);
end

fprintf('Smoke test passed.\n');

function coeffs_out = select_third_order_superharmonic(coeffs_in)
    coeffs_out = coeffs_in;

    if isfield(coeffs_out, 'G_np2m')
        mask_minus = false(size(coeffs_out.G_np2m));
        mask_minus(2:2:end) = true;
        coeffs_out = zero_fields(coeffs_out, mask_minus, ...
            {'A_np2m','B_np2m','F_np2m','G_np2m','mu_np2m','kappa_np2m','omega_np2m','kx_np2m','ky_np2m'});
    end

    if isfield(coeffs_out, 'G_2npm')
        mask_minus = false(size(coeffs_out.G_2npm));
        mask_minus(2:2:end) = true;
        coeffs_out = zero_fields(coeffs_out, mask_minus, ...
            {'A_2npm','B_2npm','F_2npm','G_2npm','mu_2npm','kappa_2npm','omega_2npm','kx_2npm','ky_2npm'});
    end

    if isfield(coeffs_out, 'G_npmpp') && isfield(coeffs_out, 'N')
        mask_sub = true(size(coeffs_out.G_npmpp));
        idx = 0;
        for n = 1:coeffs_out.N
            for m = n+1:coeffs_out.N
                for pmm = [1 -1]
                    for p = m+1:coeffs_out.N
                        for pmp = [1 -1]
                            idx = idx + 1;
                            if pmm == 1 && pmp == 1
                                mask_sub(idx) = false;
                            end
                        end
                    end
                end
            end
        end
        coeffs_out = zero_fields(coeffs_out, mask_sub, ...
            {'A_npmpp','B_npmpp','F_npmpp','G_npmpp','mu_npmpp','kappa_npmpp','omega_npmpp','kx_npmpp','ky_npmpp'});
    end
end

function coeffs_out = zero_fields(coeffs_in, mask, field_names)
    coeffs_out = coeffs_in;
    for k = 1:numel(field_names)
        name = field_names{k};
        if isfield(coeffs_out, name)
            vals = coeffs_out.(name);
            vals(mask) = 0;
            coeffs_out.(name) = vals;
        end
    end
end
