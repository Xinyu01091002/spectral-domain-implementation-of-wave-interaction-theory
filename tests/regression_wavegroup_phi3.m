clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
run(fullfile(rootDir, 'setup_paths.m'));

fprintf('=== Regression Test: Directional Wave-Group phi3 ===\n');

g = 9.81;
Tp = 12;
t = 5 * Tp;
kp = (2*pi/Tp)^2 / g;
h = 1 / kp;
Ux = 0;
Uy = 0;
Lx = 2000;
Ly = 2000;
Nx = 64;
Ny = 64;

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
sig_t1 = deg2rad(20);
sig_t2 = deg2rad(24);
w1 = 0.55;
w2 = 0.45;

Sk = exp(-0.5*((kmag_all - kp)/sig_k).^2);
St1 = exp(-0.5*(angle(exp(1i*(theta_all - theta1)))/sig_t1).^2);
St2 = exp(-0.5*(angle(exp(1i*(theta_all - theta2)))/sig_t2).^2);
W = Sk .* (w1*St1 + w2*St2);

[~, idx_sort] = sort(W, 'descend');
N = 40;
idx = idx_sort(1:N);
kx = kx_all(idx);
ky = ky_all(idx);

A0 = 0.25;
amp = A0 * W(idx) / max(W(idx));
xf = Lx/2;
yf = Ly/2;
omega_lin = sqrt(g*hypot(kx,ky).*tanh(h*hypot(kx,ky)));
phase = -(kx*xf + ky*yf) + omega_lin*t;
a = amp .* cos(phase);
b = amp .* sin(phase);

x = (0:Nx-1) * (Lx/Nx);
y = (0:Ny-1) * (Ly/Ny);
[~, iyc] = min(abs(y - Ly/2));

c2_d = mf12_direct_coefficients(2, g, h, a, b, kx, ky, Ux, Uy);
c3_d = mf12_direct_coefficients(3, g, h, a, b, kx, ky, Ux, Uy);
c2_s = mf12_spectral_coefficients(2, g, h, a, b, kx, ky, Ux, Uy);
c3_s = mf12_spectral_coefficients(3, g, h, a, b, kx, ky, Ux, Uy);
c3_s.third_order_subharmonic_mode = 'skip';
c3_d_sup = select_third_order_superharmonic(c3_d);

[X, Y] = meshgrid(x, y);
[~, phi2_d] = mf12_direct_surface(2, c2_d, X, Y, t);
[~, phi3_d] = mf12_direct_surface(3, c3_d_sup, X, Y, t);
[~, phi2_s] = mf12_spectral_surface(c2_s, Lx, Ly, Nx, Ny, t);
[~, phi3_s] = mf12_spectral_surface(c3_s, Lx, Ly, Nx, Ny, t);

phi3sup_d = phi3_d - phi2_d;
phi3sup_s = phi3_s - phi2_s;
center_err = max(abs(phi3sup_d(iyc,:) - phi3sup_s(iyc,:)));
field_err = max(abs(phi3sup_d(:) - phi3sup_s(:)));

fprintf('centerline max error = %.3e\n', center_err);
fprintf('field max error      = %.3e\n', field_err);

tol_center = 1e-10;
tol_field = 1e-9;
if center_err > tol_center || field_err > tol_field
    error('Regression test failed: center_err=%.3e, field_err=%.3e', center_err, field_err);
end

fprintf('Regression test passed.\n');

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
