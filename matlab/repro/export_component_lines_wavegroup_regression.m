clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
matlabDir = fileparts(scriptDir);
repoDir = fileparts(matlabDir);
run(fullfile(matlabDir, 'setup_paths.m'));

caseDir = fullfile(repoDir, 'cross_language_comparison', 'cases', 'wavegroup_regression');
caseData = jsondecode(fileread(fullfile(caseDir, 'case.json')));

a = readmatrix(fullfile(caseDir, caseData.arrays.a));
b = readmatrix(fullfile(caseDir, caseData.arrays.b));
kx = readmatrix(fullfile(caseDir, caseData.arrays.kx));
ky = readmatrix(fullfile(caseDir, caseData.arrays.ky));
inp = caseData.inputs;

g = inp.g;
h = inp.h;
Ux = inp.Ux;
Uy = inp.Uy;
Lx = inp.Lx;
Ly = inp.Ly;
Nx = inp.Nx;
Ny = inp.Ny;
t_eval = inp.t;

x = readmatrix(fullfile(caseDir, caseData.reference.arrays.x)).';
y = readmatrix(fullfile(caseDir, caseData.reference.arrays.y));
[X, Y] = meshgrid(x, y);

c1 = mf12_spectral_coefficients(1, g, h, a, b, kx, ky, Ux, Uy);
c2 = mf12_spectral_coefficients(2, g, h, a, b, kx, ky, Ux, Uy);
c3 = mf12_spectral_coefficients(3, g, h, a, b, kx, ky, Ux, Uy);
c3.third_order_subharmonic_mode = 'skip';

c2_d = mf12_direct_coefficients(2, g, h, a, b, kx, ky, Ux, Uy);
k2_super = [2*c2_d.kappa(:); c2_d.kappa_npm(1:2:end).'];
k2_sub = c2_d.kappa_npm(2:2:end).';
k2_super = k2_super(isfinite(k2_super) & k2_super > 0);
k2_sub = k2_sub(isfinite(k2_sub) & k2_sub >= 0);
k2_cut = 0.5 * (max(k2_sub) + min(k2_super));

phi1 = phi_surface_spectral(c1, Lx, Ly, Nx, Ny, t_eval);
phi2_total = phi_surface_spectral(c2, Lx, Ly, Nx, Ny, t_eval);
phi2_inc = phi2_total - phi1;
phi2sup = split_by_wavenumber(phi2_inc, Lx, Ly, k2_cut, 'high');
phi2sub = split_by_wavenumber(phi2_inc, Lx, Ly, k2_cut, 'low');
phi3_total = phi_surface_spectral(c3, Lx, Ly, Nx, Ny, t_eval);
phi3 = phi3_total - phi2_total;

[~, iyc] = min(abs(y - Ly/2));
Ldiag = min(Lx, Ly);
s = linspace(0, Ldiag, Nx);

center = [phi1(iyc,:); phi2sup(iyc,:); phi2sub(iyc,:); phi3(iyc,:)];
diagv = [ ...
    interp2(X, Y, phi1, s, s, 'linear'); ...
    interp2(X, Y, phi2sup, s, s, 'linear'); ...
    interp2(X, Y, phi2sub, s, s, 'linear'); ...
    interp2(X, Y, phi3, s, s, 'linear')];

outDir = fullfile(caseDir, 'reference', 'matlab_components');
if ~exist(outDir, 'dir'), mkdir(outDir); end
writematrix(center, fullfile(outDir, 'centerline_phi_components.csv'));
writematrix(diagv, fullfile(outDir, 'diagonal_phi_components.csv'));
writematrix(x(:), fullfile(outDir, 'x.csv'));
writematrix(s(:), fullfile(outDir, 's.csv'));

meta = struct();
meta.labels = {'first harmonic', 'second superharmonic', 'second subharmonic', 'third superharmonic'};
meta.reference = 'MATLAB spectral';
meta.k2_cut = k2_cut;
fid = fopen(fullfile(outDir, 'components.json'), 'w');
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', jsonencode(meta, 'PrettyPrint', true));

fprintf('Saved MATLAB component lines to: %s\n', outDir);

function phi = phi_surface_spectral(coeffs, Lx, Ly, Nx, Ny, t)
    [~, phi] = mf12_spectral_surface(coeffs, Lx, Ly, Nx, Ny, t);
end

function phi_part = split_by_wavenumber(phi, Lx, Ly, k_cut, mode)
    [Ny, Nx] = size(phi);
    dkx = 2*pi / Lx;
    dky = 2*pi / Ly;
    kx_idx = [0:(Nx/2-1), -Nx/2:-1];
    ky_idx = [0:(Ny/2-1), -Ny/2:-1];
    [KX, KY] = meshgrid(kx_idx * dkx, ky_idx * dky);
    K = hypot(KX, KY);
    spec = fft2(phi);
    switch mode
        case 'high'
            mask = K >= k_cut;
        case 'low'
            mask = K < k_cut;
        otherwise
            error('Unknown filter mode: %s', mode);
    end
    phi_part = real(ifft2(spec .* mask));
end
