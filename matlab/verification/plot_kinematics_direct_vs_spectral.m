clc; clear; close all;
scriptDir = fileparts(mfilename('fullpath'));
matlabDir = fileparts(scriptDir);
repoDir = fileparts(matlabDir);
run(fullfile(matlabDir, 'setup_paths.m'));
outDir = fullfile(repoDir, 'outputs');
if ~exist(outDir, 'dir'), mkdir(outDir); end

set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.0);
set(groot, 'defaultTextFontSize', 12);

% Config
g = 9.81;
Ux = 0;
Uy = 0;
Lx = 3000;
Ly = 3000;
Nx = 64;
Ny = 64;
t = 0.0;
z_eval = -20.0; % Constant-z horizontal plane for the spectral reconstruction.

% Wave-group-style directional case, aligned with plot_phi3_wavegroup_lines.m
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

kp = 0.0279;
h = 1 / kp;
Akp = 0.35;
kw_left = 0.004606;
kw_right = kw_left;
A_focus_total = Akp / kp;
t_focus = 0.0;

kw_vec = kw_right * ones(size(kmag_all));
kw_vec(kmag_all <= kp) = kw_left;
Sk = exp(-((kmag_all - kp).^2) ./ (2 * kw_vec.^2));

group1_heading_deg = 45;
group1_spread_deg = 5;
group1_weight = 0.5;

group2_heading_deg = -45;
group2_spread_deg = 5;
group2_weight = 0.5;

D1 = gaussian_spreading(theta_all - deg2rad(group1_heading_deg), group1_spread_deg);
D2 = gaussian_spreading(theta_all - deg2rad(group2_heading_deg), group2_spread_deg);
W = Sk .* (group1_weight * D1 + group2_weight * D2);

energy_keep_frac = 0.99;
max_components = 80; % Keep the direct third-order verification practical.
[W_sorted, idx_sort] = sort(W, 'descend');
cum_energy = cumsum(W_sorted);
N_energy = find(cum_energy >= energy_keep_frac * cum_energy(end), 1, 'first');
N = min(N_energy, max_components);
idx = idx_sort(1:N);
kx = kx_all(idx);
ky = ky_all(idx);

Wsel = W(idx);
amp = Wsel;
amp = amp * (A_focus_total / max(sum(amp), eps));
xf = Lx/2;
yf = Ly/2;
omega_lin = sqrt(g*hypot(kx,ky).*tanh(h*hypot(kx,ky)));
phase = -(kx*xf + ky*yf) + omega_lin*t_focus;
a = amp .* cos(phase);
b = amp .* sin(phase);

x = (0:Nx-1) * (Lx/Nx);
y = (0:Ny-1) * (Ly/Ny);
[X, Y] = meshgrid(x, y);

% Compute coefficients using the public workflows.
c3_d = mf12_direct_coefficients(3, g, h, a, b, kx, ky, Ux, Uy);
c3_s = mf12_spectral_coefficients(3, g, h, a, b, kx, ky, Ux, Uy);
c3_s.third_order_subharmonic_mode = 'skip';
c3_d_sup = select_third_order_superharmonic(c3_d);

% Direct and spectral kinematics on the same constant-z plane.
[u_d, v_d, w_d, p_d, phi_d, uV_d, vV_d, a_x_d, a_y_d] = kinematicsMF12(3, c3_d_sup, X, Y, z_eval, t);
[u_s, v_s, w_s, p_s, phi_s, uV_s, vV_s, a_x_s, a_y_s] = mf12_spectral_kinematics(c3_s, Lx, Ly, Nx, Ny, z_eval, t);

fields = {
    'phi', phi_d, phi_s;
    'u',   u_d,   u_s;
    'v',   v_d,   v_s;
    'w',   w_d,   w_s;
    'p',   p_d,   p_s;
    'uV',  uV_d,  uV_s;
    'vV',  vV_d,  vV_s;
    'a_x', a_x_d, a_x_s;
    'a_y', a_y_d, a_y_s;
};

fprintf('Constant-z plane comparison at z = %.3f m\n', z_eval);
fprintf('Wave-group case: retained %d components (requested %.2f%% cumulative energy, cap %d)\n', ...
    N, 100 * energy_keep_frac, max_components);
for i = 1:size(fields, 1)
    name = fields{i, 1};
    directField = fields{i, 2};
    spectralField = fields{i, 3};
    diffField = directField - spectralField;
    maxAbs = max(abs(diffField(:)));
    relAbs = maxAbs / max(max(abs(directField(:))), eps);
    fprintf('%s: max diff = %.6e, rel = %.6e\n', name, maxAbs, relAbs);
end

field_order = {'phi','u','v','w','p','uV','vV','a_x','a_y'};

% Figure 1: directional wave-group fields only
fig1 = figure('Color', 'w', 'Position', [40 40 1700 1350]);
tl1 = tiledlayout(3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numel(field_order)
    name = field_order{i};
    idxf = find(strcmp(fields(:,1), name), 1, 'first');
    spectralField = fields{idxf, 3};

    ax = nexttile;
    imagesc(x, y, spectralField);
    axis image;
    set(ax, 'YDir', 'normal');
    xlabel('x (m)');
    ylabel('y (m)');
    title(strrep(name, '_', '\_'));
    cb = colorbar;
    cb.Label.String = name;
end
title(tl1, sprintf('Directional Wave-Group Kinematics on z = %.3f m', z_eval), ...
    'FontSize', 18, 'FontWeight', 'bold');
out_fields = fullfile(outDir, 'kinematics_directional_fields.png');
exportgraphics(fig1, out_fields, 'Resolution', 220);
fprintf('Saved: %s\n', out_fields);

% Figure 2: error fields, highlighting machine-precision agreement
fig2 = figure('Color', 'w', 'Position', [40 40 1700 1350]);
tl2 = tiledlayout(3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numel(field_order)
    name = field_order{i};
    idxf = find(strcmp(fields(:,1), name), 1, 'first');
    directField = fields{idxf, 2};
    spectralField = fields{idxf, 3};
    diffField = directField - spectralField;
    maxAbs = max(abs(diffField(:)));
    ref = max(max(abs(directField(:))), eps);
    errNorm = diffField / ref;
    clim = max(abs(errNorm(:)));
    if clim <= 0
        clim = 1e-16;
    end

    ax = nexttile;
    imagesc(x, y, errNorm);
    axis image;
    set(ax, 'YDir', 'normal');
    xlabel('x (m)');
    ylabel('y (m)');
    title(sprintf('%s error', strrep(name, '_', '\_')));
    colormap(ax, parula);
    caxis(ax, [-clim, clim]);
    cb = colorbar;
    ylabel(cb, sprintf('(%s_d - %s_s) / max|%s_d|', name, name, name), 'Interpreter', 'none');
    text(ax, 0.98, 0.96, sprintf('max|err| = %.2e', maxAbs), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 2, 'FontSize', 11);
end
title(tl2, sprintf('Directional Wave-Group Error Maps on z = %.3f m', z_eval), ...
    'FontSize', 18, 'FontWeight', 'bold');
out_err = fullfile(outDir, 'kinematics_directional_errors.png');
exportgraphics(fig2, out_err, 'Resolution', 220);
fprintf('Saved: %s\n', out_err);

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

function D = gaussian_spreading(theta, spread_angle_deg)
    theta_wrapped = angle(exp(1i * theta));
    sigma = deg2rad(spread_angle_deg);
    D = exp(-0.5 * (theta_wrapped / max(sigma, eps)).^2);
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
