clc; clear; close all;
scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
addpath(genpath(fullfile(rootDir, 'irregularWavesMF12')));
outDir = fullfile(rootDir, 'outputs');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Wave-group style directional case
g = 9.81;
h = 150;
Ux = 0;
Uy = 0;
Lx = 5000;
Ly = 5000;
Nx = 128;
Ny = 128;
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

Tp = 12;
t = 5 * Tp;
kp = (2*pi/Tp)^2 / g;
h = 1 / kp;              % set k_p d = 1
lambda_p = 2*pi / kp;
% Crossing-sea wave-group setup: two directional groups with different headings.
theta1 = deg2rad(25);    % group-1 mean direction
theta2 = deg2rad(-35);   % group-2 mean direction (crossing)
sig_k = 0.10*kp;         % narrow-band => wave-group behavior
sig_t1 = deg2rad(20);
sig_t2 = deg2rad(24);
w1 = 0.55;               % energy weights
w2 = 0.45;

Sk = exp(-0.5*((kmag_all - kp)/sig_k).^2);
St1 = exp(-0.5*(angle(exp(1i*(theta_all - theta1)))/sig_t1).^2);
St2 = exp(-0.5*(angle(exp(1i*(theta_all - theta2)))/sig_t2).^2);
W = Sk .* (w1*St1 + w2*St2);

% Keep energetic modes
[~, idx_sort] = sort(W, 'descend');
N = 60;
idx = idx_sort(1:N);
kx = kx_all(idx);
ky = ky_all(idx);

% Amplitudes and focusing phase at domain center
A0 = 0.35;
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

% Keep the original benchmark pairing:
%   direct   = coeffsMF12 + surfaceMF12_new
%   spectral = coeffsMF12_superharmonic + surfaceMF12_spectral
c1_d = coeffsMF12(1, g, h, a, b, kx, ky, Ux, Uy);
c2_d = coeffsMF12(2, g, h, a, b, kx, ky, Ux, Uy);
c3_d = coeffsMF12(3, g, h, a, b, kx, ky, Ux, Uy);
c1_s = coeffsMF12_superharmonic(1, g, h, a, b, kx, ky, Ux, Uy);
c2_s = coeffsMF12_superharmonic(2, g, h, a, b, kx, ky, Ux, Uy);
c3_s = coeffsMF12_superharmonic(3, g, h, a, b, kx, ky, Ux, Uy);
c3_s.third_order_subharmonic_mode = 'skip';
c3_d_sup = select_third_order_superharmonic(c3_d);

% Second-order super/sub separation by wavenumber filtering
k2_super = [2*c2_d.kappa(:); c2_d.kappa_npm(1:2:end).'];
k2_sub = c2_d.kappa_npm(2:2:end).';
k2_super = k2_super(isfinite(k2_super) & k2_super > 0);
k2_sub = k2_sub(isfinite(k2_sub) & k2_sub >= 0);
if isempty(k2_super) || isempty(k2_sub)
    error('Unable to build second-order spectral cutoff from MF12 coefficients.');
end
k2_cut = 0.5 * (max(k2_sub) + min(k2_super));

% Decomposed components (Direct)
phi1_d = phi_surface_direct(1, c1_d, X, Y, t);
phi2_d_total = phi_surface_direct(2, c2_d, X, Y, t);
phi2_d_inc = phi2_d_total - phi1_d;
phi2sup_d = split_by_wavenumber(phi2_d_inc, Lx, Ly, k2_cut, 'high');
phi2sub_d = split_by_wavenumber(phi2_d_inc, Lx, Ly, k2_cut, 'low');
phi3_d_total = phi_surface_direct(3, c3_d_sup, X, Y, t);
phi3_d = phi3_d_total - phi2_d_total;

% Decomposed components (Spectral)
phi1_s = phi_surface_spectral(c1_s, Lx, Ly, Nx, Ny, t);
phi2_s_total = phi_surface_spectral(c2_s, Lx, Ly, Nx, Ny, t);
phi2_s_inc = phi2_s_total - phi1_s;
phi2sup_s = split_by_wavenumber(phi2_s_inc, Lx, Ly, k2_cut, 'high');
phi2sub_s = split_by_wavenumber(phi2_s_inc, Lx, Ly, k2_cut, 'low');
phi3_s_total = phi_surface_spectral(c3_s, Lx, Ly, Nx, Ny, t);
phi3_s = phi3_s_total - phi2_s_total;

% Centerline/diagonal extraction
[~, iyc] = min(abs(y - Ly/2));
center_d = {
    phi1_d(iyc,:), phi2sup_d(iyc,:), phi2sub_d(iyc,:), phi3_d(iyc,:)
};
center_s = {
    phi1_s(iyc,:), phi2sup_s(iyc,:), phi2sub_s(iyc,:), phi3_s(iyc,:)
};

Ldiag = min(Lx, Ly);
s = linspace(0, Ldiag, Nx);
diag_d = {
    interp2(X, Y, phi1_d, s, s, 'linear'), ...
    interp2(X, Y, phi2sup_d, s, s, 'linear'), ...
    interp2(X, Y, phi2sub_d, s, s, 'linear'), ...
    interp2(X, Y, phi3_d, s, s, 'linear')
};
diag_s = {
    interp2(X, Y, phi1_s, s, s, 'linear'), ...
    interp2(X, Y, phi2sup_s, s, s, 'linear'), ...
    interp2(X, Y, phi2sub_s, s, s, 'linear'), ...
    interp2(X, Y, phi3_s, s, s, 'linear')
};
labels = {'first harmonic', 'second superharmonic', 'second subharmonic', 'third superharmonic'};

% Plot window: around envelope focus, in units of a few dominant wavelengths.
half_window_lambda = 4.0;  % show +/- 4 lambda_p around focus (broader context)
xwin = [max(0, xf - half_window_lambda*lambda_p), min(Lx, xf + half_window_lambda*lambda_p)];
s_focus = min(xf, yf);     % focus projected on diagonal distance
swin = [max(0, s_focus - half_window_lambda*lambda_p), min(Ldiag, s_focus + half_window_lambda*lambda_p)];
xn = x / lambda_p;
sn = s / lambda_p;
xwin_n = xwin / lambda_p;
swin_n = swin / lambda_p;

% Report
fprintf('Wave-group decomposed compare (Direct vs Spectral):\n');
for j = 1:4
    e_c = center_d{j} - center_s{j};
    e_d = diag_d{j} - diag_s{j};
    fprintf('  %-18s: center max|diff|=%.3e, diag max|diff|=%.3e\n', labels{j}, max(abs(e_c)), max(abs(e_d)));
end

% Publication-style line plots
c_direct = [0.06, 0.33, 0.62];
c_spec = [0.82, 0.30, 0.11];
c_err = [0.20, 0.20, 0.20];
lw_main = 2.1;
lw_err = 1.6;
fs_title = 12;
fs_axis = 11;
fs_tick = 10;
fs_legend = 10;
fs_note = 9;

% Normalize error by global max |phi_s| over all plotted traces
phi_ref_max = 0;
for j = 1:4
    phi_ref_max = max(phi_ref_max, max(abs(center_d{j})));
    phi_ref_max = max(phi_ref_max, max(abs(center_s{j})));
    phi_ref_max = max(phi_ref_max, max(abs(diag_d{j})));
    phi_ref_max = max(phi_ref_max, max(abs(diag_s{j})));
end
phi_ref_max = max(phi_ref_max, eps);

% Global limits for normalized error panels
errn_max = 0;
for j = 1:4
    errn_max = max(errn_max, max(abs((center_d{j} - center_s{j}) / phi_ref_max)));
    errn_max = max(errn_max, max(abs((diag_d{j} - diag_s{j}) / phi_ref_max)));
end
if errn_max <= 0
    err_ylim = 1e-14;
else
    err_ylim = 10^(ceil(log10(errn_max)));
    err_ylim = max(err_ylim, 1e-14);
end
eps_ref = eps(1.0);

% -------- Figure 1: Direct vs Spectral lines --------
fig1 = figure('Color','w','Units','centimeters','Position',[2 2 26 11.2]);
tl1 = tiledlayout(2,4, 'Padding', 'compact', 'TileSpacing', 'compact');

for j = 1:4
    ax = nexttile;
    h1 = plot(ax, xn, center_d{j}, '-', 'Color', c_direct, 'LineWidth', lw_main); hold(ax, 'on');
    h2 = plot(ax, xn, center_s{j}, '--', 'Color', c_spec, 'LineWidth', lw_main);
    grid(ax, 'on'); grid(ax, 'minor');
    xlim(ax, xwin_n);
    title(ax, labels{j}, 'FontName', 'Times New Roman', 'FontSize', fs_title, 'FontWeight', 'normal');
    if j == 1
        ylabel(ax, 'Centerline \phi_s (m^2/s)', 'FontName', 'Times New Roman', 'FontSize', fs_axis);
        legend(ax, [h1, h2], {'Direct', 'Spectral'}, 'Location', 'northwest', ...
            'FontName', 'Times New Roman', 'FontSize', fs_legend, 'Box', 'off');
        yl = ylim(ax); d = max(eps, yl(2)-yl(1)); ylim(ax, [yl(1)-0.10*d, yl(2)+0.22*d]);
    else
        yl = ylim(ax); d = max(eps, yl(2)-yl(1)); ylim(ax, [yl(1)-0.08*d, yl(2)+0.12*d]);
    end
    set(ax, 'FontName', 'Times New Roman', 'FontSize', fs_tick, 'LineWidth', 0.8, 'Box', 'on');
    ax.GridAlpha = 0.20; ax.MinorGridAlpha = 0.10;
    ax.XTickLabel = [];
end

for j = 1:4
    ax = nexttile;
    h1 = plot(ax, sn, diag_d{j}, '-', 'Color', c_direct, 'LineWidth', lw_main); hold(ax, 'on');
    h2 = plot(ax, sn, diag_s{j}, '--', 'Color', c_spec, 'LineWidth', lw_main);
    grid(ax, 'on'); grid(ax, 'minor');
    xlim(ax, swin_n);
    if j == 1
        ylabel(ax, 'Diagonal \phi_s (m^2/s)', 'FontName', 'Times New Roman', 'FontSize', fs_axis);
        legend(ax, [h1, h2], {'Direct', 'Spectral'}, 'Location', 'northwest', ...
            'FontName', 'Times New Roman', 'FontSize', fs_legend, 'Box', 'off');
        yl = ylim(ax); d = max(eps, yl(2)-yl(1)); ylim(ax, [yl(1)-0.10*d, yl(2)+0.22*d]);
    else
        yl = ylim(ax); d = max(eps, yl(2)-yl(1)); ylim(ax, [yl(1)-0.08*d, yl(2)+0.12*d]);
    end
    set(ax, 'FontName', 'Times New Roman', 'FontSize', fs_tick, 'LineWidth', 0.8, 'Box', 'on');
    ax.GridAlpha = 0.20; ax.MinorGridAlpha = 0.10;
end
xlabel(tl1, 'x/\lambda_p', 'FontName', 'Times New Roman', 'FontSize', fs_axis);

out_cmp = fullfile(outDir, 'mf12_phi3_wavegroup_lines_comparison_pub.png');
exportgraphics(fig1, out_cmp, 'Resolution', 600);
fprintf('Saved: %s\n', out_cmp);

% -------- Figure 2: Error-only panels --------
fig2 = figure('Color','w','Units','centimeters','Position',[2 2 26 11.2]);
tl2 = tiledlayout(2,4, 'Padding', 'compact', 'TileSpacing', 'compact');

for j = 1:4
    ax = nexttile;
    e = (center_d{j} - center_s{j}) / phi_ref_max;
    plot(ax, xn, e, '-', 'Color', c_err, 'LineWidth', lw_err); hold(ax, 'on');
    yline(ax, 0, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0);
    yline(ax, eps_ref, '--', 'Color', [0.25 0.55 0.25], 'LineWidth', 0.9);
    yline(ax, -eps_ref, '--', 'Color', [0.25 0.55 0.25], 'LineWidth', 0.9);
    grid(ax, 'on'); grid(ax, 'minor');
    xlim(ax, xwin_n); ylim(ax, [-err_ylim, err_ylim]);
    title(ax, labels{j}, 'FontName', 'Times New Roman', 'FontSize', fs_title, 'FontWeight', 'normal');
    if j == 1
        ylabel(ax, 'Centerline \Delta\phi_s / max|\phi_s|', 'FontName', 'Times New Roman', 'FontSize', fs_axis);
    end
    set(ax, 'FontName', 'Times New Roman', 'FontSize', fs_tick, 'LineWidth', 0.8, 'Box', 'on');
    ax.GridAlpha = 0.20; ax.MinorGridAlpha = 0.10;
    ax.XTickLabel = [];
    text(ax, 0.97, 0.90, sprintf('max=%.1e', max(abs(e))), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'FontName', 'Times New Roman', 'FontSize', fs_note, 'BackgroundColor', 'w');
end

for j = 1:4
    ax = nexttile;
    e = (diag_d{j} - diag_s{j}) / phi_ref_max;
    plot(ax, sn, e, '-', 'Color', c_err, 'LineWidth', lw_err); hold(ax, 'on');
    yline(ax, 0, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0);
    yline(ax, eps_ref, '--', 'Color', [0.25 0.55 0.25], 'LineWidth', 0.9);
    yline(ax, -eps_ref, '--', 'Color', [0.25 0.55 0.25], 'LineWidth', 0.9);
    grid(ax, 'on'); grid(ax, 'minor');
    xlim(ax, swin_n); ylim(ax, [-err_ylim, err_ylim]);
    if j == 1
        ylabel(ax, 'Diagonal \Delta\phi_s / max|\phi_s|', 'FontName', 'Times New Roman', 'FontSize', fs_axis);
    end
    set(ax, 'FontName', 'Times New Roman', 'FontSize', fs_tick, 'LineWidth', 0.8, 'Box', 'on');
    ax.GridAlpha = 0.20; ax.MinorGridAlpha = 0.10;
    text(ax, 0.97, 0.90, sprintf('max=%.1e', max(abs(e))), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'FontName', 'Times New Roman', 'FontSize', fs_note, 'BackgroundColor', 'w');
end
xlabel(tl2, 'x/\lambda_p', 'FontName', 'Times New Roman', 'FontSize', fs_axis);
fprintf('Normalized error reference max|phi_s| = %.3e\n', phi_ref_max);
fprintf('Machine epsilon reference (double) = %.3e\n', eps_ref);

out_err = fullfile(outDir, 'mf12_phi3_wavegroup_lines_error_pub.png');
exportgraphics(fig2, out_err, 'Resolution', 600);
fprintf('Saved: %s\n', out_err);

function phi = phi_surface_direct(order, coeffs, X, Y, t)
    [~, phi] = surfaceMF12_new(order, coeffs, X, Y, t);
end

function phi = phi_surface_spectral(coeffs, Lx, Ly, Nx, Ny, t)
    [~, phi] = surfaceMF12_spectral(coeffs, Lx, Ly, Nx, Ny, t);
end

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
