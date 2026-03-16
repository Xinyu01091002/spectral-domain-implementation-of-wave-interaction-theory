clc; clear; close all;
scriptDir = fileparts(mfilename('fullpath'));
matlabDir = fileparts(scriptDir);
repoDir = fileparts(matlabDir);
run(fullfile(matlabDir, 'setup_paths.m'));
outDir = fullfile(repoDir, 'outputs');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Wave-group style directional case
g = 9.81;
h = 150;
Ux = 0;
Uy = 0;
Lx = 3000;
Ly = 3000;
Nx = 256;
Ny = 256;
t_eval = 24.0;

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

kp = 0.0279;
h = 1 / kp;
lambda_p = 2*pi / kp;
omega_p = sqrt(g*kp*tanh(h*kp));
Tp = 2*pi / omega_p;
t_focus = 0.0; % By construction, both wave groups overlap at the domain center when t_eval = 0.

Akp = 0.35;
kw_left = 0.004606;
kw_right = kw_left;
A_focus_total = Akp / kp;

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

crossing_angle_deg = abs(angle(exp(1i * deg2rad(group1_heading_deg - group2_heading_deg)))) * 180 / pi;
fprintf('Line script: crossing angle = %.2f deg.\n', crossing_angle_deg);
energy_keep_frac = 0.99;
max_components = 150; % Safety cap for the direct third-order path.
[W_sorted, idx_sort] = sort(W, 'descend');
cum_energy = cumsum(W_sorted);
N_energy = find(cum_energy >= energy_keep_frac * cum_energy(end), 1, 'first');
N = min(N_energy, max_components);
idx = idx_sort(1:N);
kx = kx_all(idx);
ky = ky_all(idx);
fprintf('Line script: retained %d components for %.4f%% cumulative energy', N, 100 * energy_keep_frac);
if N < N_energy
    fprintf(' (capped from %d for direct third-order cost)', N_energy);
end
fprintf('.\n');

% Amplitudes and focusing phase at domain center
Wsel = W(idx);
amp = Wsel;
amp = amp * (A_focus_total / max(sum(amp), eps));
xf = Lx/2;
yf = Ly/2;
omega_lin = sqrt(g*hypot(kx,ky).*tanh(h*hypot(kx,ky)));
phase = -(kx*xf + ky*yf) + omega_lin*t_focus;
a = amp .* cos(phase);
b = amp .* sin(phase);

% Grid
x = (0:Nx-1) * (Lx/Nx);
y = (0:Ny-1) * (Ly/Ny);
[X, Y] = meshgrid(x, y);

% Keep the original benchmark pairing:
%   direct   = mf12_direct_coefficients + mf12_direct_surface
%   spectral = mf12_spectral_coefficients + mf12_spectral_surface
c1_d = mf12_direct_coefficients(1, g, h, a, b, kx, ky, Ux, Uy);
c2_d = mf12_direct_coefficients(2, g, h, a, b, kx, ky, Ux, Uy);
c3_d = mf12_direct_coefficients(3, g, h, a, b, kx, ky, Ux, Uy);
c1_s = mf12_spectral_coefficients(1, g, h, a, b, kx, ky, Ux, Uy);
c2_s = mf12_spectral_coefficients(2, g, h, a, b, kx, ky, Ux, Uy);
c3_s = mf12_spectral_coefficients(3, g, h, a, b, kx, ky, Ux, Uy);
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
phi1_d = phi_surface_direct(1, c1_d, X, Y, t_eval);
phi2_d_total = phi_surface_direct(2, c2_d, X, Y, t_eval);
phi2_d_inc = phi2_d_total - phi1_d;
phi2sup_d = split_by_wavenumber(phi2_d_inc, Lx, Ly, k2_cut, 'high');
phi2sub_d = split_by_wavenumber(phi2_d_inc, Lx, Ly, k2_cut, 'low');
phi3_d_total = phi_surface_direct(3, c3_d_sup, X, Y, t_eval);
phi3_d = phi3_d_total - phi2_d_total;

% Decomposed components (Spectral)
phi1_s = phi_surface_spectral(c1_s, Lx, Ly, Nx, Ny, t_eval);
phi2_s_total = phi_surface_spectral(c2_s, Lx, Ly, Nx, Ny, t_eval);
phi2_s_inc = phi2_s_total - phi1_s;
phi2sup_s = split_by_wavenumber(phi2_s_inc, Lx, Ly, k2_cut, 'high');
phi2sub_s = split_by_wavenumber(phi2_s_inc, Lx, Ly, k2_cut, 'low');
phi3_s_total = phi_surface_spectral(c3_s, Lx, Ly, Nx, Ny, t_eval);
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

% Plot window: center on the envelope peak of the computed response.
phi_total_center_d = center_d{1} + center_d{2} + center_d{3} + center_d{4};
phi_total_center_s = center_s{1} + center_s{2} + center_s{3} + center_s{4};
phi_total_diag_d = diag_d{1} + diag_d{2} + diag_d{3} + diag_d{4};
phi_total_diag_s = diag_s{1} + diag_s{2} + diag_s{3} + diag_s{4};

[~, ix_peak_d] = max(local_envelope(phi_total_center_d));
[~, ix_peak_s] = max(local_envelope(phi_total_center_s));
[~, is_peak_d] = max(local_envelope(phi_total_diag_d));
[~, is_peak_s] = max(local_envelope(phi_total_diag_s));

x_center_plot = 0.5 * (x(ix_peak_d) + x(ix_peak_s));
s_center_plot = 0.5 * (s(is_peak_d) + s(is_peak_s));

half_window_lambda = 3.5;  % show a tighter window around the dominant envelope peak
xn = (x - x_center_plot) / lambda_p;
sn = (s - s_center_plot) / lambda_p;
xwin_n = [-half_window_lambda, half_window_lambda];
swin_n = [-half_window_lambda, half_window_lambda];

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
fs_title = 16;
fs_axis = 14;
fs_tick = 13;
fs_legend = 13;
fs_note = 11;

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
fig1 = figure('Color','w','Units','centimeters','Position',[2 2 28 11]);
ha1 = tight_subplot_local(2, 4, [0.01 0.035], [0.12 0.12], [0.10 0.03]);
h_leg_direct = [];
h_leg_spec = [];

for j = 1:4
    ax = ha1(j);
    h1 = plot(ax, xn, center_d{j}, '-', 'Color', c_direct, 'LineWidth', lw_main); hold(ax, 'on');
    h2 = plot(ax, xn, center_s{j}, '--', 'Color', c_spec, 'LineWidth', lw_main);
    grid(ax, 'on'); grid(ax, 'minor');
    xlim(ax, xwin_n);
    title(ax, labels{j}, 'FontName', 'Times New Roman', 'FontSize', fs_title, 'FontWeight', 'normal');
    if j == 1
        ylabel(ax, 'Centerline \phi_s', 'FontName', 'Times New Roman', 'FontSize', fs_axis);
        h_leg_direct = h1;
        h_leg_spec = h2;
        yl = ylim(ax); d = max(eps, yl(2)-yl(1)); ylim(ax, [yl(1)-0.04*d, yl(2)+0.08*d]);
    else
        yl = ylim(ax); d = max(eps, yl(2)-yl(1)); ylim(ax, [yl(1)-0.03*d, yl(2)+0.06*d]);
    end
    set(ax, 'FontName', 'Times New Roman', 'FontSize', fs_tick, 'LineWidth', 0.8, 'Box', 'on');
    ax.GridAlpha = 0.20; ax.MinorGridAlpha = 0.10;
    ax.XTickLabel = [];
end

for j = 1:4
    ax = ha1(4 + j);
    h1 = plot(ax, sn, diag_d{j}, '-', 'Color', c_direct, 'LineWidth', lw_main); hold(ax, 'on');
    h2 = plot(ax, sn, diag_s{j}, '--', 'Color', c_spec, 'LineWidth', lw_main);
    grid(ax, 'on'); grid(ax, 'minor');
    xlim(ax, swin_n);
    if j == 1
        ylabel(ax, 'Diagonal \phi_s', 'FontName', 'Times New Roman', 'FontSize', fs_axis);
        yl = ylim(ax); d = max(eps, yl(2)-yl(1)); ylim(ax, [yl(1)-0.04*d, yl(2)+0.08*d]);
    else
        yl = ylim(ax); d = max(eps, yl(2)-yl(1)); ylim(ax, [yl(1)-0.03*d, yl(2)+0.06*d]);
    end
    set(ax, 'FontName', 'Times New Roman', 'FontSize', fs_tick, 'LineWidth', 0.8, 'Box', 'on');
    ax.GridAlpha = 0.20; ax.MinorGridAlpha = 0.10;
end
lgd1 = legend(ha1(1), [h_leg_direct, h_leg_spec], {'Direct', 'Spectral'}, ...
    'Orientation', 'horizontal', 'Location', 'northoutside', ...
    'FontName', 'Times New Roman', 'FontSize', fs_legend, 'Box', 'off');
set(lgd1, 'Units', 'normalized');
lgd1.Position(1) = 0.39;
lgd1.Position(2) = 0.94;
lgd1.Position(3) = 0.22;
annotation(fig1, 'textbox', [0.43 0.015 0.14 0.035], 'String', 'x/\lambda_p', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontName', 'Times New Roman', 'FontSize', fs_axis);

out_cmp = fullfile(outDir, 'mf12_phi3_wavegroup_lines_comparison_pub.png');
exportgraphics(fig1, out_cmp, 'Resolution', 600);
fprintf('Saved: %s\n', out_cmp);

% -------- Figure 2: Error-only panels --------
fig2 = figure('Color','w','Units','centimeters','Position',[2 2 27 9.6]);
ha2 = tight_subplot_local(2, 4, [0.045 0.028], [0.11 0.08], [0.10 0.03]);

for j = 1:4
    ax = ha2(j);
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
    ax = ha2(4 + j);
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
annotation(fig2, 'textbox', [0.43 0.015 0.14 0.035], 'String', 'x/\lambda_p', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontName', 'Times New Roman', 'FontSize', fs_axis);
fprintf('Normalized error reference max|phi_s| = %.3e\n', phi_ref_max);
fprintf('Machine epsilon reference (double) = %.3e\n', eps_ref);

out_err = fullfile(outDir, 'mf12_phi3_wavegroup_lines_error_pub.png');
exportgraphics(fig2, out_err, 'Resolution', 600);
fprintf('Saved: %s\n', out_err);

function phi = phi_surface_direct(order, coeffs, X, Y, t)
    [~, phi] = mf12_direct_surface(order, coeffs, X, Y, t);
end

function phi = phi_surface_spectral(coeffs, Lx, Ly, Nx, Ny, t)
    [~, phi] = mf12_spectral_surface(coeffs, Lx, Ly, Nx, Ny, t);
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

function D = gaussian_spreading(theta, spread_angle_deg)
    theta_wrapped = angle(exp(1i * theta));
    sigma = deg2rad(spread_angle_deg);
    D = exp(-0.5 * (theta_wrapped / max(sigma, eps)).^2);
end

function env = local_envelope(sig)
    sig = sig(:).';
    valid = isfinite(sig);
    if ~any(valid)
        env = zeros(size(sig));
        return;
    end

    if any(~valid)
        idx = 1:numel(sig);
        sig(~valid) = interp1(idx(valid), sig(valid), idx(~valid), 'linear', 'extrap');
    end

    env = abs(hilbert(sig));
end

function ha = tight_subplot_local(Nh, Nw, gap, marg_h, marg_w)
    if nargin < 3 || isempty(gap)
        gap = .02;
    end
    if nargin < 4 || isempty(marg_h)
        marg_h = .05;
    end
    if nargin < 5 || isempty(marg_w)
        marg_w = .05;
    end

    if numel(gap) == 1
        gap = [gap gap];
    end
    if numel(marg_h) == 1
        marg_h = [marg_h marg_h];
    end
    if numel(marg_w) == 1
        marg_w = [marg_w marg_w];
    end

    axh = (1 - sum(marg_h) - (Nh - 1) * gap(1)) / Nh;
    axw = (1 - sum(marg_w) - (Nw - 1) * gap(2)) / Nw;
    py = 1 - marg_h(2) - axh;

    ha = gobjects(Nh * Nw, 1);
    idx = 0;
    for ih = 1:Nh
        px = marg_w(1);
        for ix = 1:Nw
            idx = idx + 1;
            ha(idx) = axes('Units', 'normalized', 'Position', [px py axw axh]);
            px = px + axw + gap(2);
        end
        py = py - axh - gap(1);
    end
end
