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
Lx = 3000;
Ly = 3000;
Nx = 64;
Ny = 64;
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
kp = (2*pi/Tp)^2 / g;
lambda_p = 2*pi / kp;
% Crossing-sea wave-group setup: two directional groups with different headings.
theta1 = deg2rad(25);    % group-1 mean direction
theta2 = deg2rad(-35);   % group-2 mean direction (crossing)
sig_k = 0.10*kp;         % narrow-band => wave-group behavior
sig_t1 = deg2rad(8);
sig_t2 = deg2rad(10);
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
A0 = 0.2;
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

% Order-wise decomposition coefficients.
% first harmonic -> order 1 coefficients
% second harmonics -> order 2 coefficients
% third superharmonic -> order 3 coefficients
c1 = coeffsMF12(1, g, h, a, b, kx, ky, Ux, Uy);
c2 = coeffsMF12(2, g, h, a, b, kx, ky, Ux, Uy);
c3 = coeffsMF12(3, g, h, a, b, kx, ky, Ux, Uy);

% Decomposed components (Direct)
phi1_d = phi_component_direct('first', c1, X, Y, t);
phi2sup_d = phi_component_direct('second_super', c2, X, Y, t);
phi2sub_d = phi_component_direct('second_sub', c2, X, Y, t);
phi3_d = phi_component_direct('third_term', c3, X, Y, t);

% Decomposed components (Spectral)
phi1_s = phi_component_spectral('first', c1, Lx, Ly, Nx, Ny, t);
phi2sup_s = phi_component_spectral('second_super', c2, Lx, Ly, Nx, Ny, t);
phi2sub_s = phi_component_spectral('second_sub', c2, Lx, Ly, Nx, Ny, t);
phi3_s = phi_component_spectral('third_term', c3, Lx, Ly, Nx, Ny, t);

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
    if j <= 4, ax.XTickLabel = []; end
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
    if j <= 4, ax.XTickLabel = []; end
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

function phi = phi_component_direct(component, coeffs, X, Y, t)
    phi = zeros(size(X));
    N = coeffs.N;

    switch component
        case 'first'
            for n = 1:N
                th = coeffs.omega(n)*t - coeffs.kx(n)*X - coeffs.ky(n)*Y;
                phi = phi + (coeffs.mu(n)+coeffs.muStar(n)) * ...
                    (coeffs.a(n)*sin(th) - coeffs.b(n)*cos(th));
            end

        case 'second_super'
            % Self-self (superharmonic)
            for n = 1:N
                th = coeffs.omega(n)*t - coeffs.kx(n)*X - coeffs.ky(n)*Y;
                phi = phi + coeffs.mu_2(n) * ...
                    (coeffs.A_2(n)*sin(2*th) - coeffs.B_2(n)*cos(2*th));
            end
            % Pair sum pm=+1
            cnm = 0;
            for n = 1:N
                for m = n+1:N
                    for pm = [1 -1]
                        cnm = cnm + 1;
                        if pm ~= 1, continue; end
                        thn = coeffs.omega(n)*t - coeffs.kx(n)*X - coeffs.ky(n)*Y;
                        thm = coeffs.omega(m)*t - coeffs.kx(m)*X - coeffs.ky(m)*Y;
                        th = thn + pm*thm;
                        phi = phi + coeffs.mu_npm(cnm) * ...
                            (coeffs.A_npm(cnm)*sin(th) - coeffs.B_npm(cnm)*cos(th));
                    end
                end
            end

        case 'second_sub'
            % Pair difference pm=-1 only
            cnm = 0;
            for n = 1:N
                for m = n+1:N
                    for pm = [1 -1]
                        cnm = cnm + 1;
                        if pm ~= -1, continue; end
                        thn = coeffs.omega(n)*t - coeffs.kx(n)*X - coeffs.ky(n)*Y;
                        thm = coeffs.omega(m)*t - coeffs.kx(m)*X - coeffs.ky(m)*Y;
                        th = thn + pm*thm;
                        phi = phi + coeffs.mu_npm(cnm) * ...
                            (coeffs.A_npm(cnm)*sin(th) - coeffs.B_npm(cnm)*cos(th));
                    end
                end
            end

        case 'third_term'
            % Self 3rd harmonic
            for n = 1:N
                th = coeffs.omega(n)*t - coeffs.kx(n)*X - coeffs.ky(n)*Y;
                phi = phi + coeffs.mu_3(n) * ...
                    (coeffs.A_3(n)*sin(3*th) - coeffs.B_3(n)*cos(3*th));
            end
            % Double sums (pm=+1 only)
            cnm = 0;
            for n = 1:N
                for m = n+1:N
                    for pm = [1 -1]
                        cnm = cnm + 1;
                        if pm ~= 1, continue; end
                        thn = coeffs.omega(n)*t - coeffs.kx(n)*X - coeffs.ky(n)*Y;
                        thm = coeffs.omega(m)*t - coeffs.kx(m)*X - coeffs.ky(m)*Y;
                        th1 = thn + pm*2*thm;
                        th2 = 2*thn + pm*thm;
                        phi = phi + coeffs.mu_np2m(cnm) * ...
                            (coeffs.A_np2m(cnm)*sin(th1) - coeffs.B_np2m(cnm)*cos(th1));
                        phi = phi + coeffs.mu_2npm(cnm) * ...
                            (coeffs.A_2npm(cnm)*sin(th2) - coeffs.B_2npm(cnm)*cos(th2));
                    end
                end
            end
            % Triple sums (pmm=1,pmp=1 only, with factor 2)
            c3 = 0;
            for n = 1:N
                for m = n+1:N
                    for pmm = [1 -1]
                        for p = m+1:N
                            for pmp = [1 -1]
                                c3 = c3 + 1;
                                if ~(pmm == 1 && pmp == 1), continue; end
                                thn = coeffs.omega(n)*t - coeffs.kx(n)*X - coeffs.ky(n)*Y;
                                thm = coeffs.omega(m)*t - coeffs.kx(m)*X - coeffs.ky(m)*Y;
                                thp = coeffs.omega(p)*t - coeffs.kx(p)*X - coeffs.ky(p)*Y;
                                th = thn + pmm*thm + pmp*thp;
                                phi = phi + 2*coeffs.mu_npmpp(c3) * ...
                                    (coeffs.A_npmpp(c3)*sin(th) - coeffs.B_npmpp(c3)*cos(th));
                            end
                        end
                    end
                end
            end

        otherwise
            error('Unknown component: %s', component);
    end
end

function phi = phi_component_spectral(component, coeffs, Lx, Ly, Nx, Ny, t)
    dkx = 2*pi/Lx;
    dky = 2*pi/Ly;
    spec_phi = complex(zeros(Ny, Nx));
    N = coeffs.N;

    switch component
        case 'first'
            Z = (coeffs.a(:) + 1i*coeffs.b(:)) .* exp(-1i * coeffs.omega(:) * t);
            mu = coeffs.mu(:) + coeffs.muStar(:);
            spec_phi = deposit(spec_phi, coeffs.kx(:), coeffs.ky(:), Z .* mu * (1i), dkx, dky, Nx, Ny);

        case 'second_super'
            Z2 = (coeffs.A_2(:) + 1i*coeffs.B_2(:)) .* exp(-1i * (2*coeffs.omega(:)) * t);
            spec_phi = deposit(spec_phi, 2*coeffs.kx(:), 2*coeffs.ky(:), Z2 .* coeffs.mu_2(:) * (1i), dkx, dky, Nx, Ny);
            cnm = 0;
            for n = 1:N
                for m = n+1:N
                    for pm = [1 -1]
                        cnm = cnm + 1;
                        if pm ~= 1, continue; end
                        kx_eff = coeffs.kx(n) + pm*coeffs.kx(m);
                        ky_eff = coeffs.ky(n) + pm*coeffs.ky(m);
                        om_eff = coeffs.omega(n) + pm*coeffs.omega(m); % includes dispersion correction
                        Z = (coeffs.A_npm(cnm) + 1i*coeffs.B_npm(cnm)) * exp(-1i * om_eff * t);
                        spec_phi = deposit(spec_phi, kx_eff, ky_eff, Z * coeffs.mu_npm(cnm) * (1i), dkx, dky, Nx, Ny);
                    end
                end
            end

        case 'second_sub'
            cnm = 0;
            for n = 1:N
                for m = n+1:N
                    for pm = [1 -1]
                        cnm = cnm + 1;
                        if pm ~= -1, continue; end
                        kx_eff = coeffs.kx(n) + pm*coeffs.kx(m);
                        ky_eff = coeffs.ky(n) + pm*coeffs.ky(m);
                        om_eff = coeffs.omega(n) + pm*coeffs.omega(m); % includes dispersion correction
                        Z = (coeffs.A_npm(cnm) + 1i*coeffs.B_npm(cnm)) * exp(-1i * om_eff * t);
                        spec_phi = deposit(spec_phi, kx_eff, ky_eff, Z * coeffs.mu_npm(cnm) * (1i), dkx, dky, Nx, Ny);
                    end
                end
            end

        case 'third_term'
            Z3 = (coeffs.A_3(:) + 1i*coeffs.B_3(:)) .* exp(-1i * (3*coeffs.omega(:)) * t);
            spec_phi = deposit(spec_phi, 3*coeffs.kx(:), 3*coeffs.ky(:), Z3 .* coeffs.mu_3(:) * (1i), dkx, dky, Nx, Ny);
            % Double sums: keep superharmonic pm=+1 terms only, with corrected omega combinations.
            cnm = 0;
            for n = 1:N
                for m = n+1:N
                    for pm = [1 -1]
                        cnm = cnm + 1;
                        if pm ~= 1, continue; end

                        om_np2m = coeffs.omega(n) + 2*coeffs.omega(m);
                        kx_np2m = coeffs.kx(n) + 2*coeffs.kx(m);
                        ky_np2m = coeffs.ky(n) + 2*coeffs.ky(m);
                        Znp2m = (coeffs.A_np2m(cnm) + 1i*coeffs.B_np2m(cnm)) * exp(-1i * om_np2m * t);
                        spec_phi = deposit(spec_phi, kx_np2m, ky_np2m, Znp2m * coeffs.mu_np2m(cnm) * (1i), dkx, dky, Nx, Ny);

                        om_2npm = 2*coeffs.omega(n) + coeffs.omega(m);
                        kx_2npm = 2*coeffs.kx(n) + coeffs.kx(m);
                        ky_2npm = 2*coeffs.ky(n) + coeffs.ky(m);
                        Z2npm = (coeffs.A_2npm(cnm) + 1i*coeffs.B_2npm(cnm)) * exp(-1i * om_2npm * t);
                        spec_phi = deposit(spec_phi, kx_2npm, ky_2npm, Z2npm * coeffs.mu_2npm(cnm) * (1i), dkx, dky, Nx, Ny);
                    end
                end
            end

            % Triple sums: superharmonic pmm=+1, pmp=+1 only, with corrected omega combinations.
            c3 = 0;
            for n = 1:N
                for m = n+1:N
                    for pmm = [1 -1]
                        for p = m+1:N
                            for pmp = [1 -1]
                                c3 = c3 + 1;
                                if ~(pmm == 1 && pmp == 1), continue; end
                                om = coeffs.omega(n) + coeffs.omega(m) + coeffs.omega(p);
                                kx_eff = coeffs.kx(n) + coeffs.kx(m) + coeffs.kx(p);
                                ky_eff = coeffs.ky(n) + coeffs.ky(m) + coeffs.ky(p);
                                Z = 2*(coeffs.A_npmpp(c3) + 1i*coeffs.B_npmpp(c3)) * exp(-1i * om * t);
                                spec_phi = deposit(spec_phi, kx_eff, ky_eff, Z * coeffs.mu_npmpp(c3) * (1i), dkx, dky, Nx, Ny);
                            end
                        end
                    end
                end
            end

        otherwise
            error('Unknown component: %s', component);
    end

    phi = real(ifft2(spec_phi)) * (Nx * Ny);
end

function spec = deposit(spec, kx_in, ky_in, values, dkx, dky, Nx, Ny)
    ux = (kx_in(:) / dkx);
    uy = (ky_in(:) / dky);
    vals = values(:);

    valid = isfinite(ux) & isfinite(uy) & isfinite(vals);
    ux = ux(valid); uy = uy(valid); vals = vals(valid);
    if isempty(vals), return; end

    ix0 = floor(ux); iy0 = floor(uy);
    fx = ux - ix0; fy = uy - iy0;
    tol = 1e-12;
    fx(abs(fx) < tol) = 0; fy(abs(fy) < tol) = 0;
    fx(abs(fx-1) < tol) = 1; fy(abs(fy-1) < tol) = 1;

    ix1 = ix0 + 1; iy1 = iy0 + 1;
    idx_x00 = mod(ix0, Nx) + 1; idx_y00 = mod(iy0, Ny) + 1;
    idx_x10 = mod(ix1, Nx) + 1; idx_y10 = idx_y00;
    idx_x01 = idx_x00;          idx_y01 = mod(iy1, Ny) + 1;
    idx_x11 = idx_x10;          idx_y11 = idx_y01;

    w00 = (1-fx).*(1-fy);
    w10 = fx.*(1-fy);
    w01 = (1-fx).*fy;
    w11 = fx.*fy;

    idx_x = [idx_x00; idx_x10; idx_x01; idx_x11];
    idx_y = [idx_y00; idx_y10; idx_y01; idx_y11];
    vals4 = [vals.*w00; vals.*w10; vals.*w01; vals.*w11];
    nz = abs(vals4) > 0;
    idx_x = idx_x(nz); idx_y = idx_y(nz); vals4 = vals4(nz);
    if isempty(vals4), return; end

    addv = accumarray([idx_y, idx_x], vals4, [Ny, Nx]);
    spec = spec + addv;
end
