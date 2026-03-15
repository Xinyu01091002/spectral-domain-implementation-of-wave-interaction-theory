clc; clear; close all;
% Spectral-only visualization of a directional wave group.
% This script is intentionally not a direct-vs-spectral comparison.

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
run(fullfile(rootDir, 'setup_paths.m'));
outDir = fullfile(rootDir, 'outputs');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Wave-group style directional case
g = 9.81;
h = 150;
Ux = 0;
Uy = 0;
Lx =3000;
Ly = 3000;
Nx = 256;
Ny = 256;
t_eval = 24.0; % Default visualization at focus; use t_eval < 0 for pre-focus and t_eval > 0 for after-focus.

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
fprintf('Surface script: crossing angle = %.2f deg.\n', crossing_angle_deg);

energy_keep_frac = 0.999;
[W_sorted, idx_sort] = sort(W, 'descend');
cum_energy = cumsum(W_sorted);
N = find(cum_energy >= energy_keep_frac * cum_energy(end), 1, 'first');
idx = idx_sort(1:N);
kx = kx_all(idx);
ky = ky_all(idx);
fprintf('Surface script: retained %d components for %.4f%% cumulative energy.\n', N, 100 * energy_keep_frac);

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

% Spectral path only: build order-1 and order-3 fields on a regular grid.
c1 = mf12_spectral_coefficients(1, g, h, a, b, kx, ky, Ux, Uy);
c3 = mf12_spectral_coefficients(3, g, h, a, b, kx, ky, Ux, Uy);
c3.third_order_subharmonic_mode = 'skip';

[~, phis_linear, X, Y] = mf12_spectral_surface(c1, Lx, Ly, Nx, Ny, t_eval);
[~, phis_total, ~, ~] = mf12_spectral_surface(c3, Lx, Ly, Nx, Ny, t_eval);
phis_nonlinear = phis_total - phis_linear;

X_plot = (X - xf) / lambda_p;
Y_plot = (Y - yf) / lambda_p;

half_window_x = 0.48 * (Lx / lambda_p);
half_window_y = 0.48 * (Ly / lambda_p);

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 30 15]);
ax1 = axes(fig, 'Position', [0.11 0.10 0.34 0.82]);
ax2 = axes(fig, 'Position', [0.53 0.10 0.34 0.82]);

plot_surface_panel(ax1, X_plot, Y_plot, phis_linear, ...
    '', '$\phi_s^{(1)}$ (m$^2$/s)', half_window_x, half_window_y, 'centerline');
plot_surface_panel(ax2, X_plot, Y_plot, phis_nonlinear, ...
    '', '$\phi_s^{(20+22+33)}$ (m$^2$/s)', half_window_x, half_window_y, 'diagonal');

out_png = fullfile(outDir, 'mf12_phis_wavegroup_surface_sections.png');
out_png_print = fullfile(outDir, 'mf12_phis_wavegroup_surface_sections_print.png');
exportgraphics(fig, out_png, 'Resolution', 450);
saveas(fig, out_png_print);
fprintf('Saved: %s\n', out_png);
fprintf('Saved: %s\n', out_png_print);

function plot_surface_panel(ax, Xp, Yp, field, panel_title, z_label, half_window_x, half_window_y, section_mode)
    field_abs_max = max(abs(field(:)), [], 'omitnan');
    support = abs(field) >= 0.15 * max(field_abs_max, eps);
    if any(support(:))
        x_extent = max(abs(Xp(support)));
        y_extent = max(abs(Yp(support)));
        half_window_x_plot = min(half_window_x, max(1.8, 1.55 * x_extent));
        half_window_y_plot = min(half_window_y, max(1.8, 1.55 * y_extent));
    else
        half_window_x_plot = half_window_x;
        half_window_y_plot = half_window_y;
    end

    mask = abs(Xp) <= half_window_x_plot & abs(Yp) <= half_window_y_plot;
    field_view = field;
    field_view(~mask) = NaN;
    cla(ax);
    hold(ax, 'off');

    s = surfc(ax, Xp, Yp, field_view, 'EdgeColor', 'none');
    s(1).FaceColor = 'interp';
    s(1).FaceAlpha = 0.62;
    s(1).DisplayName = panel_title;
    s(2).HandleVisibility = 'off';
    s(2).LevelList = linspace(min(field_view(:), [], 'omitnan'), max(field_view(:), [], 'omitnan'), 6);
    s(2).EdgeColor = [0.18 0.18 0.18];
    s(2).LineWidth = 1.0;
    s(2).LineStyle = '-';

    view(ax, 30, 24);
    camlight(ax, 'left');
    lighting(ax, 'gouraud');
    material(ax, 'shiny');
    shading(ax, 'interp');
    grid(ax, 'on');
    grid(ax, 'minor');
    axis(ax, 'tight');
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;
    ax.LineWidth = 1.0;
    ax.TickDir = 'out';
    ax.TickLength = [0.015, 0.02];
    ax.GridAlpha = 0.35;
    ax.MinorGridAlpha = 0.18;
    ax.Title.FontSize = 15;
    ax.Title.FontWeight = 'normal';
    ax.Clipping = 'off';
    ax.LooseInset = [0.01 0.01 0.01 0.01];

    try
        colormap(ax, mymap('bone'));
    catch
        colormap(ax, parula);
    end

    xlabel(ax, '$(x-x_f)/\lambda_p$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(ax, '$(y-y_f)/\lambda_p$', 'Interpreter', 'latex', 'FontSize', 14);
    zlabel(ax, z_label, 'Interpreter', 'latex', 'FontSize', 14);
    if ~isempty(panel_title)
        title(ax, panel_title, 'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'normal');
    else
        title(ax, '');
    end

    z_min = min(field_view(:), [], 'omitnan');
    z_max = max(field_view(:), [], 'omitnan');
    z_range = max(z_max - z_min, eps);
    z_lower = z_min - 0.20 * z_range;
    z_upper = z_max + 0.20 * z_range;
    finite_vals = field_view(isfinite(field_view));
    c_abs = prctile(abs(finite_vals), 92);
    if ~(isfinite(c_abs) && c_abs > 0)
        c_abs = max(abs(finite_vals));
    end
    if ~(isfinite(c_abs) && c_abs > 0)
        c_abs = 1;
    end
    ax.XLim = [-half_window_x_plot, half_window_x_plot];
    ax.YLim = [-half_window_y_plot, half_window_y_plot];
    ax.ZLim = [z_lower, z_upper];
    clim(ax, [-c_abs, c_abs]);
    ax.XTick = symmetric_ticks(half_window_x_plot);
    ax.YTick = symmetric_ticks(half_window_y_plot);
    ax.PlotBoxAspectRatio = [2*half_window_x_plot, 2*half_window_y_plot, 0.55*max(2*half_window_x_plot, 2*half_window_y_plot)];
    ax.LabelFontSizeMultiplier = 1.0;
    hold(ax, 'on');

    switch lower(section_mode)
        case 'centerline'
            add_centerline_section(ax, Xp, Yp, field_view, z_lower, z_upper);
        case 'diagonal'
            add_diagonal_section(ax, Xp, Yp, field_view, z_lower, z_upper);
        otherwise
            error('Unknown section mode: %s', section_mode);
    end
end

function D = gaussian_spreading(theta, spread_angle_deg)
    theta_wrapped = angle(exp(1i * theta));
    sigma = deg2rad(spread_angle_deg);
    D = exp(-0.5 * (theta_wrapped / max(sigma, eps)).^2);
end

function ticks = symmetric_ticks(half_width)
    target_count = 7;
    raw_step = max(half_width, eps) / max((target_count - 1) / 2, 1);
    pow10 = 10^floor(log10(raw_step));
    candidates = [1, 2, 5, 10] * pow10;
    [~, idx] = min(abs(candidates - raw_step));
    step = candidates(idx);
    limit = step * floor(half_width / step);
    if limit <= 0
        ticks = [-half_width, 0, half_width];
    else
        ticks = -limit:step:limit;
    end
end

function add_centerline_section(ax, Xp, Yp, field_view, z_lower, z_upper)
    [~, iy0] = min(abs(Yp(:, 1)));
    xline = Xp(iy0, :);
    yline = Yp(iy0, :);
    zline = field_view(iy0, :);

    plot3(ax, xline, yline, zline, 'Color', [0 0 0], 'LineWidth', 2.4, 'HandleVisibility', 'off');

    xlim_now = get(ax, 'XLim');
    ylim_now = get(ax, 'YLim');
    zlim_now = get(ax, 'ZLim');
    x1 = xlim_now(1);
    x2 = xlim_now(2);
    zf1 = zlim_now(1);
    zf2 = zlim_now(2);
    y0 = yline(1);
    xs = [x1, x2, x2, x1, x1];
    ys = y0 * ones(size(xs));
    zs = [zf1, zf1, zf2, zf2, zf1];
    plot3(ax, xs, ys, zs, 'Color', [0 0 0], 'LineWidth', 0.9, 'HandleVisibility', 'off');

    xmid = x1 + 0.22 * (x2 - x1);
    xhalf = 0.24 * (x2 - x1);
    zmid = zf2 - 0.11 * (zf2 - zf1);
    zhalf = 0.09 * (zf2 - zf1);
    xcorn = [xmid - xhalf, xmid + xhalf; xmid - xhalf, xmid + xhalf];
    y_offset = 0.05 * diff(ylim_now);
    ycorn = [y0 + y_offset, y0 + y_offset; y0 + y_offset, y0 + y_offset];
    zcorn = [zmid + zhalf, zmid + zhalf; zmid - zhalf, zmid - zhalf];
    place_label_texture_on_plane(ax, xcorn, ycorn, zcorn, ...
        'Centerline ($y/\lambda_p = 0$)', [0 0 0]);
end

function add_diagonal_section(ax, Xp, Yp, field_view, z_lower, z_upper)
    nDiag = min(size(Xp, 1), size(Xp, 2));
    idx = 1:nDiag;
    xline = Xp(sub2ind(size(Xp), idx, idx));
    yline = Yp(sub2ind(size(Yp), idx, idx));
    zline = field_view(sub2ind(size(field_view), idx, idx));

    keep = isfinite(zline);
    xline = xline(keep);
    yline = yline(keep);
    zline = zline(keep);
    if isempty(xline)
        return;
    end

    plot3(ax, xline, yline, zline, 'Color', [0.85 0.15 0.15], 'LineWidth', 2.4, 'HandleVisibility', 'off');

    xlim_now = get(ax, 'XLim');
    ylim_now = get(ax, 'YLim');
    zlim_now = get(ax, 'ZLim');
    x1 = xlim_now(1);
    x2 = xlim_now(2);
    y1 = ylim_now(1);
    y2 = ylim_now(2);
    zf1 = zlim_now(1);
    zf2 = zlim_now(2);
    plot3(ax, [x1, x2, x2, x1, x1], [y1, y2, y2, y1, y1], ...
        [zf1, zf1, zf2, zf2, zf1], ...
        'Color', [0.85 0.15 0.15], 'LineWidth', 0.9, 'HandleVisibility', 'off');

    tmid = 0.72;
    xmid = x1 + tmid * (x2 - x1);
    ymid = y1 + tmid * (y2 - y1);
    dx = x2 - x1;
    dy = y2 - y1;
    nrm = hypot(dx, dy);
    if nrm <= eps
        return;
    end
    tx = dx / nrm;
    ty = dy / nrm;
    xhalf = 0.24 * nrm;
    nx = -ty;
    ny = tx;
    offset = 0.04 * max(diff(xlim_now), diff(ylim_now));
    zmid = zf2 - 0.11 * (zf2 - zf1);
    zhalf = 0.09 * (zf2 - zf1);
    xcorn = [xmid - xhalf*tx + offset*nx, xmid + xhalf*tx + offset*nx; ...
             xmid - xhalf*tx + offset*nx, xmid + xhalf*tx + offset*nx];
    ycorn = [ymid - xhalf*ty + offset*ny, ymid + xhalf*ty + offset*ny; ...
             ymid - xhalf*ty + offset*ny, ymid + xhalf*ty + offset*ny];
    zcorn = [zmid + zhalf, zmid + zhalf; zmid - zhalf, zmid - zhalf];
    place_label_texture_on_plane(ax, xcorn, ycorn, zcorn, ...
        'Diagonal line ($x=y$)', [0.85 0.15 0.15]);
end

function place_label_texture_on_plane(ax, xcorn, ycorn, zcorn, labelStr, labelColor)
    figTmp = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'pixels', 'Position', [100 100 1200 260]);
    axTmp = axes('Parent', figTmp, 'Position', [0 0 1 1], 'Color', 'w');
    axis(axTmp, 'off');
    xlim(axTmp, [0 1]);
    ylim(axTmp, [0 1]);

    text(axTmp, 0.5, 0.5, labelStr, ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 44, ...
        'FontWeight', 'bold', ...
        'Color', labelColor);

    drawnow;
    frame = getframe(axTmp);
    img = frame.cdata;
    close(figTmp);

    imgDouble = double(img) / 255;
    bgColor = squeeze(imgDouble(1,1,:)).';
    diffImg = sqrt(sum((imgDouble - reshape(bgColor,1,1,3)).^2, 3));
    alphaMask = ones(size(diffImg));
    alphaMask(diffImg < 0.06) = 0;
    rows = any(alphaMask > 0, 2);
    cols = any(alphaMask > 0, 1);
    if any(rows) && any(cols)
        r1 = find(rows, 1, 'first');
        r2 = find(rows, 1, 'last');
        c1 = find(cols, 1, 'first');
        c2 = find(cols, 1, 'last');
        pad = 8;
        r1 = max(1, r1 - pad); r2 = min(size(alphaMask, 1), r2 + pad);
        c1 = max(1, c1 - pad); c2 = min(size(alphaMask, 2), c2 + pad);
        alphaMask = alphaMask(r1:r2, c1:c2);
        img = img(r1:r2, c1:c2, :);
    end
    alphaMask = 0.95 * alphaMask;

    surface(ax, xcorn, ycorn, zcorn, ...
        'CData', img, ...
        'AlphaData', alphaMask, ...
        'FaceColor', 'texturemap', ...
        'FaceAlpha', 'texturemap', ...
        'EdgeColor', 'none', ...
        'FaceLighting', 'none', ...
        'HandleVisibility', 'off');
end
