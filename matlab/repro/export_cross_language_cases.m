clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
matlabDir = fileparts(scriptDir);
repoDir = fileparts(matlabDir);
run(fullfile(matlabDir, 'setup_paths.m'));

casesDir = fullfile(repoDir, 'cross_language_comparison', 'cases');
if ~exist(casesDir, 'dir'), mkdir(casesDir); end

caseFilter = strtrim(getenv('MF12_EXPORT_CASE'));

maybe_export('minimal_small', @() export_minimal_small(casesDir), caseFilter);
maybe_export('wavegroup_regression', @() export_wavegroup_regression(casesDir), caseFilter);
maybe_export('benchmark_medium', @() export_benchmark_medium(casesDir), caseFilter);
maybe_export('benchmark_dense_100', @() export_benchmark_dense_100(casesDir), caseFilter);
maybe_export('benchmark_dense_300', @() export_benchmark_dense_300(casesDir), caseFilter);
maybe_export('benchmark_dense_600', @() export_benchmark_dense_600(casesDir), caseFilter);

fprintf('Exported shared cross-language cases under: %s\n', casesDir);

function export_minimal_small(casesDir)
cfg = struct();
cfg.case_id = 'minimal_small';
cfg.description = 'Small deterministic sanity case based on the minimal spectral example.';
cfg.purpose = 'sanity';
cfg.order = 2;
cfg.g = 9.81;
cfg.h = 40;
cfg.z_kinematics = -0.2 * cfg.h;
cfg.Ux = 0;
cfg.Uy = 0;
cfg.Lx = 400;
cfg.Ly = 400;
cfg.Nx = 64;
cfg.Ny = 64;
cfg.t = 0;
cfg.subharmonic_mode = 'skip';
cfg.a = [0.18, 0.10, 0.08];
cfg.b = [0.00, 0.04, -0.02];
cfg.kx = [0.080, 0.075, 0.070];
cfg.ky = [0.000, 0.012, -0.010];
cfg.tolerances = struct('eta_max_abs_err', 1e-9, 'phi_max_abs_err', 1e-9, 'eta_rms_err', 1e-10, 'phi_rms_err', 1e-10, 'eta_relative_l2_err', 1e-10, 'phi_relative_l2_err', 1e-10);
export_case(casesDir, cfg);
end

function export_wavegroup_regression(casesDir)
cfg = struct();
cfg.case_id = 'wavegroup_regression';
cfg.description = 'Representative directional wave-group regression case aligned with matlab/tests/regression_wavegroup_phi3.m.';
cfg.purpose = 'regression';
cfg.order = 3;
cfg.g = 9.81;
Tp = 12;
cfg.t = 5 * Tp;
kp = (2*pi/Tp)^2 / cfg.g;
cfg.h = 1 / kp;
cfg.z_kinematics = -0.2 * cfg.h;
cfg.Ux = 0;
cfg.Uy = 0;
cfg.Lx = 2000;
cfg.Ly = 2000;
cfg.Nx = 64;
cfg.Ny = 64;
cfg.subharmonic_mode = 'skip';

rng(1234);
dkx = 2*pi/cfg.Lx;
dky = 2*pi/cfg.Ly;
kx_idx_all = (-floor(cfg.Nx/2)):(ceil(cfg.Nx/2)-1);
ky_idx_all = (-floor(cfg.Ny/2)):(ceil(cfg.Ny/2)-1);
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
idx = idx_sort(1:40);
kx = kx_all(idx);
ky = ky_all(idx);
A0 = 0.25;
amp = A0 * W(idx) / max(W(idx));
xf = cfg.Lx/2;
yf = cfg.Ly/2;
omega_lin = sqrt(cfg.g*hypot(kx,ky).*tanh(cfg.h*hypot(kx,ky)));
phase = -(kx*xf + ky*yf) + omega_lin*cfg.t;
cfg.a = amp .* cos(phase);
cfg.b = amp .* sin(phase);
cfg.kx = kx;
cfg.ky = ky;
cfg.tolerances = struct('eta_max_abs_err', 1e-6, 'phi_max_abs_err', 1e-6, 'eta_rms_err', 1e-7, 'phi_rms_err', 1e-7, 'eta_relative_l2_err', 1e-6, 'phi_relative_l2_err', 1e-6);
export_case(casesDir, cfg);
end

function export_benchmark_medium(casesDir)
cfg = struct();
cfg.case_id = 'benchmark_medium';
cfg.description = 'Moderate deterministic benchmark case aligned with the direct-vs-spectral runtime study.';
cfg.purpose = 'benchmark';
cfg.order = 3;
cfg.g = 9.81;
cfg.h = 150;
cfg.z_kinematics = -0.2 * cfg.h;
cfg.Ux = 0;
cfg.Uy = 0;
cfg.t = 0.0;
cfg.Lx = 3000;
cfg.Ly = 3000;
cfg.Nx = 64;
cfg.Ny = 64;
cfg.subharmonic_mode = 'skip';
cfg.Tp = 12;
cfg.theta1 = deg2rad(25);
cfg.theta2 = deg2rad(-35);
cfg.sig_k_factor = 0.12;
cfg.sig_t1 = deg2rad(10);
cfg.sig_t2 = deg2rad(12);
cfg.w1 = 0.55;
cfg.w2 = 0.45;
cfg.Hs_target = 3.5;
cfg.component_count = 40;

cfg = populate_directional_benchmark_components(cfg);
cfg.tolerances = struct('eta_max_abs_err', 1e-6, 'phi_max_abs_err', 1e-6, 'eta_rms_err', 1e-7, 'phi_rms_err', 1e-7, 'eta_relative_l2_err', 1e-6, 'phi_relative_l2_err', 1e-6);
export_case(casesDir, cfg);
end

function export_benchmark_dense_300(casesDir)
cfg = struct();
cfg.case_id = 'benchmark_dense_300';
cfg.description = 'Dense deterministic benchmark case for cross-language speed sweeps with up to 300 retained components.';
cfg.purpose = 'benchmark';
cfg.order = 3;
cfg.g = 9.81;
cfg.h = 150;
cfg.z_kinematics = -0.2 * cfg.h;
cfg.Ux = 0;
cfg.Uy = 0;
cfg.t = 0.0;
cfg.Lx = 3000;
cfg.Ly = 3000;
cfg.Nx = 64;
cfg.Ny = 64;
cfg.subharmonic_mode = 'skip';
cfg.Tp = 12;
cfg.theta1 = deg2rad(25);
cfg.theta2 = deg2rad(-35);
cfg.sig_k_factor = 0.12;
cfg.sig_t1 = deg2rad(10);
cfg.sig_t2 = deg2rad(12);
cfg.w1 = 0.55;
cfg.w2 = 0.45;
cfg.Hs_target = 3.5;
cfg.component_count = 300;

cfg = populate_directional_benchmark_components(cfg);
cfg.tolerances = struct('eta_max_abs_err', 1e-6, 'phi_max_abs_err', 1e-6, 'eta_rms_err', 1e-7, 'phi_rms_err', 1e-7, 'eta_relative_l2_err', 1e-6, 'phi_relative_l2_err', 1e-6);
export_case(casesDir, cfg);
end

function export_benchmark_dense_100(casesDir)
cfg = struct();
cfg.case_id = 'benchmark_dense_100';
cfg.description = 'Dense deterministic benchmark case for cross-language speed sweeps with up to 100 retained components.';
cfg.purpose = 'benchmark';
cfg.order = 3;
cfg.g = 9.81;
cfg.h = 150;
cfg.z_kinematics = -0.2 * cfg.h;
cfg.Ux = 0;
cfg.Uy = 0;
cfg.t = 0.0;
cfg.Lx = 3000;
cfg.Ly = 3000;
cfg.Nx = 64;
cfg.Ny = 64;
cfg.subharmonic_mode = 'skip';
cfg.Tp = 12;
cfg.theta1 = deg2rad(25);
cfg.theta2 = deg2rad(-35);
cfg.sig_k_factor = 0.12;
cfg.sig_t1 = deg2rad(10);
cfg.sig_t2 = deg2rad(12);
cfg.w1 = 0.55;
cfg.w2 = 0.45;
cfg.Hs_target = 3.5;
cfg.component_count = 100;

cfg = populate_directional_benchmark_components(cfg);
cfg.tolerances = struct('eta_max_abs_err', 1e-6, 'phi_max_abs_err', 1e-6, 'eta_rms_err', 1e-7, 'phi_rms_err', 1e-7, 'eta_relative_l2_err', 1e-6, 'phi_relative_l2_err', 1e-6);
export_case(casesDir, cfg);
end

function export_benchmark_dense_600(casesDir)
cfg = struct();
cfg.case_id = 'benchmark_dense_600';
cfg.description = 'Dense deterministic benchmark case for cross-language speed sweeps with up to 600 retained components.';
cfg.purpose = 'benchmark';
cfg.order = 3;
cfg.g = 9.81;
cfg.h = 150;
cfg.z_kinematics = -0.2 * cfg.h;
cfg.Ux = 0;
cfg.Uy = 0;
cfg.t = 0.0;
cfg.Lx = 3000;
cfg.Ly = 3000;
cfg.Nx = 64;
cfg.Ny = 64;
cfg.subharmonic_mode = 'skip';
cfg.Tp = 12;
cfg.theta1 = deg2rad(25);
cfg.theta2 = deg2rad(-35);
cfg.sig_k_factor = 0.12;
cfg.sig_t1 = deg2rad(10);
cfg.sig_t2 = deg2rad(12);
cfg.w1 = 0.55;
cfg.w2 = 0.45;
cfg.Hs_target = 3.5;
cfg.component_count = 600;
cfg.export_reference = false;

cfg = populate_directional_benchmark_components(cfg);
cfg.tolerances = struct('eta_max_abs_err', 1e-6, 'phi_max_abs_err', 1e-6, 'eta_rms_err', 1e-7, 'phi_rms_err', 1e-7, 'eta_relative_l2_err', 1e-6, 'phi_relative_l2_err', 1e-6);
export_case(casesDir, cfg);
end

function cfg = populate_directional_benchmark_components(cfg)
rng(20260306);
dkx = 2*pi/cfg.Lx;
dky = 2*pi/cfg.Ly;
kx_idx_all = (-floor(cfg.Nx/2)):(ceil(cfg.Nx/2)-1);
ky_idx_all = (-floor(cfg.Ny/2)):(ceil(cfg.Ny/2)-1);
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
kp = (2*pi/cfg.Tp)^2 / cfg.g;
sig_k = cfg.sig_k_factor * kp;
Sk = exp(-0.5*((kmag_all - kp)/sig_k).^2);
St1 = exp(-0.5*(angle(exp(1i*(theta_all - cfg.theta1)))/cfg.sig_t1).^2);
St2 = exp(-0.5*(angle(exp(1i*(theta_all - cfg.theta2)))/cfg.sig_t2).^2);
W = Sk .* (cfg.w1*St1 + cfg.w2*St2);
valid = W > 1e-8 * max(W);
kx_all = kx_all(valid);
ky_all = ky_all(valid);
W = W(valid);
if numel(W) < cfg.component_count
    error('Requested %d retained components, but only %d candidate components are available.', cfg.component_count, numel(W));
end
[~, idx_sort] = sort(W, 'descend');
idx = idx_sort(1:cfg.component_count);
kx = kx_all(idx);
ky = ky_all(idx);
Wsel = W(idx);
amp_raw = Wsel / max(Wsel);
m0_raw = 0.5 * sum(amp_raw.^2);
Hs_raw = 4 * sqrt(max(m0_raw, eps));
amp = amp_raw * (cfg.Hs_target / Hs_raw);
xf = 0.5 * cfg.Lx;
yf = 0.5 * cfg.Ly;
omega_lin = sqrt(cfg.g*hypot(kx,ky).*tanh(cfg.h*hypot(kx,ky)));
phase = -(kx*xf + ky*yf) + omega_lin*cfg.t;
cfg.a = amp .* cos(phase);
cfg.b = amp .* sin(phase);
cfg.kx = kx;
cfg.ky = ky;
end

function export_case(casesDir, cfg)
caseDir = fullfile(casesDir, cfg.case_id);
inputsDir = fullfile(caseDir, 'inputs');
refDir = fullfile(caseDir, 'reference', 'matlab');
mkdir_if_needed(inputsDir);
mkdir_if_needed(refDir);
write_column_csv(fullfile(inputsDir, 'a.csv'), cfg.a);
write_column_csv(fullfile(inputsDir, 'b.csv'), cfg.b);
write_column_csv(fullfile(inputsDir, 'kx.csv'), cfg.kx);
write_column_csv(fullfile(inputsDir, 'ky.csv'), cfg.ky);

export_reference = ~isfield(cfg, 'export_reference') || cfg.export_reference;
if export_reference
    coeff_repeats = get_benchmark_repeats(cfg);
    coeff_times = zeros(coeff_repeats, 1);
    recon_times = zeros(coeff_repeats, 1);
    kin_times = zeros(coeff_repeats, 1);
    for r = 1:coeff_repeats
        tic;
        coeffs = mf12_spectral_coefficients(cfg.order, cfg.g, cfg.h, cfg.a, cfg.b, cfg.kx, cfg.ky, cfg.Ux, cfg.Uy, struct('enable_subharmonic', false));
        coeff_times(r) = toc;
        coeffs.third_order_subharmonic_mode = cfg.subharmonic_mode;
        tic;
        [eta, phi, X, Y] = mf12_spectral_surface(coeffs, cfg.Lx, cfg.Ly, cfg.Nx, cfg.Ny, cfg.t);
        recon_times(r) = toc;
        tic;
        [u, v, w, p, phi_vol, uV, vV, a_x, a_y] = mf12_spectral_kinematics(coeffs, cfg.Lx, cfg.Ly, cfg.Nx, cfg.Ny, cfg.z_kinematics, cfg.t);
        kin_times(r) = toc;
    end
    writematrix(eta, fullfile(refDir, 'eta.csv'));
    writematrix(phi, fullfile(refDir, 'phi.csv'));
    write_column_csv(fullfile(refDir, 'x.csv'), X(1,:));
    write_column_csv(fullfile(refDir, 'y.csv'), Y(:,1));
    writematrix(u, fullfile(refDir, 'u.csv'));
    writematrix(v, fullfile(refDir, 'v.csv'));
    writematrix(w, fullfile(refDir, 'w.csv'));
    writematrix(p, fullfile(refDir, 'p.csv'));
    writematrix(phi_vol, fullfile(refDir, 'phi_vol.csv'));
    writematrix(uV, fullfile(refDir, 'uV.csv'));
    writematrix(vV, fullfile(refDir, 'vV.csv'));
    writematrix(a_x, fullfile(refDir, 'a_x.csv'));
    writematrix(a_y, fullfile(refDir, 'a_y.csv'));

    result = struct();
    result.case_id = cfg.case_id;
    result.runtime = struct('repeats', coeff_repeats, 'mean_coefficient_s', mean(coeff_times), 'mean_reconstruction_s', mean(recon_times), 'mean_kinematics_s', mean(kin_times), 'mean_total_s', mean(coeff_times + recon_times + kin_times), 'best_total_s', min(coeff_times + recon_times + kin_times));
    result.metadata = struct('language', 'matlab', 'matlab_version', version, 'fft_backend', 'MATLAB ifft2', 'implementation', 'mf12_spectral_coefficients + mf12_spectral_surface + mf12_spectral_kinematics');
    write_json(fullfile(refDir, 'result.json'), result);
end

manifest = struct();
manifest.case_id = cfg.case_id;
manifest.description = cfg.description;
manifest.purpose = cfg.purpose;
manifest.inputs = struct('order', cfg.order, 'g', cfg.g, 'h', cfg.h, 'Ux', cfg.Ux, 'Uy', cfg.Uy, 'Lx', cfg.Lx, 'Ly', cfg.Ly, 'Nx', cfg.Nx, 'Ny', cfg.Ny, 't', cfg.t, 'z_kinematics', cfg.z_kinematics, 'subharmonic_mode', cfg.subharmonic_mode);
manifest.arrays = struct('a', fullfile('inputs', 'a.csv'), 'b', fullfile('inputs', 'b.csv'), 'kx', fullfile('inputs', 'kx.csv'), 'ky', fullfile('inputs', 'ky.csv'));
if export_reference
    manifest.reference = struct('arrays', struct( ...
        'eta', fullfile('reference', 'matlab', 'eta.csv'), ...
        'phi', fullfile('reference', 'matlab', 'phi.csv'), ...
        'x', fullfile('reference', 'matlab', 'x.csv'), ...
        'y', fullfile('reference', 'matlab', 'y.csv'), ...
        'u', fullfile('reference', 'matlab', 'u.csv'), ...
        'v', fullfile('reference', 'matlab', 'v.csv'), ...
        'w', fullfile('reference', 'matlab', 'w.csv'), ...
        'p', fullfile('reference', 'matlab', 'p.csv'), ...
        'phi_vol', fullfile('reference', 'matlab', 'phi_vol.csv'), ...
        'uV', fullfile('reference', 'matlab', 'uV.csv'), ...
        'vV', fullfile('reference', 'matlab', 'vV.csv'), ...
        'a_x', fullfile('reference', 'matlab', 'a_x.csv'), ...
        'a_y', fullfile('reference', 'matlab', 'a_y.csv')), ...
        'metadata', fullfile('reference', 'matlab', 'result.json'));
else
    manifest.reference = struct();
end
manifest.tolerances = cfg.tolerances;
manifest.tolerances.u_max_abs_err = 1e-6;
manifest.tolerances.v_max_abs_err = 1e-6;
manifest.tolerances.w_max_abs_err = 1e-6;
manifest.tolerances.p_max_abs_err = 1e-6;
manifest.tolerances.phi_vol_max_abs_err = 1e-6;
manifest.tolerances.uV_max_abs_err = 1e-6;
manifest.tolerances.vV_max_abs_err = 1e-6;
manifest.tolerances.a_x_max_abs_err = 1e-6;
manifest.tolerances.a_y_max_abs_err = 1e-6;
write_json(fullfile(caseDir, 'case.json'), manifest);
fprintf('Exported case: %s\n', caseDir);
end

function write_column_csv(pathStr, vals)
writematrix(vals(:), pathStr);
end

function write_json(pathStr, s)
fid = fopen(pathStr, 'w');
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', jsonencode(s, 'PrettyPrint', true));
end

function mkdir_if_needed(pathStr)
if ~exist(pathStr, 'dir'), mkdir(pathStr); end
end

function maybe_export(caseId, exportFn, caseFilter)
if isempty(caseFilter) || strcmp(caseFilter, caseId)
    exportFn();
end
end

function coeff_repeats = get_benchmark_repeats(cfg)
coeff_repeats = strcmp(cfg.purpose, 'benchmark') * 2 + 1;
overrideStr = strtrim(getenv('MF12_EXPORT_BENCHMARK_REPEATS'));
if isempty(overrideStr)
    return;
end
overrideVal = str2double(overrideStr);
if isfinite(overrideVal) && overrideVal >= 1
    coeff_repeats = max(1, round(overrideVal));
end
end
