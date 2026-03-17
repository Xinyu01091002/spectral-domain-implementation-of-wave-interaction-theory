function benchmark_case_from_manifest(caseDir, outputJsonPath, repeats)
if nargin < 3
    repeats = 3;
end

scriptDir = fileparts(mfilename('fullpath'));
matlabDir = fileparts(scriptDir);
run(fullfile(matlabDir, 'setup_paths.m'));

manifest = jsondecode(fileread(fullfile(caseDir, 'case.json')));
inp = manifest.inputs;

a = readmatrix(fullfile(caseDir, manifest.arrays.a));
b = readmatrix(fullfile(caseDir, manifest.arrays.b));
kx = readmatrix(fullfile(caseDir, manifest.arrays.kx));
ky = readmatrix(fullfile(caseDir, manifest.arrays.ky));

a = a(:).';
b = b(:).';
kx = kx(:).';
ky = ky(:).';

coeffTimes = zeros(repeats, 1);
reconTimes = zeros(repeats, 1);
for r = 1:repeats
    tic;
    coeffs = mf12_spectral_coefficients(inp.order, inp.g, inp.h, a, b, kx, ky, inp.Ux, inp.Uy, struct('enable_subharmonic', false));
    coeffTimes(r) = toc;
    coeffs.third_order_subharmonic_mode = inp.subharmonic_mode;
    tic;
    mf12_spectral_surface(coeffs, inp.Lx, inp.Ly, inp.Nx, inp.Ny, inp.t);
    reconTimes(r) = toc;
end

runtime = struct();
runtime.repeats = repeats;
runtime.warmup = false;
runtime.mean_coefficient_s = mean(coeffTimes);
runtime.mean_reconstruction_s = mean(reconTimes);
runtime.mean_total_s = mean(coeffTimes + reconTimes);
runtime.best_total_s = min(coeffTimes + reconTimes);

summary = struct();
summary.case_id = manifest.case_id;
summary.order = inp.order;
summary.component_count = numel(a);
summary.domain = struct('Lx', inp.Lx, 'Ly', inp.Ly);
summary.grid = struct('Nx', inp.Nx, 'Ny', inp.Ny);
summary.runtime = runtime;

outDir = fileparts(outputJsonPath);
if ~isempty(outDir) && ~exist(outDir, 'dir')
    mkdir(outDir);
end

fid = fopen(outputJsonPath, 'w');
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', jsonencode(summary, 'PrettyPrint', true));
