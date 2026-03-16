clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
matlabDir = fileparts(scriptDir);
run(fullfile(matlabDir, 'setup_paths.m'));

fprintf('=== Minimal Manuscript Reproduction ===\n');
fprintf('MATLAB root: %s\n', matlabDir);
fprintf('This script runs one representative manuscript-facing verification case.\n\n');

fprintf('Step 1/1: directional wave-group third-order potential comparison\n');
run(fullfile(matlabDir, 'verification', 'plot_phi3_wavegroup_lines.m'));

fprintf('\nDone.\n');
fprintf('Generated outputs are written under: %s\n', fullfile(fileparts(matlabDir), 'outputs'));
