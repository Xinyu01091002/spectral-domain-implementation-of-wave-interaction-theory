clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
run(fullfile(rootDir, 'setup_paths.m'));

fprintf('=== Minimal Manuscript Reproduction ===\n');
fprintf('Repository root: %s\n', rootDir);
fprintf('This script runs one representative manuscript-facing verification case.\n\n');

fprintf('Step 1/1: directional wave-group third-order potential comparison\n');
run(fullfile(rootDir, 'verification', 'plot_phi3_wavegroup_lines.m'));

fprintf('\nDone.\n');
fprintf('Generated outputs are written under: %s\n', fullfile(rootDir, 'outputs'));
