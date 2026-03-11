function setup_paths()
%SETUP_PATHS Add the bundled MF12 source tree to the MATLAB path.
%
% Run this once from the repository root before executing examples, tests,
% or manuscript reproduction scripts.

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'irregularWavesMF12', 'Source'));
rehash toolboxcache;

fprintf('MF12 source path added from: %s\n', fullfile(rootDir, 'irregularWavesMF12', 'Source'));
end
