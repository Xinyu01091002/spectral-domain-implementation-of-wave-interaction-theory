% Generate test data and check speed for coeffsMF12.m

% Add necessary paths
if exist('irregularWavesMF12', 'dir')
    addpath(genpath('irregularWavesMF12'));
else
    warning('Folder irregularWavesMF12 not found. Make sure you are in the correct directory.');
end

% Parameters
g = 9.81;
h = 20; % Depth in meters
Ux = 0;
Uy = 0;
order = 3; % 3rd order is likely the most computationally intensive

% Test cases configuration
% N is the number of wave components to simulate.
% As N increases, the computational cost grows (likely O(N^2)).
test_cases = [
    struct('N', 20, 'iterations', 10)   % Small number of components, many iterations
    struct('N', 50, 'iterations', 10)    % Medium
    struct('N', 100, 'iterations', 10)    % Large
    struct('N', 200, 'iterations', 10)
];

fprintf('--------------------------------------------------\n');
fprintf('Running speed test for coeffsMF12 (Order = %d)\n', order);
fprintf('--------------------------------------------------\n');

for i = 1:length(test_cases)
    N = test_cases(i).N;
    iters = test_cases(i).iterations;
    
    % Generate random wave components
    % Wavenumbers roughly corresponding to typical ocean waves
    % Use deterministic random numbers for reproducibility if needed, 
    % but here just random is fine.
    kx = linspace(0.1, 1.5, N)'; 
    ky = (rand(N, 1) - 0.5) * 0.2; % Small directional spread
    
    % Random amplitudes (small enough to be physically realistic sum)
    a = (rand(N, 1) - 0.5) * 0.1;
    b = (rand(N, 1) - 0.5) * 0.1;
    
    fprintf('\nCase %d: N = %d wave components\n', i, N);
    
    % Warm-up run (to handle JIT compilation overhead)
    try
        coeffsMF12(order, g, h, a, b, kx, ky, Ux, Uy);
    catch ME
        fprintf('Error during warm-up: %s\n', ME.message);
        if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
             fprintf('Make sure coeffsMF12 is on the path.\n');
        end
        continue;
    end
    
    % Timing loop
    fprintf('Running %d iterations...\n', iters);
    tic;
    for k = 1:iters
        % Suppress output if necessary, though coeffsMF12 usually returns struct
        dummy = coeffsMF12(order, g, h, a, b, kx, ky, Ux, Uy);
    end
    total_time = toc;
    
    avg_time = total_time / iters;
    
    fprintf('  [Full] Total time: %.4f s\n', total_time);
    fprintf('  [Full] Avg time per call: %.6f s\n', avg_time);
    
    % Projections for Full
    time_1k = avg_time * 1000;
    time_10k = avg_time * 10000;
    
    fprintf('  [Full] Projection for 1,000 runs:  %.2f s (%.2f min)\n', time_1k, time_1k/60);
    fprintf('  [Full] Projection for 10,000 runs: %.2f s (%.2f min)\n', time_10k, time_10k/60);
    
    % ----------------------------------------------------
    % Superharmonic Only Test
    % ----------------------------------------------------
    % Warm-up run
    try
        coeffsMF12_superharmonic(order, g, h, a, b, kx, ky, Ux, Uy);
    catch ME
        fprintf('Error during superharmonic warm-up: %s\n', ME.message);
        continue;
    end
    % linear_focus_envelope
    tic;
    for k = 1:iters
        dummy = coeffsMF12_superharmonic(order, g, h, a, b, kx, ky, Ux, Uy);
    end
    total_time_super = toc;
    avg_time_super = total_time_super / iters;
    
    fprintf('  [Super Only] Total time: %.4f s\n', total_time_super);
    fprintf('  [Super Only] Avg time per call: %.6f s\n', avg_time_super);
    
    % Speedup
    if avg_time_super > 0
        speedup = avg_time / avg_time_super;
        fprintf('  >>> Speedup factor: %.2fx\n', speedup);
    end

end

fprintf('\n--------------------------------------------------\n');
fprintf('Done.\n');
