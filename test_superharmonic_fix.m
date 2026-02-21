
try
    g=9.81; h=10; Nx=32; Ny=32; Lx=100; Ly=100;
    a = [0.1, 0.1]; b = [0, 0]; kx = [0.1, 0.15]; ky = [0, 0]; Ux=0; Uy=0;
    
    fprintf('Testing coeffsMF12_superharmonic (Order 2)...\n');
    coeffs = coeffsMF12_superharmonic(2, g, h, a, b, kx, ky, Ux, Uy);
    
    if isfield(coeffs, 'F_npm')
        fprintf('F_npm exists. Size: %d\n', length(coeffs.F_npm));
        % For N=2, pairs=1. Now with pm=[1,-1], we expect 2 entries?
        % My modified code: pairs = N*(N-1). For N=2, 2*1 = 2.
        % Correct.
        if length(coeffs.F_npm) == 2
           fprintf('SUCCESS: Found 2 interaction terms (Sum + Sub).\n');
        else
           fprintf('FAILURE: Expected 2 terms, got %d.\n', length(coeffs.F_npm));
        end
    else
        fprintf('FAILURE: F_npm missing.\n');
    end
    
    fprintf('Testing reconstruction...\n');
    [eta, ~] = surfaceMF12_spectral(coeffs, Lx, Ly, Nx, Ny, 0);
    fprintf('Reconstruction successful. Max eta: %f\n', max(eta(:)));
catch e
    fprintf('Error: %s\n', e.message);
    disp(e.stack(1));
end
