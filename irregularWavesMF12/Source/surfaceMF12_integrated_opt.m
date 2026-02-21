function [eta, phi, X, Y] = surfaceMF12_integrated_opt(order, g, h, a, b, kx, ky, Ux, Uy, Lx, Ly, Nx, Ny, t)
% SURFACEMF12_INTEGRATED_OPT
% Wrapper function for MF12 calculation combining coefficients and spectral reconstruction.
% 
% Note: Ideally, this function performs "streamed" or "buffered" calculation 
% to avoid storing the full O(N^3) coefficient matrix. 
% For this verification/comparison script on reasonable Component counts (N < 500),
% we utilize the standard verified coeffsMF12 + surfaceMF12_spectral pipeline 
% which is numerically equivalent but consumes more memory.
%
% Inputs:
%   order - Calculation order (1, 2, or 3)
%   g, h  - Gravity, Depth
%   a, b  - Linear Amplitude components (vectors)
%   kx, ky- Wavenumber components (vectors)
%   Ux, Uy- Current (scalars)
%   Lx, Ly- Domain Size
%   Nx, Ny- Grid Points
%   t     - Time

    % Ensure dependencies are on path
    p = fileparts(mfilename('fullpath'));
    addpath(p);

    % Display info
    fprintf('Integrated MF12 (Order %d) [Wrapper Mode]\n', order);
    
    % 1. Compute Coefficients
    % In a true integrated optimizer, this would encompass the loop over interactions
    % directly adding to the grid. Here we use the pre-calc structure.
    % fprintf('  > Computing Coefficients...\n');
    coeffs = coeffsMF12(order, g, h, a, b, kx, ky, Ux, Uy);
    
    % 2. Reconstruct Surface
    % fprintf('  > Reconstructing Surface...\n');
    [eta, phi, X, Y] = surfaceMF12_spectral(coeffs, Lx, Ly, Nx, Ny, t);
    
end
