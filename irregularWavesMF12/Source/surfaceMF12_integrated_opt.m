function [eta, phi, X, Y] = surfaceMF12_integrated_opt(order, g, h, a, b, kx, ky, Ux, Uy, Lx, Ly, Nx, Ny, t, opts)
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

    if nargin < 15 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'useStreaming'), opts.useStreaming = false; end
    if ~isfield(opts, 'legacyMode'),   opts.legacyMode = false; end

    if opts.legacyMode
        fprintf('Integrated MF12 (Order %d) [Legacy Wrapper: Superharmonic Coeffs]\n', order);
        coeffs = coeffsMF12_superharmonic(order, g, h, a, b, kx, ky, Ux, Uy);
        [eta, phi, X, Y] = surfaceMF12_spectral(coeffs, Lx, Ly, Nx, Ny, t);
        return;
    end

    if opts.useStreaming
        fprintf('Integrated MF12 (Order %d) [Streaming Superharmonic Mode]\n', order);
        [eta, phi, X, Y] = surfaceMF12_streaming_super(order, g, h, a, b, kx, ky, Ux, Uy, Lx, Ly, Nx, Ny, t);
    else
        fprintf('Integrated MF12 (Order %d) [Wrapper Mode: Superharmonic Coeffs]\n', order);
        coeffs = coeffsMF12_superharmonic(order, g, h, a, b, kx, ky, Ux, Uy);
        [eta, phi, X, Y] = surfaceMF12_spectral(coeffs, Lx, Ly, Nx, Ny, t);
    end
    
end
