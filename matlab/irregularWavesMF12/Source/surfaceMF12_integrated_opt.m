function [eta, phi, X, Y] = surfaceMF12_integrated_opt(order, g, h, a, b, kx, ky, Ux, Uy, Lx, Ly, Nx, Ny, t, opts)
% SURFACEMF12_INTEGRATED_OPT
% Wrapper function for MF12 calculation combining coefficients and spectral reconstruction.
% 
% Single-workflow wrapper: mf12_spectral_coefficients + mf12_spectral_surface.
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
    if ~isfield(opts, 'legacyMode'),   opts.legacyMode = false; end

    if opts.legacyMode
        fprintf('Integrated MF12 (Order %d) [Legacy Wrapper: Superharmonic Coeffs]\n', order);
        coeffs = mf12_spectral_coefficients(order, g, h, a, b, kx, ky, Ux, Uy);
        [eta, phi, X, Y] = mf12_spectral_surface(coeffs, Lx, Ly, Nx, Ny, t);
        return;
    end

    fprintf('Integrated MF12 (Order %d) [Wrapper Mode: Superharmonic Coeffs]\n', order);
    coeffs = mf12_spectral_coefficients(order, g, h, a, b, kx, ky, Ux, Uy);
    [eta, phi, X, Y] = mf12_spectral_surface(coeffs, Lx, Ly, Nx, Ny, t);
    
end
