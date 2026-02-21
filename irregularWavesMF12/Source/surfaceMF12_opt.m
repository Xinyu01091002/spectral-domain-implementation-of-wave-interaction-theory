function [eta, phi, X, Y] = surfaceMF12_opt(order, g, h, a, b, kx, ky, Ux, Uy, Lx, Ly, Nx, Ny, t, opts)
%SURFACEMF12_OPT Compatibility wrapper for optimized MF12 reconstruction.
%   This entry point delegates to surfaceMF12_integrated_opt.

    if nargin < 15
        [eta, phi, X, Y] = surfaceMF12_integrated_opt(order, g, h, a, b, kx, ky, Ux, Uy, Lx, Ly, Nx, Ny, t);
    else
        [eta, phi, X, Y] = surfaceMF12_integrated_opt(order, g, h, a, b, kx, ky, Ux, Uy, Lx, Ly, Nx, Ny, t, opts);
    end
end
