function [eta, phiS, X, Y] = surfaceMF12_spectral(coeffs, Lx, Ly, Nx, Ny, t)
% Spectral domain implementation of the surface elevation and potential reconstruction.
% This function is significantly faster than surfaceMF12 for large spatial fields 
% as it utilizes FFT (O(N log N)) instead of direct summation (O(N_waves * Nx * Ny)).
%
% Usage:
%   [eta, phiS, X, Y] = surfaceMF12_spectral(coeffs, Lx, Ly, Nx, Ny, t)
%
% Inputs:
%   coeffs - Structure returned by coeffsMF12 or coeffsMF12_superharmonic
%   Lx, Ly - Domain size in x and y directions (meters)
%   Nx, Ny - Number of grid points in x and y directions (should be power of 2 for speed)
%   t      - Time (scalar)
%
% Outputs:
%   eta    - Free surface elevation field (Ny x Nx)
%   phiS   - Velocity potential at free surface (Ny x Nx)
%   X, Y   - Meshgrid arrays of coordinates
%
% The function assumes the domain is periodic or the spectrum is band-limited within the
% Nyquist frequency defined by the grid resolution.

    % Define grid and wavenumbers
    dx = Lx / Nx;
    dy = Ly / Ny;
    x_axis = (0:Nx-1) * dx;
    y_axis = (0:Ny-1) * dy;
    [X, Y] = meshgrid(x_axis, y_axis);

    progress_enabled = strcmp(getenv('MF12_PROGRESS'), '1');
    t_total = tic;
    t_accum = tic;
    
    dkx = 2*pi/Lx;
    dky = 2*pi/Ly;

    % Initialize Spectral Grids
    % Note: MATLAB FFT origin is at (1,1). 
    % We will accumulate complex coefficients here.
    spec_eta = zeros(Ny, Nx); 
    spec_phi = zeros(Ny, Nx);

    % Accumulators (delay accumarray to the end for speed)
    eta_idx_x = [];
    eta_idx_y = [];
    eta_vals = [];
    phi_idx_x = [];
    phi_idx_y = [];
    phi_vals = [];

    % Helper variables
    Ux = coeffs.Ux;
    Uy = coeffs.Uy;
    superOnly = isfield(coeffs, 'superharmonic_only') && coeffs.superharmonic_only;
    
    % --- 1. First Order Terms ---
    % eta = sum( a*cos(theta) + b*sin(theta) )
    % phi = sum( (mu+muStar) * (a*sin(theta) - b*cos(theta)) )
    % theta = omega*t - k*x
    % Complex base: Z_linear = (a + 1i*b) * exp(-1i * omega * t)
    %
    % Note on Phi: 
    % a*sin - b*cos = Real( (a + ib) * (-i) * exp(...) ) ?
    % Let's check: (a+ib)(-i)(cos+isin) = (-ai + b)(cos+isin) = -ai cos + a sin + b cos + ib sin -> Real part: a sin + b cos. 
    % Wait, we want (a sin - b cos).
    % Try Z_phi_base = 1i * (a + 1i*b) * exp(-1i*omega*t)
    % = (ai - b) * (cos - i sin) [conjugate basis?]
    % Let's stick to the derivation:
    % We map +k to index corresponding to exp(+ikx) in IFFT.
    % Coeff for IFFT (+k): Z = (A + 1i*B) * exp(-1i*omega*t).
    % Re(Z * exp(ikx)) = Re( (A+iB)(cosO - i sinO)(cosK + i sinK) )... this is getting messy.
    %
    % Simpler: 
    % Term: A cos(wt - kx) + B sin(wt - kx).
    % Matches Re { (A + iB) * exp(-i(wt - kx)) } = Re { (A+iB)e^{-iwt} * e^{ikx} }.
    % This is exactly a coefficient C = (A+iB)e^{-iwt} for the basis function e^{ikx}.
    % MATLAB ifft(Y) computes sum( Y(k) * exp(i 2pi k n / N) ). This is the e^{ikx} basis.
    % So, Eta Coeff = (a + 1i*b) .* exp(-1i * omega * t).
    
    % For Phi: (mu) * ( A sin(theta) - B cos(theta) ).
    % We want -B cos + A sin.
    % Note that A cos + B sin corresponds to complex amplitude (A + iB).
    % -B cos + A sin is basically rotating (A + iB) by 90 degrees?
    % (A + iB) * (i) = iA - B = -B + iA.
    % Re { (-B + iA) e^{-iwt} e^{ikx} } = -B cos + A sin. Correct.
    % So Phi Coeff = (mu) * (a + 1i*b) * 1i .* exp(-1i * omega * t).
    
    % Linear Arrays
    kx = coeffs.kx(:);
    ky = coeffs.ky(:);
    omega = coeffs.omega(:);
    
    % Construct Linear Coefficients
    % Fix: Use (a - ib) to correspond to physical a*cos(theta) + b*sin(theta)
    % when reconstructed via IFFT (which uses exp(+ikx)).
    Z_lin = (coeffs.a(:) + 1i*coeffs.b(:)) .* exp(-1i * omega * t);
    
    % Accumulate Linear Eta (G=1 implicitly)
    accumulate_spectrum(kx, ky, Z_lin); 
    
    % Accumulate Linear Phi
    % Target: mu * (a*sin - b*cos)
    % Z_phi = Z_lin * (1i).
    % For Z_lin = (a+ib)e^{-iwt}, Real(i*Z_lin*e^{ikx}) = a*sin(theta) - b*cos(theta).
    mu_total = coeffs.mu(:) + coeffs.muStar(:);
    accumulate_spectrum(kx, ky, Z_lin .* mu_total * (1i), 'phi');


    % --- 2. Second Order ---
    if isfield(coeffs, 'G_2') % Self-Self
        % Terms are 2*theta. Effective k = 2*k, omega = 2*omega.
        % Eta: G_2 * (A_2 cos + B_2 sin)
        % Phi: mu_2 * (A_2 sin - B_2 cos)
        kx2 = 2*coeffs.kx(:);
        ky2 = 2*coeffs.ky(:);
        om2 = 2*coeffs.omega(:); % Note: usually it's phase doubling.
        
        Z_2 = (coeffs.A_2(:) + 1i*coeffs.B_2(:)) .* exp(-1i * om2 * t);
        
        accumulate_spectrum(kx2, ky2, Z_2 .* coeffs.G_2(:), 'eta');
        accumulate_spectrum(kx2, ky2, Z_2 .* coeffs.mu_2(:) * (1i), 'phi');
    end

    if isfield(coeffs, 'G_npm') % Sum Interactions
        % k = k_npm, omega = omega_npm
        Z_npm = (coeffs.A_npm(:) + 1i*coeffs.B_npm(:)) .* exp(-1i * coeffs.omega_npm(:) * t);
        
        accumulate_spectrum(coeffs.kx_npm(:), coeffs.ky_npm(:), Z_npm .* coeffs.G_npm(:), 'eta');
        accumulate_spectrum(coeffs.kx_npm(:), coeffs.ky_npm(:), Z_npm .* coeffs.mu_npm(:) * (1i), 'phi');
    end

    % --- 3. Third Order ---
    if isfield(coeffs, 'G_3') % Self-Self-Self (3rd Harmonic)
        kx3 = 3*coeffs.kx(:);
        ky3 = 3*coeffs.ky(:);
        om3 = 3*coeffs.omega(:); % Approximation of phase 3*theta
        
        Z_3 = (coeffs.A_3(:) + 1i*coeffs.B_3(:)) .* exp(-1i * om3 * t);
        accumulate_spectrum(kx3, ky3, Z_3 .* coeffs.G_3(:), 'eta');
        accumulate_spectrum(kx3, ky3, Z_3 .* coeffs.mu_3(:) * (1i), 'phi');
    end
    
    % --- Helper for correct B sign ---
    % eta = A*cos(theta) + B*sin(theta) = Real{ (A - iB) * exp(i*theta) }
    % theta = w*t - k*x.   exp(i*theta) = exp(iwt) * exp(-ikx).
    % FFT basis is exp(ikx). So we map to input k with coefficient Z?
    % If we map to k: term is Z * exp(ikx).
    % We want Real{ (A - iB) * exp(iwt) * exp(-ikx) }.
    % This is Real{ (A - iB) * exp(iwt) * conj(exp(ikx)) }.
    % This implies specific symmetry or just placing at -k?
    % MATLAB ifft(X) computes sum(X(k) exp(ikx)).
    % We want A*cos(wt - kx).
    % cos(wt - kx) = 0.5 * ( exp(i(wt-kx)) + exp(-i(wt-kx)) )
    %              = 0.5 * e^{iwt} e^{-ikx} + 0.5 * e^{-iwt} e^{ikx}.
    % The coefficient for e^{ikx} (spatial basis) is 0.5 * e^{-iwt}.
    % And we have amplitude A.
    % So coeff at frequency +k should be 0.5 * (A+iB_term) * e^{-iwt}.
    %
    % BUT, our spectrum is one-sided or we just fill k indices?
    % If we fill index corresponding to +k, we are setting coefficient for e^{ikx}.
    % The term we want is e^{-ikx}.
    % So we should fill index corresponding to -k? 
    % Or rely on Real() at the end? 
    % If we take Real(ifft), we are essentially summing (c_k e^{ikx} + c_{-k} e^{-ikx}).
    % If we only fill +k in spectrum with Z, and -k is 0 (or not set), ifft will limit to complex.
    % Real() will take 0.5 * (Z e^{ikx} + conj(Z) e^{-ikx}).
    % We want A cos(wt - kx).
    % If we set Z for +k? 
    % Z e^{ikx} -> Real -> 0.5( Z e^{ikx} + conj(Z) e^{-ikx} ).
    % We want cos(wt - kx).
    % cos(wt - kx) = cos(-wt + kx). (Even function)
    % = Real { e^{i(-wt + kx)} } = Real { e^{-iwt} e^{ikx} }.
    % So coefficient for e^{ikx} is e^{-iwt}.
    % So Z_k = e^{-iwt}. 
    % And we multiply by A. So Z = A * exp(-1i * omega * t).
    % This matches my code. 
    % Now what about B?
    % B * sin(wt - kx) = B * sin(-(kx - wt)) = -B * sin(kx - wt).
    % = -B * Real { -i * e^{i(kx - wt)} }   (since sin(x) = Re{-i e^ix})
    % = Real { i * B * e^{i(kx - wt)} }
    % = Real { i * B * e^{-iwt} * e^{ikx} }.
    % So coefficient is  i * B * e^{-iwt}.
    % Total Z = (A + iB) * exp(-1i * omega * t).
    % So (A + 1i*B) IS correct for the +k spatial mode.
    
    % Wait? A*cos + B*sin.
    % Term is A cos(wt - kx) + B sin(wt - kx).
    % = A cos(kx - wt) - B sin(kx - wt).
    % = Re { A e^{i(kx-wt)} } - Re { -i B e^{i(kx-wt)} }
    % = Re { (A + iB) e^{i(kx-wt)} }
    % = Re { (A + iB) e^{-iwt} e^{ikx} }.
    % Coefficient for e^{ikx} is (A + iB) e^{-iwt}.
    % So my Z calculation (A + 1i*B) was CORRECT.
    
    % Okay, so the sign is not the issue. The masking is.

    if isfield(coeffs, 'G_np2m') % Double Sum (n+2m)
        len_coeffs = length(coeffs.G_np2m);

        if superOnly
            mask = true(len_coeffs, 1);
        else
            % Full: coeffsMF12 stores [Sum, Diff, Sum, Diff...] for pm=[1, -1].
            % Keep only pm=1 (superharmonics), matching surfaceMF12_new behavior.
            mask = 1:2:len_coeffs;
        end

        Z_np2m = (coeffs.A_np2m(mask) + 1i*coeffs.B_np2m(mask)) .* exp(-1i * coeffs.omega_np2m(mask) * t);
        
        % Ensure column vectors for Z and coefficients
        Z_np2m = Z_np2m(:);
        G_val = coeffs.G_np2m(mask); G_val = G_val(:);
        mu_val = coeffs.mu_np2m(mask); mu_val = mu_val(:);
        k_x = coeffs.kx_np2m(mask); k_x = k_x(:);
        k_y = coeffs.ky_np2m(mask); k_y = k_y(:);
        
        accumulate_spectrum(k_x, k_y, Z_np2m .* G_val, 'eta');
        accumulate_spectrum(k_x, k_y, Z_np2m .* mu_val * (1i), 'phi');
    end

    if isfield(coeffs, 'G_2npm') % Double Sum (2n+m)
        len_coeffs = length(coeffs.G_2npm);

        if superOnly
            mask = true(len_coeffs, 1);
        else
            mask = 1:2:len_coeffs;
        end
        
        Z_2npm = (coeffs.A_2npm(mask) + 1i*coeffs.B_2npm(mask)) .* exp(-1i * coeffs.omega_2npm(mask) * t);
        
        Z_2npm = Z_2npm(:);
        G_val = coeffs.G_2npm(mask); G_val = G_val(:);
        mu_val = coeffs.mu_2npm(mask); mu_val = mu_val(:);
        k_x = coeffs.kx_2npm(mask); k_x = k_x(:);
        k_y = coeffs.ky_2npm(mask); k_y = k_y(:);
        
        accumulate_spectrum(k_x, k_y, Z_2npm .* G_val, 'eta');
        accumulate_spectrum(k_x, k_y, Z_2npm .* mu_val * (1i), 'phi');
    end
    
    if isfield(coeffs, 'G_npmpp') % Triple Sum (n+m+p)
        N_c = coeffs.N;
        % Triple sums
        % Count: N*(N-1)*(N-2)/6 -> unique triplets {n,m,p}.
        % Full: pmm=[1 -1], pmp=[1 -1]. 4 combinations per triplet.
        % Super only: 1 combination.
        
        len_coeffs = length(coeffs.G_npmpp);
        num_triplets = N_c*(N_c-1)*(N_c-2)/6;
        
        if superOnly
            mask = true(len_coeffs, 1);
        else
            % Reconstruct indices for pmm=1 && pmp=1
             idx = 0;
             mask_log = false(len_coeffs, 1);
            
             for n = 1:N_c
                for m = n+1:N_c
                    for pmm = [1 -1]
                        for p = m+1:N_c
                            for pmp = [1 -1]
                                idx = idx + 1;
                                if pmm == 1 && pmp == 1
                                    mask_log(idx) = true;
                                end
                            end
                        end
                    end
                end
             end
             mask = mask_log;
        end
        
        Z_npmpp = (coeffs.A_npmpp(mask) + 1i*coeffs.B_npmpp(mask)) .* exp(-1i * coeffs.omega_npmpp(mask) * t);
        
        % Factor of 2 required for triple sums (as per surfaceMF12_new)
        Z_npmpp = Z_npmpp * 2; 
 
        
        Z_npmpp = Z_npmpp(:);
        G_val = coeffs.G_npmpp(mask); G_val = G_val(:);
        mu_val = coeffs.mu_npmpp(mask); mu_val = mu_val(:);
        k_x = coeffs.kx_npmpp(mask); k_x = k_x(:);
        k_y = coeffs.ky_npmpp(mask); k_y = k_y(:);
        
        accumulate_spectrum(k_x, k_y, Z_npmpp .* G_val, 'eta');
        accumulate_spectrum(k_x, k_y, Z_npmpp .* mu_val * (1i), 'phi');
    end

    % --- Final Spectral Accumulation ---
    if progress_enabled
        fprintf('surfaceMF12_spectral: accumulation done in %.2fs\n', toc(t_accum));
    end
    t_grid = tic;
    if ~isempty(eta_vals)
        spec_eta = accumarray([eta_idx_y, eta_idx_x], eta_vals, [Ny, Nx]);
    end
    if ~isempty(phi_vals)
        spec_phi = accumarray([phi_idx_y, phi_idx_x], phi_vals, [Ny, Nx]);
    end
    if progress_enabled
        fprintf('surfaceMF12_spectral: accumarray done in %.2fs\n', toc(t_grid));
    end

    % --- Reconstruction ---
    % Perform Inverse FFT
    % Multiply by (Nx*Ny) because MATLAB's ifft includes a 1/N scaling, but our coefficients are already magnitudes.
    t_ifft = tic;
    eta = real(ifft2(spec_eta)) * (Nx * Ny);
    
    phi_wave = real(ifft2(spec_phi)) * (Nx * Ny);
    phiS = phi_wave + Ux*X + Uy*Y;
    if progress_enabled
        fprintf('surfaceMF12_spectral: ifft done in %.2fs\n', toc(t_ifft));
        fprintf('surfaceMF12_spectral: total time %.2fs\n', toc(t_total));
    end
    
    % --- Nested Helper Function for Binning ---
    function accumulate_spectrum(k_in_x, k_in_y, values, type)
        if nargin < 4, type = 'eta'; end
        
        % Map wavenumbers to grid indices
        % Index 1: k=0
        % Index 2..N/2: Positive k
        % Index N/2+1..N: Negative k

        ux = (k_in_x(:) / dkx);
        uy = (k_in_y(:) / dky);
        vals = values(:);

        valid = isfinite(ux) & isfinite(uy) & isfinite(vals);
        ux = ux(valid);
        uy = uy(valid);
        vals = vals(valid);
        if isempty(vals), return; end

        % Bilinear deposition in spectral-index space.
        % For exact-grid k, this collapses to single-bin deposition.
        ix0 = floor(ux);
        iy0 = floor(uy);
        fx = ux - ix0;
        fy = uy - iy0;

        % Snap near-integer values to avoid tiny floating-point leakage.
        tol = 1e-12;
        fx(abs(fx) < tol) = 0;
        fy(abs(fy) < tol) = 0;
        fx(abs(fx-1) < tol) = 1;
        fy(abs(fy-1) < tol) = 1;

        ix1 = ix0 + 1;
        iy1 = iy0 + 1;

        idx_x00 = mod(ix0, Nx) + 1;
        idx_y00 = mod(iy0, Ny) + 1;
        idx_x10 = mod(ix1, Nx) + 1;
        idx_y10 = idx_y00;
        idx_x01 = idx_x00;
        idx_y01 = mod(iy1, Ny) + 1;
        idx_x11 = idx_x10;
        idx_y11 = idx_y01;

        w00 = (1-fx).*(1-fy);
        w10 = fx.*(1-fy);
        w01 = (1-fx).*fy;
        w11 = fx.*fy;

        idx_x = [idx_x00; idx_x10; idx_x01; idx_x11];
        idx_y = [idx_y00; idx_y10; idx_y01; idx_y11];
        vals4 = [vals.*w00; vals.*w10; vals.*w01; vals.*w11];

        nz = (abs(vals4) > 0);
        idx_x = idx_x(nz);
        idx_y = idx_y(nz);
        vals4 = vals4(nz);

        if strcmp(type, 'eta')
            eta_idx_x = [eta_idx_x; idx_x];
            eta_idx_y = [eta_idx_y; idx_y];
            eta_vals = [eta_vals; vals4];
        else
            phi_idx_x = [phi_idx_x; idx_x];
            phi_idx_y = [phi_idx_y; idx_y];
            phi_vals = [phi_vals; vals4];
        end
    end

end
