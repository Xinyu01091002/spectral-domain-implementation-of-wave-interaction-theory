function [u, v, w, p, phi, uV, vV, a_x, a_y, X, Y] = mf12_spectral_kinematics(coeffs, Lx, Ly, Nx, Ny, z, t)
%MF12_SPECTRAL_KINEMATICS Spectral-domain reconstruction of MF12 kinematics on a constant-z plane.
% This is the spectral analogue of kinematicsMF12 for a regular horizontal grid.
%
% Usage:
%   [u, v, w, p, phi, uV, vV, a_x, a_y, X, Y] = ...
%       mf12_spectral_kinematics(coeffs, Lx, Ly, Nx, Ny, z, t)
%
% Inputs:
%   coeffs - Structure returned by mf12_direct_coefficients or mf12_spectral_coefficients
%   Lx, Ly - Domain size in x and y directions (meters)
%   Nx, Ny - Number of grid points in x and y directions
%   z      - Scalar vertical coordinate for the horizontal evaluation plane
%   t      - Time (scalar)
%
% Outputs:
%   u, v, w - Velocity components on the z-plane
%   p       - Dynamic pressure p^+/rho
%   phi     - Velocity potential
%   uV, vV  - Drag-style derived fields, accumulated to match kinematicsMF12
%   a_x, a_y- Total horizontal accelerations, accumulated to match kinematicsMF12
%   X, Y    - Meshgrid arrays of coordinates
%
% Notes:
%   - This implementation assumes a constant-z horizontal plane. That is the
%     spectral counterpart of evaluating kinematics on a regular x-y grid.
%   - The primary harmonic fields (phi, u, v, w, p) are reconstructed
%     spectrally.
%   - The derived fields (uV, vV, a_x, a_y) follow the original
%     kinematicsMF12 accumulation logic so they can be compared directly.

    if ~isscalar(z) || ~isfinite(z)
        error('mf12_spectral_kinematics:InvalidZ', ...
            'Input z must be a finite scalar. This function reconstructs a constant-z plane.');
    end

    dx = Lx / Nx;
    dy = Ly / Ny;
    x_axis = (0:Nx-1) * dx;
    y_axis = (0:Ny-1) * dy;
    [X, Y] = meshgrid(x_axis, y_axis);

    dkx = 2*pi/Lx;
    dky = 2*pi/Ly;
    Z = z + coeffs.h;

    spec_phi = zeros(Ny, Nx);
    spec_u   = zeros(Ny, Nx);
    spec_v   = zeros(Ny, Nx);
    spec_w   = zeros(Ny, Nx);
    spec_p   = zeros(Ny, Nx);

    spec_ut  = zeros(Ny, Nx);
    spec_vt  = zeros(Ny, Nx);
    spec_ux  = zeros(Ny, Nx);
    spec_uy  = zeros(Ny, Nx);
    spec_uz  = zeros(Ny, Nx);
    spec_vx  = zeros(Ny, Nx);
    spec_vy  = zeros(Ny, Nx);
    spec_vz  = zeros(Ny, Nx);

    superOnly = isfield(coeffs, 'superharmonic_only') && coeffs.superharmonic_only;

    mode3 = 'auto';
    if isfield(coeffs, 'third_order_subharmonic_mode')
        vmode = lower(string(coeffs.third_order_subharmonic_mode));
        if vmode == "include" || vmode == "skip" || vmode == "auto"
            mode3 = char(vmode);
        end
    end

    if superOnly
        waveGroupLike = false;
        nearResonant = false;
        skip3Sub = true;
    else
        [waveGroupLike, nearResonant] = detect_wavegroup_or_resonance(coeffs);
        if strcmp(mode3, 'include')
            skip3Sub = false;
        elseif strcmp(mode3, 'skip')
            skip3Sub = true;
        else
            skip3Sub = waveGroupLike || nearResonant;
        end
    end

    mask_np2m = [];
    mask_2npm = [];
    mask_npmpp = [];
    if isfield(coeffs, 'G_np2m')
        len_np2m = length(coeffs.G_np2m);
        if superOnly
            mask_np2m = true(len_np2m, 1);
        elseif skip3Sub
            mask_np2m = false(len_np2m, 1);
            mask_np2m(1:2:len_np2m) = true;
        else
            mask_np2m = true(len_np2m, 1);
        end
    end
    if isfield(coeffs, 'G_2npm')
        len_2npm = length(coeffs.G_2npm);
        if superOnly
            mask_2npm = true(len_2npm, 1);
        elseif skip3Sub
            mask_2npm = false(len_2npm, 1);
            mask_2npm(1:2:len_2npm) = true;
        else
            mask_2npm = true(len_2npm, 1);
        end
    end
    if isfield(coeffs, 'G_npmpp')
        len_npmpp = length(coeffs.G_npmpp);
        if superOnly
            mask_npmpp = true(len_npmpp, 1);
        elseif skip3Sub
            N_c = coeffs.N;
            mask_npmpp = false(len_npmpp, 1);
            idxm = 0;
            for n = 1:N_c
                for m = n+1:N_c
                    for pmm = [1 -1]
                        for pidx = m+1:N_c
                            for pmp = [1 -1]
                                idxm = idxm + 1;
                                if pmm == 1 && pmp == 1
                                    mask_npmpp(idxm) = true;
                                end
                            end
                        end
                    end
                end
            end
        else
            mask_npmpp = true(len_npmpp, 1);
        end
    end

    deposit_kinematics(coeffs.kx(:), coeffs.ky(:), coeffs.omega(:), coeffs.kappa(:), ...
        coeffs.F(:), coeffs.a(:), coeffs.b(:));

    if isfield(coeffs, 'G_2')
        deposit_kinematics(2*coeffs.kx(:), 2*coeffs.ky(:), 2*coeffs.omega(:), coeffs.kappa_2(:), ...
            coeffs.F_2(:), coeffs.A_2(:), coeffs.B_2(:));
    end

    if isfield(coeffs, 'G_npm')
        deposit_kinematics(coeffs.kx_npm(:), coeffs.ky_npm(:), coeffs.omega_npm(:), coeffs.kappa_npm(:), ...
            coeffs.F_npm(:), coeffs.A_npm(:), coeffs.B_npm(:));
    end

    if isfield(coeffs, 'F13')
        deposit_kinematics(coeffs.kx(:), coeffs.ky(:), coeffs.omega(:), coeffs.kappa(:), ...
            coeffs.F13(:), coeffs.a(:), coeffs.b(:));
    end

    if isfield(coeffs, 'G_3')
        deposit_kinematics(3*coeffs.kx(:), 3*coeffs.ky(:), 3*coeffs.omega(:), coeffs.kappa_3(:), ...
            coeffs.F_3(:), coeffs.A_3(:), coeffs.B_3(:));
    end

    if isfield(coeffs, 'G_np2m')
        mask = mask_np2m;
        deposit_kinematics(coeffs.kx_np2m(mask), coeffs.ky_np2m(mask), coeffs.omega_np2m(mask), coeffs.kappa_np2m(mask), ...
            coeffs.F_np2m(mask), coeffs.A_np2m(mask), coeffs.B_np2m(mask));
    end

    if isfield(coeffs, 'G_2npm')
        mask = mask_2npm;
        deposit_kinematics(coeffs.kx_2npm(mask), coeffs.ky_2npm(mask), coeffs.omega_2npm(mask), coeffs.kappa_2npm(mask), ...
            coeffs.F_2npm(mask), coeffs.A_2npm(mask), coeffs.B_2npm(mask));
    end

    if isfield(coeffs, 'G_npmpp')
        mask = mask_npmpp;
        deposit_kinematics(coeffs.kx_npmpp(mask), coeffs.ky_npmpp(mask), coeffs.omega_npmpp(mask), coeffs.kappa_npmpp(mask), ...
            coeffs.F_npmpp(mask), coeffs.A_npmpp(mask), coeffs.B_npmpp(mask));
    end

    phi_wave = real(ifft2(spec_phi)) * (Nx * Ny);
    u_wave   = real(ifft2(spec_u))   * (Nx * Ny);
    v_wave   = real(ifft2(spec_v))   * (Nx * Ny);
    w        = real(ifft2(spec_w))   * (Nx * Ny);
    p        = real(ifft2(spec_p))   * (Nx * Ny);

    u_t      = real(ifft2(spec_ut))  * (Nx * Ny);
    v_t      = real(ifft2(spec_vt))  * (Nx * Ny);
    u_x      = real(ifft2(spec_ux))  * (Nx * Ny);
    u_y      = real(ifft2(spec_uy))  * (Nx * Ny);
    u_z      = real(ifft2(spec_uz))  * (Nx * Ny);
    v_x      = real(ifft2(spec_vx))  * (Nx * Ny);
    v_y      = real(ifft2(spec_vy))  * (Nx * Ny);
    v_z      = real(ifft2(spec_vz))  * (Nx * Ny);

    phi = phi_wave + coeffs.Ux*X + coeffs.Uy*Y;
    u = u_wave + coeffs.Ux;
    v = v_wave + coeffs.Uy;

    [uV, vV, a_x, a_y] = accumulate_derived_outputs();

    function deposit_kinematics(kx_in, ky_in, omega_in, kappa_in, F_in, A_in, B_in)
        kxv = kx_in(:);
        kyv = ky_in(:);
        omegav = omega_in(:);
        kappav = kappa_in(:);
        Fv = F_in(:);
        Av = A_in(:);
        Bv = B_in(:);

        valid = isfinite(kxv) & isfinite(kyv) & isfinite(omegav) & ...
                isfinite(kappav) & isfinite(Fv) & isfinite(Av) & isfinite(Bv);
        kxv = kxv(valid);
        kyv = kyv(valid);
        omegav = omegav(valid);
        kappav = kappav(valid);
        Fv = Fv(valid);
        Av = Av(valid);
        Bv = Bv(valid);
        if isempty(kxv)
            return;
        end

        base = (Av + 1i*Bv) .* exp(-1i * omegav * t);
        coshZ = cosh(kappav * Z);
        sinhZ = sinh(kappav * Z);
        phi_amp = Fv .* coshZ;
        w_amp = Fv .* sinhZ .* kappav;

        phi_coeff = phi_amp .* base * (1i);
        u_coeff = -kxv .* phi_amp .* base;
        v_coeff = -kyv .* phi_amp .* base;
        w_coeff = w_amp .* base * (1i);
        p_coeff = -omegav .* phi_amp .* base;

        ut_coeff = (kxv .* omegav) .* phi_coeff;
        vt_coeff = (kyv .* omegav) .* phi_coeff;
        ux_coeff = -(kxv.^2) .* phi_coeff;
        uy_coeff = -(kxv .* kyv) .* phi_coeff;
        uz_coeff = -(kxv .* kappav) .* (Fv .* sinhZ) .* base;
        vx_coeff = -(kxv .* kyv) .* phi_coeff;
        vy_coeff = -(kyv.^2) .* phi_coeff;
        vz_coeff = -(kyv .* kappav) .* (Fv .* sinhZ) .* base;

        spec_phi = add_to_spectrum(spec_phi, kxv, kyv, phi_coeff);
        spec_u   = add_to_spectrum(spec_u,   kxv, kyv, u_coeff);
        spec_v   = add_to_spectrum(spec_v,   kxv, kyv, v_coeff);
        spec_w   = add_to_spectrum(spec_w,   kxv, kyv, w_coeff);
        spec_p   = add_to_spectrum(spec_p,   kxv, kyv, p_coeff);

        spec_ut  = add_to_spectrum(spec_ut,  kxv, kyv, ut_coeff);
        spec_vt  = add_to_spectrum(spec_vt,  kxv, kyv, vt_coeff);
        spec_ux  = add_to_spectrum(spec_ux,  kxv, kyv, ux_coeff);
        spec_uy  = add_to_spectrum(spec_uy,  kxv, kyv, uy_coeff);
        spec_uz  = add_to_spectrum(spec_uz,  kxv, kyv, uz_coeff);
        spec_vx  = add_to_spectrum(spec_vx,  kxv, kyv, vx_coeff);
        spec_vy  = add_to_spectrum(spec_vy,  kxv, kyv, vy_coeff);
        spec_vz  = add_to_spectrum(spec_vz,  kxv, kyv, vz_coeff);
    end

    function spec = add_to_spectrum(spec, k_in_x, k_in_y, values)
        ux = (k_in_x(:) / dkx);
        uy = (k_in_y(:) / dky);
        vals = values(:);

        valid = isfinite(ux) & isfinite(uy) & isfinite(vals);
        ux = ux(valid);
        uy = uy(valid);
        vals = vals(valid);
        if isempty(vals)
            return;
        end

        ix0 = floor(ux);
        iy0 = floor(uy);
        fx = ux - ix0;
        fy = uy - iy0;

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

        nz = abs(vals4) > 0;
        if any(nz)
            spec = spec + accumarray([idx_y(nz), idx_x(nz)], vals4(nz), [Ny, Nx]);
        end
    end

    function [uV_loc, vV_loc, a_x_loc, a_y_loc] = accumulate_derived_outputs()
        u_loc = coeffs.Ux + zeros(Ny, Nx);
        v_loc = coeffs.Uy + zeros(Ny, Nx);
        w_loc = zeros(Ny, Nx);
        V0 = sqrt(coeffs.Ux^2 + coeffs.Uy^2) * coeffs.Ux;
        uV_loc = V0 + zeros(Ny, Nx);
        vV_loc = coeffs.Uy * sqrt(coeffs.Ux^2 + coeffs.Uy^2) + zeros(Ny, Nx);
        a_x_loc = zeros(Ny, Nx);
        a_y_loc = zeros(Ny, Nx);

        for n = 1:coeffs.N
            accumulate_one(coeffs.omega(n), coeffs.kx(n), coeffs.ky(n), coeffs.kappa(n), ...
                coeffs.F(n), coeffs.a(n), coeffs.b(n));
        end

        if isfield(coeffs, 'G_2')
            for n = 1:coeffs.N
                accumulate_one(2*coeffs.omega(n), coeffs.kx_2(n), coeffs.ky_2(n), coeffs.kappa_2(n), ...
                    coeffs.F_2(n), coeffs.A_2(n), coeffs.B_2(n));
            end

            cnm = 0;
            for n = 1:coeffs.N
                for m = n+1:coeffs.N
                    for pm = [1 -1]
                        cnm = cnm + 1;
                        accumulate_one(coeffs.omega_npm(cnm), coeffs.kx_npm(cnm), coeffs.ky_npm(cnm), coeffs.kappa_npm(cnm), ...
                            coeffs.F_npm(cnm), coeffs.A_npm(cnm), coeffs.B_npm(cnm));
                    end
                end
            end
        end

        if isfield(coeffs, 'F13')
            for n = 1:coeffs.N
                accumulate_one(coeffs.omega(n), coeffs.kx(n), coeffs.ky(n), coeffs.kappa(n), ...
                    coeffs.F13(n), coeffs.a(n), coeffs.b(n));

                accumulate_one(3*coeffs.omega(n), 3*coeffs.kx(n), 3*coeffs.ky(n), coeffs.kappa_3(n), ...
                    coeffs.F_3(n), coeffs.A_3(n), coeffs.B_3(n));
            end

            cnm = 0;
            pair_super = 0;
            for n = 1:coeffs.N
                for m = n+1:coeffs.N
                    for pm = [1 -1]
                        cnm = cnm + 1;
                        if superOnly
                            if pm == 1
                                pair_super = pair_super + 1;
                                accumulate_one(coeffs.omega_np2m(pair_super), coeffs.kx_np2m(pair_super), coeffs.ky_np2m(pair_super), coeffs.kappa_np2m(pair_super), ...
                                    coeffs.F_np2m(pair_super), coeffs.A_np2m(pair_super), coeffs.B_np2m(pair_super));
                                accumulate_one(coeffs.omega_2npm(pair_super), coeffs.kx_2npm(pair_super), coeffs.ky_2npm(pair_super), coeffs.kappa_2npm(pair_super), ...
                                    coeffs.F_2npm(pair_super), coeffs.A_2npm(pair_super), coeffs.B_2npm(pair_super));
                            end
                        else
                            if isempty(mask_np2m) || mask_np2m(cnm)
                                accumulate_one(coeffs.omega_np2m(cnm), coeffs.kx_np2m(cnm), coeffs.ky_np2m(cnm), coeffs.kappa_np2m(cnm), ...
                                    coeffs.F_np2m(cnm), coeffs.A_np2m(cnm), coeffs.B_np2m(cnm));
                            end
                            if isempty(mask_2npm) || mask_2npm(cnm)
                                accumulate_one(coeffs.omega_2npm(cnm), coeffs.kx_2npm(cnm), coeffs.ky_2npm(cnm), coeffs.kappa_2npm(cnm), ...
                                    coeffs.F_2npm(cnm), coeffs.A_2npm(cnm), coeffs.B_2npm(cnm));
                            end
                        end
                    end
                end
            end

            c3 = 0;
            trip_super = 0;
            for n = 1:coeffs.N
                for m = n+1:coeffs.N
                    for pmm = [1 -1]
                        for q = m+1:coeffs.N
                            for pmp = [1 -1]
                                c3 = c3 + 1;
                                if superOnly
                                    if pmm == 1 && pmp == 1
                                        trip_super = trip_super + 1;
                                        accumulate_one(coeffs.omega_npmpp(trip_super), coeffs.kx_npmpp(trip_super), coeffs.ky_npmpp(trip_super), coeffs.kappa_npmpp(trip_super), ...
                                            coeffs.F_npmpp(trip_super), coeffs.A_npmpp(trip_super), coeffs.B_npmpp(trip_super));
                                    end
                                else
                                    if isempty(mask_npmpp) || mask_npmpp(c3)
                                        accumulate_one(coeffs.omega_npmpp(c3), coeffs.kx_npmpp(c3), coeffs.ky_npmpp(c3), coeffs.kappa_npmpp(c3), ...
                                            coeffs.F_npmpp(c3), coeffs.A_npmpp(c3), coeffs.B_npmpp(c3));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        function accumulate_one(omega_val, kx_val, ky_val, kappa_val, F_val, A_val, B_val)
                if ~(isfinite(omega_val) && isfinite(kx_val) && isfinite(ky_val) && ...
                     isfinite(kappa_val) && isfinite(F_val) && isfinite(A_val) && isfinite(B_val))
                    return;
                end

                theta = omega_val*t - kx_val*X - ky_val*Y;
                factorZ = F_val * cosh(kappa_val*Z);
                phiAdd = factorZ .* (A_val*sin(theta) - B_val*cos(theta));
                uAdd = kx_val * factorZ .* (-A_val*cos(theta) - B_val*sin(theta));
                vAdd = ky_val * factorZ .* (-A_val*cos(theta) - B_val*sin(theta));
                wAdd = F_val * sinh(kappa_val*Z) * kappa_val .* ...
                    (A_val*sin(theta) - B_val*cos(theta));
                pAdd = -factorZ * omega_val .* (A_val*cos(theta) + B_val*sin(theta));

                V = sqrt(uAdd.^2 + vAdd.^2);
                uV_loc = uV_loc + uAdd .* V;
                vV_loc = vV_loc + vAdd .* V;

                u_loc = u_loc + uAdd;
                v_loc = v_loc + vAdd;
                w_loc = w_loc + wAdd;

                omega_safe = omega_val;
                if omega_safe == 0
                    omega_safe = eps;
                end
                T = tanh(kappa_val*Z) .* pAdd;
                u_t_loc = kx_val * omega_val .* phiAdd;
                v_t_loc = ky_val * omega_val .* phiAdd;
                u_x_loc = -(kx_val^2) .* phiAdd;
                v_x_loc = -(kx_val*ky_val) .* phiAdd;
                u_y_loc = -(kx_val*ky_val) .* phiAdd;
                v_y_loc = -(ky_val^2) .* phiAdd;
                u_z_loc = kx_val * kappa_val / omega_safe .* T;
                v_z_loc = ky_val * kappa_val / omega_safe .* T;

                a_x_loc = a_x_loc + u_t_loc + u_loc.*u_x_loc + v_loc.*u_y_loc + w_loc.*u_z_loc;
                a_y_loc = a_y_loc + v_t_loc + u_loc.*v_x_loc + v_loc.*v_y_loc + w_loc.*v_z_loc;
        end
    end

    function [isWG, isNearRes] = detect_wavegroup_or_resonance(c)
        isWG = false;
        isNearRes = false;

        try
            if isfield(c, 'kappa') && ~isempty(c.kappa) && isfield(c, 'kx') && isfield(c, 'ky')
                kap = c.kappa(:);
                kap = kap(isfinite(kap) & kap > 0);
                if ~isempty(kap)
                    cv_k = std(kap) / max(mean(kap), eps);
                    th = atan2(c.ky(:), c.kx(:));
                    R = abs(mean(exp(1i*th)));
                    circ_std = sqrt(max(0, -2*log(max(R, eps))));
                    isWG = (cv_k < 0.22) && (circ_std < deg2rad(20));
                end
            end

            if isfield(c, 'omega') && ~isempty(c.omega)
                w_ref = median(abs(c.omega(:)));
                if ~isfinite(w_ref) || w_ref <= 0
                    w_ref = 1;
                end
                w_tol = max(1e-10, 1e-3*w_ref);

                resFlags = false(3,1);
                if isfield(c, 'omega_np2m') && numel(c.omega_np2m) >= 2
                    wv = c.omega_np2m(2:2:end);
                    resFlags(1) = any(isfinite(wv) & abs(wv) < w_tol);
                end
                if isfield(c, 'omega_2npm') && numel(c.omega_2npm) >= 2
                    wv = c.omega_2npm(2:2:end);
                    resFlags(2) = any(isfinite(wv) & abs(wv) < w_tol);
                end
                if isfield(c, 'omega_npmpp') && isfield(c, 'N')
                    wv = c.omega_npmpp(:);
                    if ~isempty(wv) && isscalar(c.N) && c.N >= 3
                        idx = 0;
                        subMask = false(numel(wv),1);
                        for n = 1:c.N
                            for m = n+1:c.N
                                for pmm = [1 -1]
                                    for pidx = m+1:c.N
                                        for pmp = [1 -1]
                                            idx = idx + 1;
                                            if ~(pmm == 1 && pmp == 1)
                                                subMask(idx) = true;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        ww = wv(subMask);
                        resFlags(3) = any(isfinite(ww) & abs(ww) < w_tol);
                    end
                end
                isNearRes = any(resFlags);
            end
        catch
            isWG = false;
            isNearRes = false;
        end
    end
end
