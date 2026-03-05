function [eta, phiS, X, Y] = surfaceMF12_streaming_super(order, g, h, a, b, kx, ky, Ux, Uy, Lx, Ly, Nx, Ny, t)
%SURFACEMF12_STREAMING_SUPER Streaming MF12 reconstruction (2nd +/- and 3rd +++).
%   This implementation avoids storing large O(N^3) coefficient arrays by
%   accumulating directly into spectral bins.

% Normalize vector shapes.
a = a(:).';
b = b(:).';
kx = kx(:).';
ky = ky(:).';
N = numel(a);

% Grid and spectral spacing.
dx = Lx / Nx;
dy = Ly / Ny;
x_axis = (0:Nx-1) * dx;
y_axis = (0:Ny-1) * dy;
[X, Y] = meshgrid(x_axis, y_axis);

dkx = 2*pi/Lx;
dky = 2*pi/Ly;

spec_eta = complex(zeros(Ny, Nx));
spec_phi = complex(zeros(Ny, Nx));
use_scalar_fast = ~strcmp(getenv('MF12_DISABLE_SCALAR_FAST'), '1');

% Linear quantities.
kappa = hypot(kx, ky);                                % Eq. 3.4
omega1 = sqrt(g*kappa.*tanh(h*kappa));               % Eq. 3.5b
omega = kx*Ux + ky*Uy + omega1;                      % Eq. 3.5a
F = -omega1./(kappa.*sinh(h*kappa));                 % Eq. 3.6
mu = F.*cosh(h*kappa);                               % Eq. 3.77
kappa_2 = 2*kappa;
c = hypot(a, b);
muStar = zeros(1, N);

if order >= 2
    G_2 = 0.5*h*kappa.*(2 + cosh(2*h*kappa)).*coth(h*kappa)./(sinh(h*kappa).^2);  % Eq. 3.26a
    F_2 = -0.75*h*omega1./(sinh(h*kappa).^4);                                      % Eq. 3.26b
    A_2 = (a.^2 - b.^2)/(2*h);                                                     % Eq. 3.11
    B_2 = (a.*b)/h;                                                                 % Eq. 3.11
    mu_2 = F_2.*cosh(h*kappa_2) - h*omega1;                                        % Eq. 3.79
end

if order == 3
    % Third-order corrections to first-order quantities.
    Upsilon = omega1.*kappa.*(-13 + 24*cosh(2*h*kappa) + cosh(4*h*kappa))./(64*sinh(h*kappa).^5); % Eq. 3.67
    Xi = (omega1.*G_2 + F_2.*kappa_2.*sinh(h*kappa_2) - g*h*kappa.^2./(2*omega1))/(4*h);          % Eq. 3.85
    F13 = c.^2.*Upsilon;                                                                            % Eq. 3.66 (part)
    muStar = c.^2.*Xi;                                                                              % Eq. 3.84 (part)
    Omega = (8 + cosh(4*h*kappa))./(16*sinh(h*kappa).^4);                                          % Eq. 3.74
    omega3 = c.^2.*kappa.^2.*Omega;                                                                 % Eq. 3.73 (part)

    for n = 1:N
        for m = [1:n-1, n+1:N]
            [~, ~, ~, kappa_np, ~, gamma_np, ~, F_np, G_np, ~] = pair_terms(1, n, m);
            [~, ~, ~, kappa_nm, ~, gamma_nm, ~, F_nm, G_nm, ~] = pair_terms(-1, n, m);

            ups = Upsilon_nm(omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), ...
                F_np,F_nm,G_np,G_nm, kappa_np,kappa_nm, g,h);
            xi = Xi_nm(omega1(n),kappa(n), omega1(m), G_np,F_np,gamma_np, G_nm,F_nm,gamma_nm, h,g);
            om = Omega_nm(omega1(n),kx(n),ky(n), omega1(m),kx(m),ky(m),kappa(m), F_np,F_nm,G_np,G_nm, kappa_np,kappa_nm, g,h);

            F13(n) = F13(n) + c(m)^2*ups;
            muStar(n) = muStar(n) + c(m)^2*xi;
            omega3(n) = omega3(n) + c(m)^2*kappa(m)^2*om;
        end
    end

    omega = omega + omega3.*omega1;               % Eq. 3.72
    muStar = muStar + F13.*cosh(h*kappa);         % Eq. 3.84 correction
end

% 1) Linear terms.
Z_lin = (a + 1i*b) .* exp(-1i * omega * t);
add_to_spec(kx, ky, Z_lin, 'eta');
add_to_spec(kx, ky, Z_lin .* (mu + muStar) * (1i), 'phi');

% 2) Second-order terms.
if order >= 2
    Z_2 = (A_2 + 1i*B_2) .* exp(-1i * (2*omega) * t);
    add_to_spec(2*kx, 2*ky, Z_2 .* G_2, 'eta');
    add_to_spec(2*kx, 2*ky, Z_2 .* mu_2 * (1i), 'phi');

    for n = 1:N
        for m = (n+1):N
            for pm = [1, -1]
                [omega_npm, kx_npm, ky_npm, ~, ~, ~, ~, ~, G_npm, mu_npm] = pair_terms(pm, n, m);
                A_npm = (a(n)*a(m) - pm*b(n)*b(m))/h;  % Eq. 3.10a
                B_npm = (a(m)*b(n) + pm*a(n)*b(m))/h;  % Eq. 3.10b
                Z_npm = (A_npm + 1i*B_npm) * exp(-1i * omega_npm * t);
                add_to_spec(kx_npm, ky_npm, Z_npm * G_npm, 'eta');
                add_to_spec(kx_npm, ky_npm, Z_npm * mu_npm * (1i), 'phi');
            end
        end
    end
end

% 3) Third-order superharmonic terms only.
if order == 3
    % Single summations.
    kappa_3 = 3*kappa;
    gamma_2 = kappa_2.*sinh(h*kappa_2);
    A_3 = zeros(1, N);
    B_3 = zeros(1, N);
    F_3 = zeros(1, N);
    G_3 = zeros(1, N);
    mu_3 = zeros(1, N);
    for n = 1:N
        A_3(n) = 0.5*ThetaA(a(n),b(n), a(n),b(n), a(n),b(n), h); % Eq. 3.38
        B_3(n) = 0.5*ThetaB(a(n),b(n), a(n),b(n), a(n),b(n), h); % Eq. 3.39
        F_3(n) = (h^2*kappa(n)*omega1(n)/(32*sinh(h*kappa(n))^7))*(-11 + 2*cosh(2*h*kappa(n))); % Eq. 3.65
        G_3(n) = (3*h^2*kappa(n)^2/(128*sinh(h*kappa(n))^6))*(14 + 15*cosh(2*h*kappa(n)) + 6*cosh(4*h*kappa(n)) + cosh(6*h*kappa(n))); % Eq. 3.64
        mu_3(n) = F_3(n)*cosh(h*kappa_3(n)) - g*h^2*kappa(n)^2/(4*omega1(n)) + 0.5*h*(F_2(n)*gamma_2(n) - omega1(n)*G_2(n)); % Eq. 3.80
    end

    Z_3 = (A_3 + 1i*B_3) .* exp(-1i * (3*omega) * t);
    add_to_spec(3*kx, 3*ky, Z_3 .* G_3, 'eta');
    add_to_spec(3*kx, 3*ky, Z_3 .* mu_3 * (1i), 'phi');

    % Double summations: n+2m and 2n+m (pm=+1 only).
    for n = 1:N
        for m = (n+1):N
            [kappa_nm_p, ~, gamma_nm_p, ~, F_nm_p, G_nm_p] = pair_plus_pack(n, m);

            omega_np2m = omega1(n) + 2*omega1(m); % Eq. 3.44b, pm=+1
            kx_np2m = kx(n) + 2*kx(m);
            ky_np2m = ky(n) + 2*ky(m);
            kappa_np2m = hypot(kx_np2m, ky_np2m);
            alpha_np2m = omega_np2m*cosh(h*kappa_np2m);
            gamma_np2m = kappa_np2m*sinh(h*kappa_np2m);
            beta_np2m = omega_np2m^2*cosh(h*kappa_np2m) - g*kappa_np2m*sinh(h*kappa_np2m);

            A_np2m = 0.5*ThetaA(a(n),b(n), a(m),b(m), a(m),b(m), h); % Eq. 3.36, pm=+1
            B_np2m = 0.5*ThetaB(a(n),b(n), a(m),b(m), a(m),b(m), h);

            G_np2m = Lambda3(omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), omega1(m),kx(m),ky(m),kappa(m), ...
                kappa_nm_p,gamma_nm_p,G_nm_p,F_nm_p, ...
                kappa_nm_p,gamma_nm_p,G_nm_p,F_nm_p, ...
                kappa_2(m),gamma_2(m),G_2(m),F_2(m), ...
                omega_np2m,alpha_np2m,gamma_np2m,beta_np2m, g,h);

            F_np2m = Gamma3(omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), omega1(m),kx(m),ky(m),kappa(m), ...
                kappa_nm_p,gamma_nm_p,G_nm_p,F_nm_p, ...
                kappa_nm_p,gamma_nm_p,G_nm_p,F_nm_p, ...
                kappa_2(m),gamma_2(m),G_2(m),F_2(m), ...
                omega_np2m,beta_np2m, g,h);

            mu_np2m = Pi(omega1(n),kappa(n), omega1(m),kappa(m), omega1(m),kappa(m), ...
                gamma_nm_p,G_nm_p,F_nm_p, gamma_nm_p,G_nm_p,F_nm_p, gamma_2(m),G_2(m),F_2(m), F_np2m,kappa_np2m, g,h);

            Z_np2m = (A_np2m + 1i*B_np2m) * exp(-1i * omega_np2m * t);
            add_to_spec(kx_np2m, ky_np2m, Z_np2m * G_np2m, 'eta');
            add_to_spec(kx_np2m, ky_np2m, Z_np2m * mu_np2m * (1i), 'phi');

            omega_2npm = 2*omega1(n) + omega1(m); % Eq. 3.44c, pm=+1
            kx_2npm = 2*kx(n) + kx(m);
            ky_2npm = 2*ky(n) + ky(m);
            kappa_2npm = hypot(kx_2npm, ky_2npm);
            alpha_2npm = omega_2npm*cosh(h*kappa_2npm);
            gamma_2npm = kappa_2npm*sinh(h*kappa_2npm);
            beta_2npm = omega_2npm^2*cosh(h*kappa_2npm) - g*kappa_2npm*sinh(h*kappa_2npm);

            A_2npm = 0.5*ThetaA(a(n),b(n), a(n),b(n), a(m),b(m), h);
            B_2npm = 0.5*ThetaB(a(n),b(n), a(n),b(n), a(m),b(m), h);

            G_2npm = Lambda3(omega1(n),kx(n),ky(n),kappa(n), omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), ...
                kappa_2(n),gamma_2(n),G_2(n),F_2(n), ...
                kappa_nm_p,gamma_nm_p,G_nm_p,F_nm_p, ...
                kappa_nm_p,gamma_nm_p,G_nm_p,F_nm_p, ...
                omega_2npm,alpha_2npm,gamma_2npm,beta_2npm, g,h);

            F_2npm = Gamma3(omega1(n),kx(n),ky(n),kappa(n), omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), ...
                kappa_2(n),gamma_2(n),G_2(n),F_2(n), ...
                kappa_nm_p,gamma_nm_p,G_nm_p,F_nm_p, ...
                kappa_nm_p,gamma_nm_p,G_nm_p,F_nm_p, ...
                omega_2npm,beta_2npm, g,h);

            mu_2npm = Pi(omega1(n),kappa(n), omega1(n),kappa(n), omega1(m),kappa(m), ...
                gamma_2(n),G_2(n),F_2(n), gamma_nm_p,G_nm_p,F_nm_p, gamma_nm_p,G_nm_p,F_nm_p, F_2npm,kappa_2npm, g,h);

            Z_2npm = (A_2npm + 1i*B_2npm) * exp(-1i * omega_2npm * t);
            add_to_spec(kx_2npm, ky_2npm, Z_2npm * G_2npm, 'eta');
            add_to_spec(kx_2npm, ky_2npm, Z_2npm * mu_2npm * (1i), 'phi');
        end
    end

    % Triple summations: n+m+p only.
    for n = 1:N
        for m = (n+1):N
            [kappa_nm, ~, gamma_nm, ~, F_nm, G_nm] = pair_plus_pack(n, m);
            for p = (m+1):N
                [kappa_np, ~, gamma_np, ~, F_np, G_np] = pair_plus_pack(n, p);
                [kappa_mp, ~, gamma_mp, ~, F_mp, G_mp] = pair_plus_pack(m, p);

                omega_npmpp = omega1(n) + omega1(m) + omega1(p); % Eq. 3.44a, +++
                kx_npmpp = kx(n) + kx(m) + kx(p);
                ky_npmpp = ky(n) + ky(m) + ky(p);
                kappa_npmpp = hypot(kx_npmpp, ky_npmpp);
                alpha_npmpp = omega_npmpp*cosh(h*kappa_npmpp);
                beta_npmpp = omega_npmpp^2*cosh(h*kappa_npmpp) - g*kappa_npmpp*sinh(h*kappa_npmpp);
                gamma_npmpp = kappa_npmpp*sinh(h*kappa_npmpp);

                A_npmpp = 0.5*ThetaA(a(n),b(n), a(m),b(m), a(p),b(p), h);
                B_npmpp = 0.5*ThetaB(a(n),b(n), a(m),b(m), a(p),b(p), h);

                G_npmpp = Lambda3(omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), omega1(p),kx(p),ky(p),kappa(p), ...
                    kappa_nm,gamma_nm,G_nm,F_nm, ...
                    kappa_np,gamma_np,G_np,F_np, ...
                    kappa_mp,gamma_mp,G_mp,F_mp, ...
                    omega_npmpp,alpha_npmpp,gamma_npmpp,beta_npmpp, g,h);

                F_npmpp = Gamma3(omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), omega1(p),kx(p),ky(p),kappa(p), ...
                    kappa_nm,gamma_nm,G_nm,F_nm, ...
                    kappa_np,gamma_np,G_np,F_np, ...
                    kappa_mp,gamma_mp,G_mp,F_mp, ...
                    omega_npmpp,beta_npmpp, g,h);

                mu_npmpp = Pi(omega1(n),kappa(n), omega1(m),kappa(m), omega1(p),kappa(p), ...
                    gamma_nm,G_nm,F_nm, gamma_np,G_np,F_np, gamma_mp,G_mp,F_mp, F_npmpp,kappa_npmpp, g,h);

                Z_npmpp = 2 * (A_npmpp + 1i*B_npmpp) * exp(-1i * omega_npmpp * t); % Factor 2 matches MF12/new
                add_to_spec(kx_npmpp, ky_npmpp, Z_npmpp * G_npmpp, 'eta');
                add_to_spec(kx_npmpp, ky_npmpp, Z_npmpp * mu_npmpp * (1i), 'phi');
            end
        end
    end
end

eta = real(ifft2(spec_eta)) * (Nx * Ny);
phi_wave = real(ifft2(spec_phi)) * (Nx * Ny);
phiS = phi_wave + Ux*X + Uy*Y;

    function [omega_out,kx_out,ky_out,kappa_out,alpha_out,gamma_out,beta_out,F_out,G_out,mu_out] = pair_terms(pm,n,m)
        omega_out = omega1(n) + pm*omega1(m);
        kx_out = kx(n) + pm*kx(m);
        ky_out = ky(n) + pm*ky(m);
        kappa_out = hypot(kx_out, ky_out);
        alpha_out = omega_out*cosh(h*kappa_out);
        gamma_out = kappa_out*sinh(h*kappa_out);
        beta_out = omega_out^2*cosh(h*kappa_out) - g*kappa_out*sinh(h*kappa_out);
        F_out = Gamma2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), omega_out,beta_out, g,h);
        G_out = Lambda2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), omega_out,alpha_out,gamma_out,beta_out, g,h);
        mu_out = F_out*cosh(h*kappa_out) - 0.5*h*(omega1(n) + pm*omega1(m));
    end

    function [kappa_out, alpha_out, gamma_out, beta_out, F_out, G_out] = pair_plus_pack(n,m)
        [~, ~, ~, kappa_out, alpha_out, gamma_out, beta_out, F_out, G_out, ~] = pair_terms(1, n, m);
    end

    function add_to_spec(kx_in, ky_in, values, which_field)
        ux = (kx_in(:) / dkx);
        uy = (ky_in(:) / dky);
        vals = values(:);

        valid = isfinite(ux) & isfinite(uy) & isfinite(vals);
        ux = ux(valid);
        uy = uy(valid);
        vals = vals(valid);
        if isempty(vals)
            return;
        end

        % Fast scalar path: avoid O(Ng) accumarray for inner-loop scalar calls.
        if use_scalar_fast && numel(vals) == 1
            ix0 = floor(ux);
            iy0 = floor(uy);
            fx = ux - ix0;
            fy = uy - iy0;

            tol = 1e-12;
            if abs(fx) < tol, fx = 0; elseif abs(fx-1) < tol, fx = 1; end
            if abs(fy) < tol, fy = 0; elseif abs(fy-1) < tol, fy = 1; end

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

            v = vals;
            w00 = (1-fx)*(1-fy);
            w10 = fx*(1-fy);
            w01 = (1-fx)*fy;
            w11 = fx*fy;

            if strcmp(which_field, 'eta')
                spec_eta(idx_y00, idx_x00) = spec_eta(idx_y00, idx_x00) + v*w00;
                spec_eta(idx_y10, idx_x10) = spec_eta(idx_y10, idx_x10) + v*w10;
                spec_eta(idx_y01, idx_x01) = spec_eta(idx_y01, idx_x01) + v*w01;
                spec_eta(idx_y11, idx_x11) = spec_eta(idx_y11, idx_x11) + v*w11;
            else
                spec_phi(idx_y00, idx_x00) = spec_phi(idx_y00, idx_x00) + v*w00;
                spec_phi(idx_y10, idx_x10) = spec_phi(idx_y10, idx_x10) + v*w10;
                spec_phi(idx_y01, idx_x01) = spec_phi(idx_y01, idx_x01) + v*w01;
                spec_phi(idx_y11, idx_x11) = spec_phi(idx_y11, idx_x11) + v*w11;
            end
            return;
        end

        ix0 = floor(ux);
        iy0 = floor(uy);
        fx = ux - ix0;
        fy = uy - iy0;

        % Keep exact-grid points exact to avoid tiny floating-point leakage.
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
        if isempty(vals4)
            return;
        end

        lin = sub2ind([Ny, Nx], idx_y, idx_x);
        addv = accumarray(lin, vals4, [Ny*Nx, 1]);
        if strcmp(which_field, 'eta')
            spec_eta = spec_eta + reshape(addv, [Ny, Nx]);
        else
            spec_phi = spec_phi + reshape(addv, [Ny, Nx]);
        end
    end
end

% Second-order internal functions
function out = Lambda2(omega1n,knx,kny,kappan, omega1m,kmx,kmy,kappam, omega_npm,alpha_npm,gamma_npm,beta_npm,g,h)
    knkm = knx*kmx + kny*kmy;
    out = h/(2*omega1n*omega1m*beta_npm)*(g*alpha_npm*(omega1n*(kappam^2 + knkm) + omega1m*(kappan^2 + knkm)) ...
        + gamma_npm*(g^2*knkm + omega1n^2*omega1m^2 - omega1n*omega1m*omega_npm^2));
end

function out = Gamma2(omega1n,knx,kny,kappan, omega1m,kmx,kmy,kappam, omega_npm,beta_npm,g,h)
    knkm = knx*kmx + kny*kmy;
    out = h/(2*omega1n*omega1m*beta_npm)*(omega1n*omega1m*omega_npm*(omega_npm^2 - omega1n*omega1m) ...
        - g^2*omega1n*(kappam^2 + 2*knkm) - g^2*omega1m*(kappan^2 + 2*knkm));
end

% Third-order internal functions
function out = Upsilon_nm(omega1n,knx,kny,kappan, omega1m, kmx,kmy,kappam, Fnpm,Fnmm,Gnpm,Gnmm, kappanpm,kappanmm, g,h)
    knkm = knx*kmx + kny*kmy;
    out = g/(4*omega1n*omega1m*cosh(h*kappan))*(omega1m*(kappan^2 - kappam^2) - omega1n*knkm) ...
        + (Gnpm + Gnmm)/(4*h*omega1n^2*omega1m*cosh(h*kappan))*(g^2*knkm + omega1m^3*omega1n) ...
        - 1/(4*h*cosh(h*kappan))*(Fnpm*kappanpm*sinh(h*kappanpm) + Fnmm*kappanmm*sinh(h*kappanmm)) ...
        + g*Fnpm*cosh(h*kappanpm)/(4*h*omega1n^2*omega1m*cosh(h*kappan))*((omega1n + omega1m)*(knkm + kappam^2) - omega1m*kappanpm^2) ...
        + g*Fnmm*cosh(h*kappanmm)/(4*h*omega1n^2*omega1m*cosh(h*kappan))*((omega1n - omega1m)*(knkm - kappam^2) - omega1m*kappanmm^2);
end

function out = Xi_nm(omega1n,kappan, omega1m, Gnpm,Fnpm,gamma_npm, Gnmm,Fnmm,gamma_nmm, h,g)
    out = 1/(2*h)*(omega1m*(Gnpm - Gnmm) + Fnpm*gamma_npm + Fnmm*gamma_nmm - g*h*kappan^2/(2*omega1n));
end

function out = ThetaA(an,bn, am,bm, ap,bp, h)
    out = (an*am*ap - bn*bm*ap - bn*am*bp - an*bm*bp)/(h^2);
end

function out = ThetaB(an,bn, am,bm, ap,bp, h)
    out = (bn*am*ap + an*bm*ap + an*am*bp - bn*bm*bp)/(h^2);
end

function out = Omega_nm(omega1n,knx,kny, omega1m,kmx,kmy,kappam, Fnpm,Fnmm,Gnpm,Gnmm, kappanpm,kappanmm, g,h)
    knkm = knx*kmx + kny*kmy;
    out = 1/(kappam^2)*((2*omega1m^2 + omega1n^2)/(4*omega1n*omega1m)*knkm + 0.25*kappam^2) ...
        + (Gnpm + Gnmm)/(kappam^2)*(g*knkm/(4*h*omega1n*omega1m) - omega1m^2/(4*g*h)) ...
        + omega1n/(4*g*h*kappam^2)*(Fnpm*kappanpm*sinh(h*kappanpm) + Fnmm*kappanmm*sinh(h*kappanmm)) ...
        - Fnpm*cosh(h*kappanpm)/(4*h*omega1n*omega1m*kappam^2)*((omega1n - omega1m)*(kappam^2 + knkm) + omega1m*kappanpm^2) ...
        + Fnmm*cosh(h*kappanmm)/(4*h*omega1n*omega1m*kappam^2)*((omega1n + omega1m)*(kappam^2 - knkm) - omega1m*kappanmm^2);
end

function out = Lambda3(omega1n,knx,kny,kappan, omega1m,kmx,kmy,kappam, omega1p,kpx,kpy,kappap, ...
                       kappanpm,gammanpm,Gnpm,Fnpm, kappanpp,gammanpp,Gnpp,Fnpp, kappampp,gammampp,Gmpp,Fmpp, ...
                       omega_npmpp,alpha_npmpp,gamma_npmpp,beta_npmpp, g,h)
    knkm = knx*kmx + kny*kmy;
    knkp = knx*kpx + kny*kpy;
    kmkp = kmx*kpx + kmy*kpy;
    out = h^2/(4*beta_npmpp)*(alpha_npmpp*(omega1n*(knkm + knkp + kappan^2) + omega1m*(knkm + kmkp + kappam^2) + omega1p*(knkp + kmkp + kappap^2)) ...
        + gamma_npmpp*(g/omega1n*(omega1m*knkm + omega1p*knkp - omega_npmpp*kappan^2) + ...
                       g/omega1m*(omega1n*knkm + omega1p*kmkp - omega_npmpp*kappam^2) + ...
                       g/omega1p*(omega1n*knkp + omega1m*kmkp - omega_npmpp*kappap^2))) ...
        - h*Fnpm/(2*beta_npmpp)*(alpha_npmpp*cosh(h*kappanpm)*(knkp + kmkp + kappanpm^2) + gamma_npmpp*(g/omega1p*(knkp + kmkp)*cosh(h*kappanpm) - gammanpm*omega_npmpp)) ...
        - h*Fnpp/(2*beta_npmpp)*(alpha_npmpp*cosh(h*kappanpp)*(knkm + kmkp + kappanpp^2) + gamma_npmpp*(g/omega1m*(knkm + kmkp)*cosh(h*kappanpp) - gammanpp*omega_npmpp)) ...
        - h*Fmpp/(2*beta_npmpp)*(alpha_npmpp*cosh(h*kappampp)*(knkm + knkp +kappampp^2) + gamma_npmpp*(g/omega1n*(knkm + knkp)*cosh(h*kappampp) - gammampp*omega_npmpp)) ...
        + h*Gnpm/(2*beta_npmpp)*(alpha_npmpp*g/omega1p*(knkp + kmkp + kappap^2) - gamma_npmpp*omega1p^2) ...
        + h*Gnpp/(2*beta_npmpp)*(alpha_npmpp*g/omega1m*(knkm + kmkp + kappam^2) - gamma_npmpp*omega1m^2) ...
        + h*Gmpp/(2*beta_npmpp)*(alpha_npmpp*g/omega1n*(knkm + knkp + kappan^2) - gamma_npmpp*omega1n^2);
end

function out = Gamma3(omega1n,knx,kny,kappan, omega1m,kmx,kmy,kappam, omega1p,kpx,kpy,kappap, ...
                      kappanpm,gammanpm,Gnpm,Fnpm, kappanpp,gammanpp,Gnpp,Fnpp, kappampp,gammampp,Gmpp,Fmpp, ...
                      omega_npmpp,beta_npmpp, g,h)
    knkm = knx*kmx + kny*kmy;
    knkp = knx*kpx + kny*kpy;
    kmkp = kmx*kpx + kmy*kpy;
    out = -g*h^2/(4*beta_npmpp)*(omega1n*(knkm + knkp + kappan^2) + omega1m*(knkm + kmkp + kappam^2) + omega1p*(knkp + kmkp + kappap^2) ...
        + omega_npmpp/omega1n*(omega1m*knkm + omega1p*knkp - omega_npmpp*kappan^2) ...
        + omega_npmpp/omega1m*(omega1n*knkm + omega1p*kmkp - omega_npmpp*kappam^2) ...
        + omega_npmpp/omega1p*(omega1n*knkp + omega1m*kmkp - omega_npmpp*kappap^2)) ...
        + h*Fnpm/(2*beta_npmpp)*(g*cosh(h*kappanpm)*((knkp + kmkp + kappanpm^2) + omega_npmpp/omega1p*(knkp + kmkp)) - gammanpm*omega_npmpp^2) ...
        + h*Fnpp/(2*beta_npmpp)*(g*cosh(h*kappanpp)*((knkm + kmkp + kappanpp^2) + omega_npmpp/omega1m*(knkm + kmkp)) - gammanpp*omega_npmpp^2) ...
        + h*Fmpp/(2*beta_npmpp)*(g*cosh(h*kappampp)*((knkm + knkp + kappampp^2) + omega_npmpp/omega1n*(knkm + knkp)) - gammampp*omega_npmpp^2) ...
        + h*Gnpm/(2*beta_npmpp)*(omega1p^2*omega_npmpp - g^2/omega1p*(knkp + kmkp + kappap^2)) ...
        + h*Gnpp/(2*beta_npmpp)*(omega1m^2*omega_npmpp - g^2/omega1m*(knkm + kmkp + kappam^2)) ...
        + h*Gmpp/(2*beta_npmpp)*(omega1n^2*omega_npmpp - g^2/omega1n*(knkm + knkp + kappan^2));
end

function out = Pi(omega1n,kappan, omega1m,kappam, omega1p,kappap, ...
    gamma_npm,Gnpm,Fnpm, gamma_npp,Gnpp,Fnpp, gamma_mpp,Gmpp,Fmpp, Fnpmpp,kappa_npmpp, g,h)
    out = Fnpmpp*cosh(h*kappa_npmpp) ...
        - g*h^2/4*(kappan^2/omega1n + kappam^2/omega1m + kappap^2/omega1p) ...
        - h/2*(omega1n*Gmpp + omega1m*Gnpp + omega1p*Gnpm) ...
        + h/2*(Fnpm*gamma_npm + Fnpp*gamma_npp + Fmpp*gamma_mpp);
end
