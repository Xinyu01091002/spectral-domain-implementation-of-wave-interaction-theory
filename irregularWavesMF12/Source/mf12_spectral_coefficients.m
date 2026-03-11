function [coeffs] = mf12_spectral_coefficients(order,g,h,a,b,kx,ky,Ux,Uy,varargin)
%MF12_SPECTRAL_COEFFICIENTS MF12 coefficients with selective interaction retention.
% This is the preferred spectral-coefficient implementation name in this repository.
%   Keeps:
%   - Second-order self, sum, and difference interactions.
%   - Third-order first-order corrections (dispersion and potential corrections).
%   - Third-order superharmonics only: (n+2m), (2n+m), and (n+m+p).
%
%   Inputs match mf12_direct_coefficients.
%
%   Optional trailing inputs:
%     dispCoeffs (scalar, default 0)
%     opts (struct, optional):
%       opts.enable_subharmonic (default false)
%   Notes:
%     Disk streaming/out-of-core paths are intentionally disabled in this
%     single-workflow version; this function always computes in RAM.

dispCoeffs = 0;
opts = struct();
if ~isempty(varargin)
    if isnumeric(varargin{1}) || islogical(varargin{1})
        dispCoeffs = varargin{1};
        if numel(varargin) >= 2 && isstruct(varargin{2})
            opts = varargin{2};
        end
    elseif isstruct(varargin{1})
        opts = varargin{1};
    end
end

if ~isfield(opts, 'enable_subharmonic'), opts.enable_subharmonic = false; end
if opts.enable_subharmonic
    % Isolated full-subharmonic path. No subharmonic memory/compute in default mode.
    coeffs = mf12_direct_coefficients(order,g,h,a,b,kx,ky,Ux,Uy,dispCoeffs);
    coeffs.superharmonic_only = false;
    coeffs.enable_subharmonic = true;
    coeffs.opts = opts;
    return;
end

% Ensure row vectors for consistent indexing/storage.
a = a(:).';
b = b(:).';
kx = kx(:).';
ky = ky(:).';

N = numel(a);
muStar = zeros(1,N);

% Linear quantities
kappa = sqrt(kx.^2 + ky.^2);                 % Eq. 3.4
omega1 = sqrt(g*kappa.*tanh(h*kappa));       % Eq. 3.5b
omega = kx*Ux + ky*Uy + omega1;              % Eq. 3.5a
F = -omega1./(kappa.*sinh(h*kappa));         % Eq. 3.6
mu = F.*cosh(h*kappa);                       % Eq. 3.77
kappa_2 = 2*kappa;
kx_2 = 2*kx;
ky_2 = 2*ky;
c = sqrt(a.^2 + b.^2);

% Second-order self terms
if order >= 2
    G_2 = 0.5*h*kappa.*(2 + cosh(2*h*kappa)).*coth(h*kappa)./(sinh(h*kappa).^2); % Eq. 3.26a
    F_2 = -0.75*h*omega1./(sinh(h*kappa).^4);                                     % Eq. 3.26b
    A_2 = (a.^2 - b.^2)/(2*h);                                                    % Eq. 3.11
    B_2 = (a.*b)/h;                                                                % Eq. 3.11
    mu_2 = F_2.*cosh(h*kappa_2) - h*omega1;                                       % Eq. 3.79

    % Mass flux coefficient retained for consistency with the direct coefficient structure.
    M = c.^2.*omega1./(2*kappa).*coth(h*kappa);                                   % Eq. 3.70 factor

    % Second-order pair interactions: keep both n+m and n-m
    numPairs = N*(N-1)/2;
    len2 = 2*numPairs;

    omega_npm = zeros(1,len2);
    kx_npm = zeros(1,len2); ky_npm = zeros(1,len2);
    kappa_npm = zeros(1,len2);
    alpha_npm = zeros(1,len2);
    gamma_npm = zeros(1,len2);
    beta_npm = zeros(1,len2);
    F_npm = zeros(1,len2);
    G_npm = zeros(1,len2);
    mu_npm = zeros(1,len2);
    A_npm = zeros(1,len2);
    B_npm = zeros(1,len2);

    pairCount = 0;
    for n = 1:N
        for m = (n+1):N
            pairCount = pairCount + 1;
            idxPlus = 2*pairCount - 1;
            idxMinus = idxPlus + 1;

            [omega_npm(idxPlus), kx_npm(idxPlus), ky_npm(idxPlus), kappa_npm(idxPlus), ...
                alpha_npm(idxPlus), gamma_npm(idxPlus), beta_npm(idxPlus), ...
                F_npm(idxPlus), G_npm(idxPlus), mu_npm(idxPlus)] = pair_terms(1,n,m);

            [omega_npm(idxMinus), kx_npm(idxMinus), ky_npm(idxMinus), kappa_npm(idxMinus), ...
                alpha_npm(idxMinus), gamma_npm(idxMinus), beta_npm(idxMinus), ...
                F_npm(idxMinus), G_npm(idxMinus), mu_npm(idxMinus)] = pair_terms(-1,n,m);

            A_npm(idxPlus) = (a(n)*a(m) - b(n)*b(m))/h;    % Eq. 3.10a, pm=+1
            B_npm(idxPlus) = (a(m)*b(n) + a(n)*b(m))/h;    % Eq. 3.10b, pm=+1
            A_npm(idxMinus) = (a(n)*a(m) + b(n)*b(m))/h;   % Eq. 3.10a, pm=-1
            B_npm(idxMinus) = (a(m)*b(n) - a(n)*b(m))/h;   % Eq. 3.10b, pm=-1
        end
    end
end

% Third-order terms
if order == 3
    % Third-order corrections to first-order quantities.
    Upsilon = omega1.*kappa.*(-13 + 24*cosh(2*h*kappa) + cosh(4*h*kappa))./(64*sinh(h*kappa).^5); % Eq. 3.67
    Xi = (omega1.*G_2 + F_2.*kappa_2.*sinh(h*kappa_2) - g*h*kappa.^2./(2*omega1))/(4*h);           % Eq. 3.85
    F13 = c.^2.*Upsilon;                                                                             % Eq. 3.66 (part)
    muStar = c.^2.*Xi;                                                                               % Eq. 3.84 (part)
    Omega = (8 + cosh(4*h*kappa))./(16*sinh(h*kappa).^4);                                           % Eq. 3.74
    omega3 = c.^2.*kappa.^2.*Omega;                                                                  % Eq. 3.73 (part)

    for n = 1:N
        for m = [1:n-1, n+1:N]
            % + interaction
            [omega_np, kx_np, ky_np, kappa_np, alpha_np, gamma_np, beta_np, F_np, G_np, ~] = pair_terms(1,n,m);
            % - interaction
            [omega_nm, kx_nm, ky_nm, kappa_nm, alpha_nm, gamma_nm, beta_nm, F_nm, G_nm, ~] = pair_terms(-1,n,m);

            ups = Upsilon_nm(omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), ...
                F_np,F_nm,G_np,G_nm, kappa_np,kappa_nm, g,h);
            xi = Xi_nm(omega1(n),kappa(n), omega1(m), G_np,F_np,gamma_np, G_nm,F_nm,gamma_nm, h,g);
            om = Omega_nm(omega1(n),kx(n),ky(n), omega1(m),kx(m),ky(m),kappa(m), F_np,F_nm,G_np,G_nm, kappa_np,kappa_nm, g,h);

            F13(n) = F13(n) + c(m)^2*ups;
            muStar(n) = muStar(n) + c(m)^2*xi;
            omega3(n) = omega3(n) + c(m)^2*kappa(m)^2*om;

        end
    end

    omega = omega + omega3.*omega1;              % Eq. 3.72
    muStar = muStar + F13.*cosh(h*kappa);        % Eq. 3.84 correction

    % Third-order single summations (self-self-self harmonic)
    kappa_3 = 3*kappa;
    gamma_2 = kappa_2.*sinh(h*kappa_2);

    A_3 = zeros(1,N);
    B_3 = zeros(1,N);
    F_3 = zeros(1,N);
    G_3 = zeros(1,N);
    mu_3 = zeros(1,N);

    for n = 1:N
        A_3(n) = 0.5*ThetaA(a(n),b(n), a(n),b(n), a(n),b(n), h); % Eq. 3.38
        B_3(n) = 0.5*ThetaB(a(n),b(n), a(n),b(n), a(n),b(n), h); % Eq. 3.39
        F_3(n) = (h^2*kappa(n)*omega1(n)/(32*sinh(h*kappa(n))^7))*(-11 + 2*cosh(2*h*kappa(n))); % Eq. 3.65
        G_3(n) = (3*h^2*kappa(n)^2/(128*sinh(h*kappa(n))^6))*(14 + 15*cosh(2*h*kappa(n)) + 6*cosh(4*h*kappa(n)) + cosh(6*h*kappa(n))); % Eq. 3.64
        mu_3(n) = F_3(n)*cosh(h*kappa_3(n)) - g*h^2*kappa(n)^2/(4*omega1(n)) + 0.5*h*(F_2(n)*gamma_2(n) - omega1(n)*G_2(n)); % Eq. 3.80
    end

    % Two pair-index maps are needed to mirror the retained direct-index layout:
    % 1) Row-wise odd map: matches cnm index for pmm=+1 in triple loops.
    M_nm_row_odd = zeros(N);
    pairCount = 0;
    for n = 1:N
        for m = (n+1):N
            pairCount = pairCount + 1;
            M_nm_row_odd(n,m) = 2*pairCount - 1;
        end
    end
    % 2) Column-major odd map: matches the historical M_nm construction.
    M_nm_col_odd = zeros(N);
    nm_indices = 1:(N*(N-1)/2);
    M_nm_col_odd(triu(ones(N),1)==1) = nm_indices;
    M_nm_col_odd = 2*M_nm_col_odd - 1;

    numPairs = N*(N-1)/2;
    omega_np2m = zeros(1,numPairs);
    kx_np2m = zeros(1,numPairs); ky_np2m = zeros(1,numPairs);
    kappa_np2m = zeros(1,numPairs); alpha_np2m = zeros(1,numPairs);
    gamma_np2m = zeros(1,numPairs); beta_np2m = zeros(1,numPairs);
    F_np2m = zeros(1,numPairs); G_np2m = zeros(1,numPairs); mu_np2m = zeros(1,numPairs);
    A_np2m = zeros(1,numPairs); B_np2m = zeros(1,numPairs);

    omega_2npm = zeros(1,numPairs);
    kx_2npm = zeros(1,numPairs); ky_2npm = zeros(1,numPairs);
    kappa_2npm = zeros(1,numPairs); alpha_2npm = zeros(1,numPairs);
    gamma_2npm = zeros(1,numPairs); beta_2npm = zeros(1,numPairs);
    F_2npm = zeros(1,numPairs); G_2npm = zeros(1,numPairs); mu_2npm = zeros(1,numPairs);
    A_2npm = zeros(1,numPairs); B_2npm = zeros(1,numPairs);

    pairCount = 0;
    for n = 1:N
        for m = (n+1):N
            pairCount = pairCount + 1;
            idxSum_nm = 2*pairCount - 1; % n+m entry in second-order arrays

            % n+2m
            omega_np2m(pairCount) = omega1(n) + 2*omega1(m);                    % Eq. 3.44b (pm=+1)
            kx_np2m(pairCount) = kx(n) + 2*kx(m);
            ky_np2m(pairCount) = ky(n) + 2*ky(m);
            kappa_np2m(pairCount) = hypot(kx_np2m(pairCount), ky_np2m(pairCount));
            alpha_np2m(pairCount) = omega_np2m(pairCount)*cosh(h*kappa_np2m(pairCount));
            gamma_np2m(pairCount) = kappa_np2m(pairCount)*sinh(h*kappa_np2m(pairCount));
            beta_np2m(pairCount) = omega_np2m(pairCount)^2*cosh(h*kappa_np2m(pairCount)) - g*kappa_np2m(pairCount)*sinh(h*kappa_np2m(pairCount));

            A_np2m(pairCount) = 0.5*ThetaA(a(n),b(n), a(m),b(m), a(m),b(m), h);  % Eq. 3.36, pm=+1
            B_np2m(pairCount) = 0.5*ThetaB(a(n),b(n), a(m),b(m), a(m),b(m), h);

            G_np2m(pairCount) = Lambda3(omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), omega1(m),kx(m),ky(m),kappa(m), ...
                kappa_npm(idxSum_nm),gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                kappa_npm(idxSum_nm),gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                kappa_2(m),gamma_2(m),G_2(m),F_2(m), ...
                omega_np2m(pairCount),alpha_np2m(pairCount),gamma_np2m(pairCount),beta_np2m(pairCount), g,h);

            F_np2m(pairCount) = Gamma3(omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), omega1(m),kx(m),ky(m),kappa(m), ...
                kappa_npm(idxSum_nm),gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                kappa_npm(idxSum_nm),gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                kappa_2(m),gamma_2(m),G_2(m),F_2(m), ...
                omega_np2m(pairCount),beta_np2m(pairCount), g,h);

            mu_np2m(pairCount) = Pi(omega1(n),kappa(n), omega1(m),kappa(m), omega1(m),kappa(m), ...
                gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                gamma_2(m),G_2(m),F_2(m), F_np2m(pairCount),kappa_np2m(pairCount), g,h);

            % 2n+m
            omega_2npm(pairCount) = 2*omega1(n) + omega1(m);                    % Eq. 3.44c (pm=+1)
            kx_2npm(pairCount) = 2*kx(n) + kx(m);
            ky_2npm(pairCount) = 2*ky(n) + ky(m);
            kappa_2npm(pairCount) = hypot(kx_2npm(pairCount), ky_2npm(pairCount));
            alpha_2npm(pairCount) = omega_2npm(pairCount)*cosh(h*kappa_2npm(pairCount));
            gamma_2npm(pairCount) = kappa_2npm(pairCount)*sinh(h*kappa_2npm(pairCount));
            beta_2npm(pairCount) = omega_2npm(pairCount)^2*cosh(h*kappa_2npm(pairCount)) - g*kappa_2npm(pairCount)*sinh(h*kappa_2npm(pairCount));

            A_2npm(pairCount) = 0.5*ThetaA(a(n),b(n), a(n),b(n), a(m),b(m), h);  % From Eq. 3.32/3.34
            B_2npm(pairCount) = 0.5*ThetaB(a(n),b(n), a(n),b(n), a(m),b(m), h);  % From Eq. 3.33/3.35

            G_2npm(pairCount) = Lambda3(omega1(n),kx(n),ky(n),kappa(n), omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), ...
                kappa_2(n),gamma_2(n),G_2(n),F_2(n), ...
                kappa_npm(idxSum_nm),gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                kappa_npm(idxSum_nm),gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                omega_2npm(pairCount),alpha_2npm(pairCount),gamma_2npm(pairCount),beta_2npm(pairCount), g,h);

            F_2npm(pairCount) = Gamma3(omega1(n),kx(n),ky(n),kappa(n), omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), ...
                kappa_2(n),gamma_2(n),G_2(n),F_2(n), ...
                kappa_npm(idxSum_nm),gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                kappa_npm(idxSum_nm),gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                omega_2npm(pairCount),beta_2npm(pairCount), g,h);

            mu_2npm(pairCount) = Pi(omega1(n),kappa(n), omega1(n),kappa(n), omega1(m),kappa(m), ...
                gamma_2(n),G_2(n),F_2(n), ...
                gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                F_2npm(pairCount),kappa_2npm(pairCount), g,h);
        end
    end

    % Triple summations: superharmonic n+m+p only
    numTriplets = N*(N-1)*(N-2)/6;
    omega_npmpp = zeros(1,numTriplets);
    kx_npmpp = zeros(1,numTriplets); ky_npmpp = zeros(1,numTriplets);
    kappa_npmpp = zeros(1,numTriplets);
    alpha_npmpp = zeros(1,numTriplets);
    beta_npmpp = zeros(1,numTriplets);
    gamma_npmpp = zeros(1,numTriplets);
    F_npmpp = zeros(1,numTriplets);
    G_npmpp = zeros(1,numTriplets);
    mu_npmpp = zeros(1,numTriplets);
    A_npmpp = zeros(1,numTriplets);
    B_npmpp = zeros(1,numTriplets);

    c3 = 0;
    for n = 1:N
        for m = (n+1):N
            idxSum_nm = M_nm_row_odd(n,m);
            for p = (m+1):N
                c3 = c3 + 1;
                idxSum_np = M_nm_col_odd(n,p);
                idxSum_mp = M_nm_col_odd(m,p);

                omega_npmpp(c3) = omega1(n) + omega1(m) + omega1(p);            % Eq. 3.44a, +++
                kx_npmpp(c3) = kx(n) + kx(m) + kx(p);
                ky_npmpp(c3) = ky(n) + ky(m) + ky(p);
                kappa_npmpp(c3) = hypot(kx_npmpp(c3), ky_npmpp(c3));
                alpha_npmpp(c3) = omega_npmpp(c3)*cosh(h*kappa_npmpp(c3));
                beta_npmpp(c3) = omega_npmpp(c3)^2*cosh(h*kappa_npmpp(c3)) - g*kappa_npmpp(c3)*sinh(h*kappa_npmpp(c3));
                gamma_npmpp(c3) = kappa_npmpp(c3)*sinh(h*kappa_npmpp(c3));

                A_npmpp(c3) = 0.5*ThetaA(a(n),b(n), a(m),b(m), a(p),b(p), h);
                B_npmpp(c3) = 0.5*ThetaB(a(n),b(n), a(m),b(m), a(p),b(p), h);

                G_npmpp(c3) = Lambda3(omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), omega1(p),kx(p),ky(p),kappa(p), ...
                    kappa_npm(idxSum_nm),gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                    kappa_npm(idxSum_np),gamma_npm(idxSum_np),G_npm(idxSum_np),F_npm(idxSum_np), ...
                    kappa_npm(idxSum_mp),gamma_npm(idxSum_mp),G_npm(idxSum_mp),F_npm(idxSum_mp), ...
                    omega_npmpp(c3),alpha_npmpp(c3),gamma_npmpp(c3),beta_npmpp(c3), g,h);

                F_npmpp(c3) = Gamma3(omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), omega1(p),kx(p),ky(p),kappa(p), ...
                    kappa_npm(idxSum_nm),gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                    kappa_npm(idxSum_np),gamma_npm(idxSum_np),G_npm(idxSum_np),F_npm(idxSum_np), ...
                    kappa_npm(idxSum_mp),gamma_npm(idxSum_mp),G_npm(idxSum_mp),F_npm(idxSum_mp), ...
                    omega_npmpp(c3),beta_npmpp(c3), g,h);

                mu_npmpp(c3) = Pi(omega1(n),kappa(n), omega1(m),kappa(m), omega1(p),kappa(p), ...
                    gamma_npm(idxSum_nm),G_npm(idxSum_nm),F_npm(idxSum_nm), ...
                    gamma_npm(idxSum_np),G_npm(idxSum_np),F_npm(idxSum_np), ...
                    gamma_npm(idxSum_mp),G_npm(idxSum_mp),F_npm(idxSum_mp), ...
                    F_npmpp(c3),kappa_npmpp(c3), g,h);
            end
        end
    end
end

% Output coefficients
coeffs.g = g;
coeffs.h = h;
coeffs.N = N;
coeffs.a = a;
coeffs.b = b;
coeffs.kx = kx;
coeffs.ky = ky;
coeffs.Ux = Ux;
coeffs.Uy = Uy;
coeffs.kappa = kappa;
coeffs.omega1 = omega1;
coeffs.omega = omega;
coeffs.mu = mu;
coeffs.muStar = muStar;
coeffs.F = F;
coeffs.c = c;
coeffs.kappa_2 = kappa_2;
coeffs.superharmonic_only = true;

if order >= 2
    coeffs.A_2 = A_2;
    coeffs.B_2 = B_2;
    coeffs.F_2 = F_2;
    coeffs.G_2 = G_2;
    coeffs.mu_2 = mu_2;

    coeffs.F_npm = F_npm;
    coeffs.G_npm = G_npm;
    coeffs.A_npm = A_npm;
    coeffs.B_npm = B_npm;
    coeffs.mu_npm = mu_npm;
    coeffs.kappa_npm = kappa_npm;
    coeffs.omega_npm = omega_npm;

    coeffs.kx_2 = kx_2;
    coeffs.ky_2 = ky_2;
    coeffs.kx_npm = kx_npm;
    coeffs.ky_npm = ky_npm;
    coeffs.M = M;
end

if order >= 3
    omega_npm_corr = omega_npm;
    pairCount = 0;
    for n = 1:N
        for m = (n+1):N
            pairCount = pairCount + 1;
            idxSum = 2*pairCount - 1;
            idxDiff = idxSum + 1;
            omega_npm_corr(idxSum) = omega(n) + omega(m);
            omega_npm_corr(idxDiff) = omega(n) - omega(m);
            omega_np2m(pairCount) = omega(n) + 2*omega(m);
            omega_2npm(pairCount) = 2*omega(n) + omega(m);
        end
    end

    c3corr = 0;
    for n = 1:N
        for m = (n+1):N
            for p = (m+1):N
                c3corr = c3corr + 1;
                omega_npmpp(c3corr) = omega(n) + omega(m) + omega(p);
            end
        end
    end

    coeffs.A_3 = A_3;
    coeffs.B_3 = B_3;
    coeffs.F_3 = F_3;
    coeffs.G_3 = G_3;
    coeffs.mu_3 = mu_3;
    coeffs.kappa_3 = kappa_3;
    coeffs.F13 = F13;

    coeffs.A_np2m = A_np2m;
    coeffs.B_np2m = B_np2m;
    coeffs.F_np2m = F_np2m;
    coeffs.G_np2m = G_np2m;
    coeffs.mu_np2m = mu_np2m;
    coeffs.kappa_np2m = kappa_np2m;

    coeffs.A_2npm = A_2npm;
    coeffs.B_2npm = B_2npm;
    coeffs.F_2npm = F_2npm;
    coeffs.G_2npm = G_2npm;
    coeffs.mu_2npm = mu_2npm;
    coeffs.kappa_2npm = kappa_2npm;

    coeffs.omega_np2m = omega_np2m;
    coeffs.omega_2npm = omega_2npm;

    coeffs.A_npmpp = A_npmpp;
    coeffs.B_npmpp = B_npmpp;
    coeffs.F_npmpp = F_npmpp;
    coeffs.G_npmpp = G_npmpp;
    coeffs.mu_npmpp = mu_npmpp;
    coeffs.kappa_npmpp = kappa_npmpp;
    coeffs.omega_npmpp = omega_npmpp;

    coeffs.kx_np2m = kx_np2m;
    coeffs.ky_np2m = ky_np2m;
    coeffs.kx_2npm = kx_2npm;
    coeffs.ky_2npm = ky_2npm;
    coeffs.kx_npmpp = kx_npmpp;
    coeffs.ky_npmpp = ky_npmpp;
    coeffs.gamma_2 = gamma_2;
end

if order >= 2
    coeffs.omega_npm = omega_npm;
end
if order >= 3
    coeffs.omega_npm = omega_npm_corr;
end

if dispCoeffs == 1
    disp(' ');
    disp('Transfer function coefficients:');
    if order >= 2
        disp(['G_2n = ' num2str(G_2)]);
        disp(['F_2n = ' num2str(F_2)]);
        disp(['mu_2n = ' num2str(mu_2)]);
        disp(['G_npm = ' num2str(G_npm)]);
        disp(['F_npm = ' num2str(F_npm)]);
        disp(['mu_npm = ' num2str(mu_npm)]);
    end
    if order == 3
        disp(['G_3n = ' num2str(G_3)]);
        disp(['F_3n = ' num2str(F_3)]);
        disp(['mu_3n = ' num2str(mu_3)]);
        disp(['G_npmpp = ' num2str(G_npmpp)]);
        disp(['F_npmpp = ' num2str(F_npmpp)]);
        disp(['mu_npmpp = ' num2str(mu_npmpp)]);
        disp(['G_np2m = ' num2str(G_np2m)]);
        disp(['F_np2m = ' num2str(F_np2m)]);
        disp(['mu_np2m = ' num2str(mu_np2m)]);
        disp(['G_2npm = ' num2str(G_2npm)]);
        disp(['F_2npm = ' num2str(F_2npm)]);
        disp(['mu_2npm = ' num2str(mu_2npm)]);
    end
end

    function [omega_out,kx_out,ky_out,kappa_out,alpha_out,gamma_out,beta_out,F_out,G_out,mu_out] = pair_terms(pm,n,m)
        omega_out = omega1(n) + pm*omega1(m);
        kx_out = kx(n) + pm*kx(m);
        ky_out = ky(n) + pm*ky(m);
        kappa_out = hypot(kx_out, ky_out);
        alpha_out = omega_out*cosh(h*kappa_out);
        gamma_out = kappa_out*sinh(h*kappa_out);
        beta_out = omega_out^2*cosh(h*kappa_out) - g*kappa_out*sinh(h*kappa_out);

        F_out = Gamma2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
            omega_out,beta_out, g,h);
        G_out = Lambda2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
            omega_out,alpha_out,gamma_out,beta_out, g,h);
        mu_out = F_out*cosh(h*kappa_out) - 0.5*h*(omega1(n) + pm*omega1(m));
    end
end

% Lambda2 function, Eq. 3.18, p. 310
function out = Lambda2(omega1n,knx,kny,kappan, omega1m,kmx,kmy,kappam, omega_npm,alpha_npm,gamma_npm,beta_npm,g,h)
    knkm = knx*kmx + kny*kmy;
    out = h/(2*omega1n*omega1m*beta_npm)*(g*alpha_npm*(omega1n*(kappam^2 + knkm) ...
                                                     + omega1m*(kappan^2 + knkm)) ...
        + gamma_npm*(g^2*knkm + omega1n^2*omega1m^2 - omega1n*omega1m*omega_npm^2));
end

% Gamma2 function, Eq. 3.21, p. 310
function out = Gamma2(omega1n,knx,kny,kappan, omega1m,kmx,kmy,kappam, omega_npm,beta_npm,g,h)
    knkm = knx*kmx + kny*kmy;
    out = h/(2*omega1n*omega1m*beta_npm)*(omega1n*omega1m*omega_npm*(omega_npm^2 - omega1n*omega1m) ...
        - g^2*omega1n*(kappam^2 + 2*knkm) - g^2*omega1m*(kappan^2 + 2*knkm));
end

% Upsilon_nm function, Eq. 3.68, p. 316
function out = Upsilon_nm(omega1n,knx,kny,kappan, omega1m, kmx,kmy,kappam, Fnpm,Fnmm,Gnpm,Gnmm, kappanpm,kappanmm, g,h)
    knkm = knx*kmx + kny*kmy;
    out = g/(4*omega1n*omega1m*cosh(h*kappan))*(omega1m*(kappan^2 - kappam^2) - omega1n*knkm) ...
        + (Gnpm + Gnmm)/(4*h*omega1n^2*omega1m*cosh(h*kappan))*(g^2*knkm + omega1m^3*omega1n) ...
        - 1/(4*h*cosh(h*kappan))*(Fnpm*kappanpm*sinh(h*kappanpm) + Fnmm*kappanmm*sinh(h*kappanmm)) ...
        + g*Fnpm*cosh(h*kappanpm)/(4*h*omega1n^2*omega1m*cosh(h*kappan))*((omega1n + omega1m)*(knkm + kappam^2) - omega1m*kappanpm^2) ...
        + g*Fnmm*cosh(h*kappanmm)/(4*h*omega1n^2*omega1m*cosh(h*kappan))*((omega1n - omega1m)*(knkm - kappam^2) - omega1m*kappanmm^2);
end

% Xi_nm function, Eq. 3.86, p. 319
function out = Xi_nm(omega1n,kappan, omega1m, Gnpm,Fnpm,gamma_npm, Gnmm,Fnmm,gamma_nmm, h,g)
    out = 1/(2*h)*(omega1m*(Gnpm - Gnmm) + Fnpm*gamma_npm + Fnmm*gamma_nmm - g*h*kappan^2/(2*omega1n));
end

% ThetaA function, Eq. 3.34, p. 312
function out = ThetaA(an,bn, am,bm, ap,bp, h)
    out = (an*am*ap - bn*bm*ap - bn*am*bp - an*bm*bp)/(h^2);
end

% ThetaB function, Eq. 3.35, p. 312
function out = ThetaB(an,bn, am,bm, ap,bp, h)
    out = (bn*am*ap + an*bm*ap + an*am*bp - bn*bm*bp)/(h^2);
end

% Omega_nm function, Eq. 3.75, p. 317
function out = Omega_nm(omega1n,knx,kny, omega1m,kmx,kmy,kappam, Fnpm,Fnmm,Gnpm,Gnmm, kappanpm,kappanmm, g,h)
    knkm = knx*kmx + kny*kmy;
    out = 1/(kappam^2)*((2*omega1m^2 + omega1n^2)/(4*omega1n*omega1m)*knkm + 0.25*kappam^2) ...
        + (Gnpm + Gnmm)/(kappam^2)*(g*knkm/(4*h*omega1n*omega1m) - omega1m^2/(4*g*h)) ...
        + omega1n/(4*g*h*kappam^2)*(Fnpm*kappanpm*sinh(h*kappanpm) + Fnmm*kappanmm*sinh(h*kappanmm)) ...
        - Fnpm*cosh(h*kappanpm)/(4*h*omega1n*omega1m*kappam^2)*((omega1n - omega1m)*(kappam^2 + knkm) + omega1m*kappanpm^2) ...
        + Fnmm*cosh(h*kappanmm)/(4*h*omega1n*omega1m*kappam^2)*((omega1n + omega1m)*(kappam^2 - knkm) - omega1m*kappanmm^2);
end

% Lambda3 function, Eq. 3.53, p. 313-314
function out = Lambda3(omega1n,knx,kny,kappan, omega1m,kmx,kmy,kappam, omega1p,kpx,kpy,kappap, ...
                       kappanpm,gammanpm,Gnpm,Fnpm, kappanpp,gammanpp,Gnpp,Fnpp, kappampp,gammampp,Gmpp,Fmpp, ...
                       omega_npmpp,alpha_npmpp,gamma_npmpp,beta_npmpp, g,h)
    knkm = knx*kmx + kny*kmy;
    knkp = knx*kpx + kny*kpy;
    kmkp = kmx*kpx + kmy*kpy;
    out = h^2/(4*beta_npmpp)*(alpha_npmpp*(omega1n*(knkm + knkp + kappan^2) + ...
        omega1m*(knkm + kmkp + kappam^2) + omega1p*(knkp + kmkp + kappap^2)) ...
        + gamma_npmpp*(g/omega1n*(omega1m*knkm + omega1p*knkp - omega_npmpp*kappan^2) + ...
                       g/omega1m*(omega1n*knkm + omega1p*kmkp - omega_npmpp*kappam^2) + ...
                       g/omega1p*(omega1n*knkp + omega1m*kmkp - omega_npmpp*kappap^2))) ...
        - h*Fnpm/(2*beta_npmpp)*(alpha_npmpp*cosh(h*kappanpm)*(knkp + kmkp + kappanpm^2) + ...
          gamma_npmpp*(g/omega1p*(knkp + kmkp)*cosh(h*kappanpm) - gammanpm*omega_npmpp)) ...
        - h*Fnpp/(2*beta_npmpp)*(alpha_npmpp*cosh(h*kappanpp)*(knkm + kmkp + kappanpp^2) + ...
          gamma_npmpp*(g/omega1m*(knkm + kmkp)*cosh(h*kappanpp) - gammanpp*omega_npmpp)) ...
        - h*Fmpp/(2*beta_npmpp)*(alpha_npmpp*cosh(h*kappampp)*(knkm + knkp +kappampp^2) + ...
          gamma_npmpp*(g/omega1n*(knkm + knkp)*cosh(h*kappampp) - gammampp*omega_npmpp)) ...
        + h*Gnpm/(2*beta_npmpp)*(alpha_npmpp*g/omega1p*(knkp + kmkp + kappap^2) - gamma_npmpp*omega1p^2) ...
        + h*Gnpp/(2*beta_npmpp)*(alpha_npmpp*g/omega1m*(knkm + kmkp + kappam^2) - gamma_npmpp*omega1m^2) ...
        + h*Gmpp/(2*beta_npmpp)*(alpha_npmpp*g/omega1n*(knkm + knkp + kappan^2) - gamma_npmpp*omega1n^2);
end

% Gamma3 function, Eq. 3.56, p. 314-315
function out = Gamma3(omega1n,knx,kny,kappan, omega1m,kmx,kmy,kappam, omega1p,kpx,kpy,kappap, ...
                      kappanpm,gammanpm,Gnpm,Fnpm, kappanpp,gammanpp,Gnpp,Fnpp, kappampp,gammampp,Gmpp,Fmpp, ...
                      omega_npmpp,beta_npmpp, g,h)
    knkm = knx*kmx + kny*kmy;
    knkp = knx*kpx + kny*kpy;
    kmkp = kmx*kpx + kmy*kpy;
    out = -g*h^2/(4*beta_npmpp)*(omega1n*(knkm + knkp + kappan^2) ...
        + omega1m*(knkm + kmkp + kappam^2) + omega1p*(knkp + kmkp + kappap^2) ...
        + omega_npmpp/omega1n*(omega1m*knkm + omega1p*knkp - omega_npmpp*kappan^2) ...
        + omega_npmpp/omega1m*(omega1n*knkm + omega1p*kmkp - omega_npmpp*kappam^2) ...
        + omega_npmpp/omega1p*(omega1n*knkp + omega1m*kmkp - omega_npmpp*kappap^2)) ...
        + h*Fnpm/(2*beta_npmpp)*(g*cosh(h*kappanpm)*((knkp + kmkp + kappanpm^2) + ...
            omega_npmpp/omega1p*(knkp + kmkp)) - gammanpm*omega_npmpp^2) ...
        + h*Fnpp/(2*beta_npmpp)*(g*cosh(h*kappanpp)*((knkm + kmkp + kappanpp^2) + ...
            omega_npmpp/omega1m*(knkm + kmkp)) - gammanpp*omega_npmpp^2) ...
        + h*Fmpp/(2*beta_npmpp)*(g*cosh(h*kappampp)*((knkm + knkp + kappampp^2) + ...
            omega_npmpp/omega1n*(knkm + knkp)) - gammampp*omega_npmpp^2) ...
        + h*Gnpm/(2*beta_npmpp)*(omega1p^2*omega_npmpp - g^2/omega1p*(knkp + kmkp + kappap^2)) ...
        + h*Gnpp/(2*beta_npmpp)*(omega1m^2*omega_npmpp - g^2/omega1m*(knkm + kmkp + kappam^2)) ...
        + h*Gmpp/(2*beta_npmpp)*(omega1n^2*omega_npmpp - g^2/omega1n*(knkm + knkp + kappan^2));
end

% Pi function, Eq. 3.81, p. 318
function out = Pi(omega1n,kappan, omega1m,kappam, omega1p,kappap, ...
    gamma_npm,Gnpm,Fnpm, gamma_npp,Gnpp,Fnpp, gamma_mpp,Gmpp,Fmpp, Fnpmpp,kappa_npmpp, g,h)
    out = Fnpmpp*cosh(h*kappa_npmpp) ...
        - g*h^2/4*(kappan^2/omega1n + kappam^2/omega1m + kappap^2/omega1p) ...
        - h/2*(omega1n*Gmpp + omega1m*Gnpp + omega1p*Gnpm) ...
        + h/2*(Fnpm*gamma_npm + Fnpp*gamma_npp + Fmpp*gamma_mpp);
end
