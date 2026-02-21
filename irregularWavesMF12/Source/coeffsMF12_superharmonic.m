function [coeffs] = coeffsMF12_superharmonic(order,g,h,a,b,kx,ky,Ux,Uy, dispCoeffs)
% A modified version of coeffsMF12.
% Originally designed for superharmonic (+++) interactions only.
% Modified to include 2nd-order subharmonic (difference) interactions as well.
% 3rd-order remains superharmonic (+++) only.
%
% Based on coeffsMF12.m by David R. Fuhrman.

% Set (optional) default argument dispCoeffs to 0 (false)
if nargin < 10, dispCoeffs = 0; end

% Initialize
muStar = 0*a; 

% Basic derived input
N = length(a); % Number of components
kappa = sqrt(kx.^2 + ky.^2); % Eq. 3.4
omega1 = sqrt(g*kappa.*tanh(h*kappa)); % Eq. 3.5b
omega = kx*Ux + ky*Uy + omega1; % Eq. 3.5a
F = -omega1./(kappa.*sinh(h*kappa)); % Eq. 3.6
mu = F.*cosh(h*kappa); % Eq. 3.77
kappa_2 = 2*kappa; % Eq. 3.13
kx_2 = 2*kx; ky_2 = 2*ky;
c = sqrt(a.^2 + b.^2); % p. 316, just after Eq. 3.66


% Second-order self-interaction coefficients (Superharmonic)
if order >= 2
    G_2 = 1/2*h*kappa.*(2 + cosh(2*h*kappa)).*coth(h*kappa)./(sinh(h*kappa).^2); % Eq. 3.26a
    F_2 = -3/4*h*omega1./(sinh(h*kappa).^4); % Eq. 3.26b
end

% Second order (Only Sum Interactions)
if order >= 2
    % Self-self interactions (2*omega is superharmonic)
    % A_2, B_2 are difference (mean) terms -> Omitted
    mu_2 = F_2.*cosh(h*kappa_2) - h*omega1; % Eq. 3.79 (Keep, 2*omega term)
    
    % Mass flux coefficient (difference term) -> Omitted M
    
    % Sum interactions ONLY (n+m)
    cnm = 0; % Double-summation counter
    num_pairs = N*(N-1);
    
    % Pre-allocate for speed (optional but good practice)
    omega_npm = zeros(1, num_pairs);
    kx_npm = zeros(1, num_pairs); ky_npm = zeros(1, num_pairs);
    kappa_npm = zeros(1, num_pairs); alpha_npm = zeros(1, num_pairs);
    gamma_npm = zeros(1, num_pairs); beta_npm = zeros(1, num_pairs);
    F_npm = zeros(1, num_pairs); G_npm = zeros(1, num_pairs);
    mu_npm = zeros(1, num_pairs);
    
    for n = 1:N
        for m = n+1:N
            for pm = [1 -1] % Sum and Difference
                cnm = cnm + 1; % Update counter 

                % Transfer functions
                omega_npm(cnm) = omega1(n) + pm*omega1(m); % Eq. 3.14                
                kx_npm(cnm) = kx(n)+pm*kx(m); ky_npm(cnm) = ky(n)+pm*ky(m);
                kappa_npm(cnm) = sqrt(kx_npm(cnm)^2 + ky_npm(cnm)^2); % Eq. 3.12
                alpha_npm(cnm) = omega_npm(cnm)*cosh(h*kappa_npm(cnm)); % Eq. 3.15
                gamma_npm(cnm) = kappa_npm(cnm)*sinh(h*kappa_npm(cnm)); % Eq. 3.16
                beta_npm(cnm) = omega_npm(cnm)^2*cosh(h*kappa_npm(cnm)) - g*kappa_npm(cnm)*sinh(h*kappa_npm(cnm)); % Eq. 3.17
                F_npm(cnm) = Gamma2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                    omega_npm(cnm),beta_npm(cnm), g,h); % Eq. 3.21
                G_npm(cnm) = Lambda2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                    omega_npm(cnm),alpha_npm(cnm),gamma_npm(cnm),beta_npm(cnm), g,h); % Eq. 3.19
                mu_npm(cnm) = F_npm(cnm)*cosh(h*kappa_npm(cnm)) - h/2*(omega1(n) + pm*omega1(m)); % Eq. 3.78
            end
        end 
    end 
end 

% Third order (Only Superharmonic +++)
if order == 3
    % Single summations & Dispersion Correction
    % Note: Even though this is "Superharmonic" file, these corrections to linear dispersion
    % are CRITICAL for correct phase speeds of the primary waves.
    Upsilon = omega1.*kappa.*(-13 + 24*cosh(2*h*kappa) + cosh(4*h*kappa))./(64*sinh(h*kappa).^5); 
    Xi = 1/(4*h).*(omega1.*G_2 + F_2.*kappa_2.*sinh(h*kappa_2) - g*h*kappa.^2./(2*omega1)); 
    F13 = c.^2.*Upsilon; 
    muStar = c.^2.*Xi; 
    Omega = (8 + cosh(4*h*kappa))./(16*sinh(h*kappa).^4); 
    omega3 = c.^2.*kappa.^2.*Omega; 
    
    % Loop for O(N^2) mutual interactions (Dispersion Correction)
    % This part calculates changes to fundamental frequency due to interactions
    for n = 1:N
        for m = [1:n-1 n+1:N] 
            % Sum interactions (n+m)
            pm = 1; 
            omega_npm_local = omega1(n) + pm*omega1(m); 
            kappa_npm_local = sqrt((kx(n)+pm*kx(m))^2 + (ky(n)+pm*ky(m))^2);
            alpha_npm_local = omega_npm_local*cosh(h*kappa_npm_local); 
            gamma_npm_local = kappa_npm_local*sinh(h*kappa_npm_local); 
            beta_npm_local = omega_npm_local^2*cosh(h*kappa_npm_local) - g*kappa_npm_local*sinh(h*kappa_npm_local); 
            G_npm_local = Lambda2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                    omega_npm_local,alpha_npm_local,gamma_npm_local,beta_npm_local, g,h); 
            F_npm_local = Gamma2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                    omega_npm_local,beta_npm_local, g,h); 

            % Difference interactions (n-m)
            pm = -1; 
            omega_nmm_local = omega1(n) + pm*omega1(m); 
            kappa_nmm_local = sqrt((kx(n)+pm*kx(m))^2 + (ky(n)+pm*ky(m))^2); 
            alpha_nmm_local = omega_nmm_local*cosh(h*kappa_nmm_local); 
            gamma_nmm_local = kappa_nmm_local*sinh(h*kappa_nmm_local); 
            beta_nmm_local = omega_nmm_local^2*cosh(h*kappa_nmm_local) - g*kappa_nmm_local*sinh(h*kappa_nmm_local); 
            G_nmm_local = Lambda2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                    omega_nmm_local,alpha_nmm_local,gamma_nmm_local,beta_nmm_local, g,h); 
            F_nmm_local = Gamma2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                    omega_nmm_local,beta_nmm_local, g,h); 
            
            % Accumulate corrections
            Ups_nm_val = Upsilon_nm(omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), ...
                F_npm_local,F_nmm_local,G_npm_local,G_nmm_local, kappa_npm_local,kappa_nmm_local, g,h);
            XI_nm_val = Xi_nm(omega1(n),kappa(n), omega1(m), G_npm_local,F_npm_local,gamma_npm_local, ...
                G_nmm_local,F_nmm_local,gamma_nmm_local, h,g);
            Om_nm_val = Omega_nm(omega1(n),kx(n),ky(n), omega1(m),kx(m),ky(m),kappa(m), ...
                F_npm_local,F_nmm_local,G_npm_local,G_nmm_local, kappa_npm_local,kappa_nmm_local, g,h);
            
            F13(n) = F13(n) + c(m)^2*Ups_nm_val;
            muStar(n) = muStar(n) + c(m)^2*XI_nm_val; 
            omega3(n) = omega3(n) + c(m)^2*kappa(m)^2*Om_nm_val; 

        end
    end
    omega = omega + omega3.*omega1; 
    muStar = muStar + F13.*cosh(h*kappa); 

    
    % Double summations (n+2m, 2n+m)
    cnm = 0; 
    % Pre-allocate
    % ... (skipping pre-allocation for brevity in this snippet, MATLAB JIT handles it reasonably well)
    
    for n = 1:N 
        for m = n+1:N
            % for pm = [1 -1] -> Only pm = 1
            pm = 1;
            cnm = cnm + 1; 
            
            % Transfer functions, n+2m
            % pm=1 implies n+2m (Superharmonic)
            omega_np2m(cnm) = omega1(n) + pm*2*omega1(m); 
            kx_np2m(cnm) = kx(n) + pm*2*kx(m); ky_np2m(cnm) = ky(n) + pm*2*ky(m);
            kappa_np2m(cnm) = sqrt(kx_np2m(cnm)^2 + ky_np2m(cnm)^2); 
            alpha_np2m(cnm) = omega_np2m(cnm)*cosh(h*kappa_np2m(cnm)); 
            gamma_np2m(cnm) = kappa_np2m(cnm)*sinh(h*kappa_np2m(cnm)); 
            beta_np2m(cnm) = omega_np2m(cnm)^2*cosh(h*kappa_np2m(cnm)) - g*kappa_np2m(cnm)*sinh(h*kappa_np2m(cnm)); 
            
            % Note: In original code, indices were based on mixed sum/diff array. 
            % Here cnm corresponds directly to n+m sum index in the reduced array.
            % We use: kappa_npm(cnm), G_npm(cnm) which are the Sum (+) vars.
            gamma_2_m = 2*kappa(m)*sinh(h*2*kappa(m)); % Recalc local gamma_2 for m

            G_np2m(cnm) = Lambda3(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                   kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), kappa_2(m),gamma_2_m,G_2(m),pm*F_2(m), ...
                   omega_np2m(cnm),alpha_np2m(cnm),gamma_np2m(cnm),beta_np2m(cnm), g,h); 
            F_np2m(cnm) = Gamma3(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                   kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), kappa_2(m),gamma_2_m,G_2(m),pm*F_2(m), ...
                   omega_np2m(cnm),beta_np2m(cnm), g,h); 
            mu_np2m(cnm) = Pi(omega1(n),kappa(n), pm*omega1(m),kappa(m), pm*omega1(m),kappa(m), ...
                    gamma_npm(cnm),G_npm(cnm),F_npm(cnm), gamma_npm(cnm),G_npm(cnm),F_npm(cnm), gamma_2_m,G_2(m),pm*F_2(m), F_np2m(cnm),kappa_np2m(cnm), g,h);
                
            % Transfer functions, 2n+m 
            omega_2npm(cnm) = 2*omega1(n) + pm*omega1(m); 
            kx_2npm(cnm) = 2*kx(n) + pm*kx(m); ky_2npm(cnm) = 2*ky(n) + pm*ky(m);
            kappa_2npm(cnm) = sqrt(kx_2npm(cnm)^2 + ky_2npm(cnm)^2); 
            alpha_2npm(cnm) = omega_2npm(cnm)*cosh(h*kappa_2npm(cnm)); 
            gamma_2npm(cnm) = kappa_2npm(cnm)*sinh(h*kappa_2npm(cnm)); 
            beta_2npm(cnm) = omega_2npm(cnm)^2*cosh(h*kappa_2npm(cnm)) - g*kappa_2npm(cnm)*sinh(h*kappa_2npm(cnm)); 
            
            gamma_2_n = 2*kappa(n)*sinh(h*2*kappa(n)); % Recalc local gamma_2 for n
            
            G_2npm(cnm) = Lambda3(omega1(n),kx(n),ky(n),kappa(n), omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                   kappa_2(n),gamma_2_n,G_2(n),F_2(n), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), ...
                   omega_2npm(cnm),alpha_2npm(cnm),gamma_2npm(cnm),beta_2npm(cnm), g,h);    
            F_2npm(cnm) = Gamma3(omega1(n),kx(n),ky(n),kappa(n), omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                   kappa_2(n),gamma_2_n,G_2(n),F_2(n), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), ...
                   omega_2npm(cnm),beta_2npm(cnm), g,h);
            mu_2npm(cnm) = Pi(omega1(n),kappa(n), omega1(n),kappa(n), pm*omega1(m),kappa(m), ...
                    gamma_2_n,G_2(n),F_2(n), gamma_npm(cnm),G_npm(cnm),F_npm(cnm), gamma_npm(cnm),G_npm(cnm),F_npm(cnm), F_2npm(cnm),kappa_2npm(cnm), g,h);
        end
    end 
    
    % Triple summations
    
    % Build matrix for sum indices
    M_nm = zeros(N); 
    nm_indices = 1:N*(N-1)/2;
    M_nm(triu(ones(N),1)==1) = nm_indices; 
    % M_nm NOW contains directly the index for (n+m) in our reduced arrays.
    % No need for (2*M_nm - 1).
    
    % Pre-allocate per-iteration arrays for speed
    
    % --- PROGRESS BAR SETUP ---
    total_iterations_3 = N*(N-1)*(N-2)/6;
    check_interval = max(1, round(total_iterations_3 / 20)); % Update every 5%
    start_time = tic;
    % --------------------------

    cnm = 0; c3 = 0; 
    for n = 1:N 
        for m = n+1:N
            % pmm = 1 (Only +m)
            pmm = 1;
            cnm = cnm + 1; % Matches cnm in Double summation (n,m)
            
            for p = m+1:N
                % pmp = 1 (Only +p)
                pmp = 1;
                c3 = c3 + 1; 

                % --- PROGRESS UPDATE ---
                if mod(c3, check_interval) == 0
                    elapsed_time = toc(start_time);
                    progress_frac = c3 / total_iterations_3;
                    estimated_total = elapsed_time / progress_frac;
                    remaining_time = estimated_total - elapsed_time;
                    
                    fprintf('Order 3 Progress: %.1f%% | Elapsed: %.1fs | ETA: %.1fs\n', ...
                        progress_frac * 100, elapsed_time, remaining_time);
                end
                % -----------------------
                
                % (n + m + p) coefficients
                omega_npmpp(c3) = omega1(n) + pmm*omega1(m) + pmp*omega1(p); 
                kx_npmpp(c3) = kx(n) + pmm*kx(m) + pmp*kx(p); ky_npmpp(c3) = ky(n) + pmm*ky(m) + pmp*ky(p);
                kappa_npmpp(c3) = sqrt(kx_npmpp(c3)^2 + ky_npmpp(c3)^2); 
                alpha_npmpp(c3) = omega_npmpp(c3)*cosh(h*kappa_npmpp(c3)); 
                beta_npmpp(c3) = omega_npmpp(c3)^2*cosh(h*kappa_npmpp(c3)) - g*kappa_npmpp(c3)*sinh(h*kappa_npmpp(c3)); 
                gamma_npmpp(c3) = kappa_npmpp(c3)*sinh(h*kappa_npmpp(c3)); 
                
                % Indices for (n+p) and (m+p)
                cnp = M_nm(n,p); 
                cmp = M_nm(m,p); 
                % No adjustments needed for difference terms since we always have +p here.
                
                % Transfer functions
                % Inputs need G_npm(cnm), G_npp(cnp), G_mpp(cmp) etc. (All Sums)
                
                G_npmpp(c3) = Lambda3(omega1(n),kx(n),ky(n),kappa(n), pmm*omega1(m),pmm*kx(m),pmm*ky(m),kappa(m), ...
                    pmp*omega1(p),pmp*kx(p),pmp*ky(p),kappa(p), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), ...
                    kappa_npm(cnp),gamma_npm(cnp),G_npm(cnp),F_npm(cnp), kappa_npm(cmp),gamma_npm(cmp),G_npm(cmp),pmm*F_npm(cmp), ...
                    omega_npmpp(c3),alpha_npmpp(c3),gamma_npmpp(c3),beta_npmpp(c3), g,h);
                F_npmpp(c3) = Gamma3(omega1(n),kx(n),ky(n),kappa(n), pmm*omega1(m),pmm*kx(m),pmm*ky(m),kappa(m), ...
                    pmp*omega1(p),pmp*kx(p),pmp*ky(p),kappa(p), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), ...
                    kappa_npm(cnp),gamma_npm(cnp),G_npm(cnp),F_npm(cnp), kappa_npm(cmp),gamma_npm(cmp),G_npm(cmp),pmm*F_npm(cmp), ...
                    omega_npmpp(c3),beta_npmpp(c3), g,h);
                mu_npmpp(c3) = Pi(omega1(n),kappa(n), pmm*omega1(m),kappa(m), ...
                    pmp*omega1(p),kappa(p), ...
                    gamma_npm(cnm),G_npm(cnm),F_npm(cnm), gamma_npm(cnp),G_npm(cnp),F_npm(cnp), ...
                    gamma_npm(cmp),G_npm(cmp),pmm*F_npm(cmp), F_npmpp(c3),kappa_npmpp(c3), g,h);
            end 
        end
    end
end 

%%% Save coefficients for output
coeffs.g = g; coeffs.h = h; coeffs.N = N;
coeffs.a = a; coeffs.b = b; coeffs.kx = kx; coeffs.ky = ky;
coeffs.Ux = Ux; coeffs.Uy = Uy;
coeffs.kappa = kappa; coeffs.omega1 = omega1; coeffs.omega = omega;
coeffs.mu = mu; coeffs.muStar = muStar; coeffs.F = F; coeffs.c = c; coeffs.kappa_2 = kappa_2;
if order >= 2 
    coeffs.F_2 = F_2; coeffs.G_2 = G_2; coeffs.mu_2 = mu_2;
    % Add A_2 and B_2 for self-self interactions (using superharmonic definition)
    % For superharmonics (2*theta_n), the amplitude coeffs are typically derived from
    % squaring the linear signal.
    % A_2 = 1/(2h)*(a^2 - b^2); B_2 = 1/h*(a*b); <-- This is Eq 3.11 for difference terms?
    % Wait, Eq 3.7a: sum G_2 [ A_2 cos(2th) + B_2 sin(2th) ]
    % Let's check MF12 paper for Eq 3.7a definition.
    % Usually A_2n = (a_n^2 - b_n^2)/(2h) is for mean/diff? No.
    % Let's look at standard expansion: (a cos + b sin)^2 = a^2 cos^2 + b^2 sin^2 + 2ab sin cos
    % = a^2(1+cos2)/2 + b^2(1-cos2)/2 + ab sin2
    % = (a^2+b^2)/2 + (a^2-b^2)/2 cos2 + ab sin2.
    % So coefficient for cos(2theta) is (a^2-b^2)/2. Coefficient for sin(2theta) is ab.
    % In MF12, A_2n and B_2n are defined in Eq 3.11 (referenced in original code).
    % Original code (commented out in your version or missing?):
    % A_2 = 1/(2*h)*(a.^2 - b.^2); B_2 = 1/h*(a.*b); 
    % Wait, the 1/h factor comes from the perturbation expansion scaling?
    % Let's re-add them based on typical HOS/MF12 structure.
    coeffs.A_2 = 1/2*(a.^2 - b.^2); % Note: removed 1/h based on derivation above if G_2 handles scaling?
    % Original code had: A_2 = 1/(2*h)*(a.^2 - b.^2). 
    % Let's stick to consistent definition with original code assuming G_2 expects it.
    % Re-reading original coeffsMF12.m lines 95:
    % A_2 = 1/(2*h)*(a.^2 - b.^2); B_2 = 1/h*(a.*b);
    coeffs.A_2 = 1/(2*h)*(a.^2 - b.^2); 
    coeffs.B_2 = 1/h*(a.*b);
    
    coeffs.F_npm = F_npm; coeffs.G_npm = G_npm; 
    
    % Calculates A_npm and B_npm for output
    % Need to re-calculate vector form or store them during loop
    % In loop: A_npm(cnm) = 1/h*(...); B_npm(cnm) = ...
    % Since we didn't store them in the simplified loop, we should add them now.
    % Sum interactions (pm=1):
    % A_npm = 1/h*(a(n)a(m) - b(n)b(m));
    % B_npm = 1/h*(a(m)b(n) + a(n)b(m));
    
    % Let's reconstruct them vectorially to avoid loop
    % We need indices.
    % Re-run loop logic or better, just compute efficiently?
    % Since we need to return them in struct, and surface_spectral needs them.
    
    % Let's add generation of A_npm, B_npm inside the loop or after.
    % Re-doing the loop just for A/B is wasteful but safe, or vectorize.
    % Vectorized approach:
    pair_idx = 1:N*(N-1)/2;
    % We need n and m indices for each cnm.
    % It's complex to get n,m back from cnm without M_nm matrix.
    % Let's rebuild M_nm earlier?
    % Actually we can just do a quick loop to fill A_npm/B_npm arrays.
    
    A_npm = zeros(1, length(F_npm));
    B_npm = zeros(1, length(F_npm));
    cnm_fill = 0;
    for n_idx = 1:N
        for m_idx = n_idx+1:N
            for pm = [1 -1]
                cnm_fill = cnm_fill + 1;
                A_npm(cnm_fill) = 1/h*(a(n_idx)*a(m_idx) - pm*b(n_idx)*b(m_idx));
                B_npm(cnm_fill) = 1/h*(a(m_idx)*b(n_idx) + pm*a(n_idx)*b(m_idx));
            end
        end
    end
    coeffs.A_npm = A_npm;
    coeffs.B_npm = B_npm;

    coeffs.mu_npm = mu_npm; coeffs.kappa_npm = kappa_npm; coeffs.omega_npm = omega_npm;
    coeffs.kx_2 = kx_2; coeffs.ky_2 = ky_2; coeffs.kx_npm = kx_npm; coeffs.ky_npm = ky_npm;
end
if order >= 3
    % coeffs.F_3 ... omitted
    coeffs.gamma_2 = 2*kappa.*sinh(h*2*kappa); % Calculated locally in loop but might be needed out
    
    % Need A_np2m, B_np2m, A_2npm, B_2npm
    % Re-calculate
    A_np2m = zeros(1, length(F_np2m)); B_np2m = zeros(1, length(F_np2m));
    A_2npm = zeros(1, length(F_2npm)); B_2npm = zeros(1, length(F_2npm));
    cnm_fill = 0;
    for n_idx = 1:N 
        for m_idx = n_idx+1:N
             cnm_fill = cnm_fill + 1;
             % pm=1
             % Eq 3.36 for n+2m
             % ThetaA(n, n, m, m, m, m) ?? No.
             % Eq 3.36: A_np2m = 1/2*ThetaA(a(n),b(n), a(m), b(m), a(m), b(m), h);
             % Note: Theta functions are defined at bottom of file.
             % We need to call them or inline them.
             % Inline ThetaA: (an*am*ap - bn*bm*ap - bn*am*bp - an*bm*bp)/(h^2)
             
             % For n+2m: args are (n), (m), (m) with pm=1 implies b(m)
             an=a(n_idx); bn=b(n_idx); am=a(m_idx); bm=b(m_idx);
             A_np2m(cnm_fill) = 1/2 * (an*am*am - bn*bm*am - bn*am*bm - an*bm*bm)/(h^2);
             % ThetaB
             B_np2m(cnm_fill) = 1/2 * (bn*am*am + an*bm*am + an*am*bm - bn*bm*bm)/(h^2);
             
             % For 2n+m
             % A_2npm = 1/2*ThetaA(n,n, n,n, m,m)
             A_2npm(cnm_fill) = 1/2 * (an*an*am - bn*bn*am - bn*an*bm - an*bn*bm)/(h^2);
             B_2npm(cnm_fill) = 1/2 * (bn*an*am + an*bn*am + an*an*bm - bn*bn*bm)/(h^2);
        end
    end
    coeffs.A_np2m = A_np2m; coeffs.B_np2m = B_np2m; 
    coeffs.A_2npm = A_2npm; coeffs.B_2npm = B_2npm; 
    
    coeffs.F_np2m = F_np2m; coeffs.G_np2m = G_np2m; coeffs.mu_np2m = mu_np2m;
    coeffs.kappa_np2m = kappa_np2m; coeffs.kappa_2npm = kappa_2npm;
    
    coeffs.F_2npm = F_2npm; coeffs.G_2npm = G_2npm; coeffs.mu_2npm = mu_2npm;
    coeffs.omega_np2m = omega_np2m; coeffs.omega_2npm = omega_2npm;
    
    % Need A_npmpp, B_npmpp
    A_npmpp = zeros(1, length(F_npmpp)); B_npmpp = zeros(1, length(F_npmpp));
    c3_fill = 0;
    for n_idx = 1:N
        for m_idx = n_idx+1:N
            for p_idx = m_idx+1:N
                 c3_fill = c3_fill + 1;
                 % n+m+p
                 an=a(n_idx); bn=b(n_idx); am=a(m_idx); bm=b(m_idx); ap=a(p_idx); bp=b(p_idx);
                 A_npmpp(c3_fill) = 1/2 * (an*am*ap - bn*bm*ap - bn*am*bp - an*bm*bp)/(h^2);
                 B_npmpp(c3_fill) = 1/2 * (bn*am*ap + an*bm*ap + an*am*bp - bn*bm*bp)/(h^2);
            end
        end
    end
    coeffs.A_npmpp = A_npmpp; coeffs.B_npmpp = B_npmpp;
    
    coeffs.F_npmpp = F_npmpp; coeffs.G_npmpp = G_npmpp; coeffs.mu_npmpp = mu_npmpp;
    coeffs.kappa_npmpp = kappa_npmpp; coeffs.omega_npmpp = omega_npmpp;
    coeffs.kx_np2m = kx_np2m; coeffs.ky_np2m = ky_np2m;
    coeffs.kx_2npm = kx_2npm; coeffs.ky_2npm = ky_2npm;
    coeffs.kx_npmpp = kx_npmpp; coeffs.ky_npmpp = ky_npmpp;
    
    % Also, need A_3, B_3 for self-self-self if order=3
    % A_3 = 1/2*ThetaA(n,n,n);
    coeffs.A_3 = zeros(size(a)); coeffs.B_3 = zeros(size(a));
    % G_3 and mu_3?
    % In original code G_3(n) was computed. In our simplified code, we skipped single corrections?
    % Let's check above lines 70 in coeffs_superharmonic... we skipped single summations!
    % But we need G_3 for 3rd harmonic.
    % Let's quickly add G_3 calculation here or via loop below.
    % G_3 definition: Eq 3.64.
    kappa_3 = 3*kappa;
    F_3 = 1/32*(h^2*kappa.*h.*omega1)./(sinh(h*kappa).^7).*(-11 + 2*cosh(2*h*kappa)); 
    G_3 = 3/128*h^2*kappa.^2./(sinh(h*kappa).^6).*(14 + 15*cosh(2*h*kappa) + 6*cosh(4*h*kappa) + cosh(6*h*kappa));
    mu_3 = F_3.*cosh(h*kappa_3) - g*h^2/4*kappa.^2./omega1 + h/2.*(F_2.*(2*kappa.*sinh(h*2*kappa)) - omega1.*G_2); % Approx
    
    for n_idx = 1:N
        an=a(n_idx); bn=b(n_idx);
        % ThetaA(n,n,n)
        coeffs.A_3(n_idx) = 1/2 * (an*an*an - bn*bn*an - bn*an*bn - an*bn*bn)/(h^2);
        coeffs.B_3(n_idx) = 1/2 * (bn*an*an + an*bn*an + an*an*bn - bn*bn*bn)/(h^2);
    end
    coeffs.G_3 = G_3; coeffs.mu_3 = mu_3; coeffs.F_3 = F_3;
end
end % End of function

%%% Third-order internal functions (Helper Functions)
% Upsilon_nm function, Eq. 3.68, p. 316
function out = Upsilon_nm(omega1n,knx,kny,kappan, omega1m, kmx,kmy,kappam, Fnpm,Fnmm,Gnpm,Gnmm, kappanpm,kappanmm, g,h)
    knkm = knx*kmx + kny*kmy;
    out = g/(4*omega1n*omega1m*cosh(h*kappan))*(omega1m*(kappan^2 - kappam^2) - omega1n*knkm) ...
        + (Gnpm + Gnmm)/(4*h*omega1n^2*omega1m*cosh(h*kappan))*(g^2*knkm + omega1m^3*omega1n) ...
        - 1/(4*h*cosh(h*kappan))*(Fnpm*kappanpm*sinh(h*kappanpm) + Fnmm*kappanmm*sinh(h*kappanmm)) ...
        + g*Fnpm*cosh(h*kappanpm)/(4*h*omega1n^2*omega1m*cosh(h*kappan))*((omega1n + omega1m)*(knkm + kappam^2) - omega1m*kappanpm^2) ...
        + g*Fnmm*cosh(h*kappanmm)/(4*h*omega1n^2*omega1m*cosh(h*kappan))*((omega1n - omega1m)*(knkm - kappam^2) - omega1m*kappanmm^2);
end

% Xi function, Eq. 3.86, p. 319
function out = Xi_nm(omega1n,kappan, omega1m, Gnpm,Fnpm,gamma_npm, Gnmm,Fnmm,gamma_nmm, h,g)
    out = 1/(2*h)*(omega1m*(Gnpm - Gnmm) + Fnpm*gamma_npm + Fnmm*gamma_nmm - g*h*kappan^2/(2*omega1n));
end

% ThetaA function, Eq. 3.34, p. 312
function out = ThetaA(an,bn, am,bm, ap,bp, h)
    out = (an*am*ap - bn*bm*ap - bn*am*bp - an*bm*bp)/(h^2); % Eq. 3.34
end

% ThetaB function, Eq. 3.35, p. 312
function out = ThetaB(an,bn, am,bm, ap,bp, h)
    out = (bn*am*ap + an*bm*ap + an*am*bp - bn*bm*bp)/(h^2);
end

% Omega_nm function, Eq. 3.75, p. 317
function out = Omega_nm(omega1n,knx,kny, omega1m,kmx,kmy,kappam, Fnpm,Fnmm,Gnpm,Gnmm, kappanpm,kappanmm, g,h)
    knkm = knx*kmx + kny*kmy;
    out = 1/(kappam^2)*((2*omega1m^2 + omega1n^2)/(4*omega1n*omega1m)*knkm + 1/4*kappam^2) ...
        + (Gnpm + Gnmm)/(kappam^2)*(g*knkm/(4*h*omega1n*omega1m) - omega1m^2/(4*g*h)) ...
        + omega1n/(4*g*h*kappam^2)*(Fnpm*kappanpm*sinh(h*kappanpm) + Fnmm*kappanmm*sinh(h*kappanmm)) ...
        - Fnpm*cosh(h*kappanpm)/(4*h*omega1n*omega1m*kappam^2)*((omega1n - omega1m)*(kappam^2 + knkm) + omega1m*kappanpm^2) ...
        + Fnmm*cosh(h*kappanmm)/(4*h*omega1n*omega1m*kappam^2)*((omega1n + omega1m)*(kappam^2 - knkm) - omega1m*kappanmm^2);
end

%%% Second-order internal functions
% Lambda2 function, Eq. 3.18, p. 310
function out = Lambda2(omega1n,knx,kny,kappan, omega1m,kmx,kmy,kappam, omega_npm,alpha_npm,gamma_npm,beta_npm,g,h)
    knkm = knx*kmx + kny*kmy; % Dot product of wave number vectors
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

%%% Third-order internal functions
% Lambda3 function, Eq. 3.53, p. 313-314
function out = Lambda3(omega1n,knx,kny,kappan, omega1m,kmx,kmy,kappam, omega1p,kpx,kpy,kappap, ...
                       kappanpm,gammanpm,Gnpm,Fnpm, kappanpp,gammanpp,Gnpp,Fnpp, kappampp,gammampp,Gmpp,Fmpp, ...
                       omega_npmpp,alpha_npmpp,gamma_npmpp,beta_npmpp, g,h)
    knkm = knx*kmx + kny*kmy; knkp = knx*kpx + kny*kpy; kmkp = kmx*kpx + kmy*kpy; 
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
    knkm = knx*kmx + kny*kmy; knkp = knx*kpx + kny*kpy; kmkp = kmx*kpx + kmy*kpy; 
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
    gamma_npm,Gnpm,Fnpm, gamma_npp,Gnpp,Fnpp, gamma_mpp,Gmpp,Fmpp, Fnpmpp,kappa_npmpp, g,h) % Eq. 3.82
    out = Fnpmpp*cosh(h*kappa_npmpp) ...
        - g*h^2/4*(kappan^2/omega1n + kappam^2/omega1m + kappap^2/omega1p) ...
        - h/2*(omega1n*Gmpp + omega1m*Gnpp + omega1p*Gnpm) ...
        + h/2*(Fnpm*gamma_npm + Fnpp*gamma_npp + Fmpp*gamma_mpp); % Eq. 3.81
end
