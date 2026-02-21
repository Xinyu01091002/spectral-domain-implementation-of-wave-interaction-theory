function [coeffs] = coeffsMF12(order,g,h,a,b,kx,ky,Ux,Uy, dispCoeffs)
% A function to evaluate the third-order multi-directional irregular wave 
% theory of Madsen & Fuhrman (2012) (MF12), also utilizing corrections from 
% the appendix of Fuhrman et al. (2023).  Equation numbers in comments 
% below correspond to those of MF12, with corrections also indicated when
% relevant.
%
% References:
%
% Fuhrman et al. (2023) A new probability density function for the surface
% elevation in irregular seas. J. Fluid Mech. 970, A38. 
% https://doi.org/10.1017/jfm.2023.669.
%
% Madsen, P.A. & Fuhrman, D.R. (2012) Third-order theory for
% multi-directional irregular waves. J. Fluid Mech. 698, 304-334.
%
% Programmed by David R. Fuhrman, November, 2022

% Set (optional) default argument dispCoeffs to 0 (false)
if nargin < 10, dispCoeffs = 0; end

% Initialize
muStar = 0*a; 

% Basic derived input
N = length(a); % Number of components
progress_enabled = false;
progress_pct_step = 5;
if strcmp(getenv('MF12_PROGRESS'), '1')
    progress_enabled = N >= 200;
end
kappa = sqrt(kx.^2 + ky.^2); % Eq. 3.4
omega1 = sqrt(g*kappa.*tanh(h*kappa)); % Eq. 3.5b
omega = kx*Ux + ky*Uy + omega1; % Eq. 3.5a
F = -omega1./(kappa.*sinh(h*kappa)); % Eq. 3.6
mu = F.*cosh(h*kappa); % Eq. 3.77
kappa_2 = 2*kappa; % Eq. 3.13
kx_2 = 2*kx; ky_2 = 2*ky;
c = sqrt(a.^2 + b.^2); % p. 316, just after Eq. 3.66

% Third-order corrections to first-order quantities
if order >= 2
    G_2 = 1/2*h*kappa.*(2 + cosh(2*h*kappa)).*coth(h*kappa)./(sinh(h*kappa).^2); % Eq. 3.26a
    F_2 = -3/4*h*omega1./(sinh(h*kappa).^4); % Eq. 3.26b
end
if order == 3
    Upsilon = omega1.*kappa.*(-13 + 24*cosh(2*h*kappa) + cosh(4*h*kappa))./(64*sinh(h*kappa).^5); % Eq. 3.67
    Xi = 1/(4*h).*(omega1.*G_2 + F_2.*kappa_2.*sinh(h*kappa_2) - g*h*kappa.^2./(2*omega1)); % Eq. 3.85
    F13 = c.^2.*Upsilon; % First part of Eq. 3.66
    muStar = c.^2.*Xi; % Second part of Eq. 3.84
    Omega = (8 + cosh(4*h*kappa))./(16*sinh(h*kappa).^4); % Eq. 3.74
    omega3 = c.^2.*kappa.^2.*Omega; % First part of Eq. 3.73
    cnm = 0; % Counter
    if progress_enabled
        total_cnm_o3 = N*(N-1);
        step_cnm_o3 = max(1, floor(total_cnm_o3*progress_pct_step/100));
        t_cnm_o3 = tic;
    end
    for n = 1:N
        for m = [1:n-1 n+1:N] % m = 1:N, m~=n
            cnm = cnm + 1; % Update counter
            if progress_enabled && (mod(cnm, step_cnm_o3) == 0 || cnm == total_cnm_o3)
                elapsed = toc(t_cnm_o3);
                pct = 100*cnm/total_cnm_o3;
                eta = elapsed*(100/pct - 1);
                fprintf('MF12 order3 (n,m): %d/%d (%.1f%%) elapsed %.1fs, ETA %.1fs\n', ...
                    cnm, total_cnm_o3, pct, elapsed, eta);
            end
            
            % Sum interactions
            pm = 1; % n+m
            omega_npm = omega1(n) + pm*omega1(m); % Eq. 3.14
            kappa_npm = sqrt((kx(n)+pm*kx(m))^2 + (ky(n)+pm*ky(m))^2); % Eq. 3.12
            alpha_npm = omega_npm*cosh(h*kappa_npm); % Eq. 3.15
            gamma_npm = kappa_npm*sinh(h*kappa_npm); % Eq. 3.16
            beta_npm = omega_npm^2*cosh(h*kappa_npm) - g*kappa_npm*sinh(h*kappa_npm); % Eq. 3.17
            G_npm = Lambda2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                    omega_npm,alpha_npm,gamma_npm,beta_npm, g,h); % Eq. 3.19
            F_npm = Gamma2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                    omega_npm,beta_npm, g,h); % Eq. 3.21
            
            % Difference interactions
            pm = -1; % n-m
            omega_nmm = omega1(n) + pm*omega1(m); % Eq. 3.14
            kappa_nmm = sqrt((kx(n)+pm*kx(m))^2 + (ky(n)+pm*ky(m))^2); % Eq. 3.12
            alpha_nmm = omega_nmm*cosh(h*kappa_nmm); % Eq. 3.15
            gamma_nmm = kappa_nmm*sinh(h*kappa_nmm); % Eq. 3.16
            beta_nmm = omega_nmm^2*cosh(h*kappa_nmm) - g*kappa_nmm*sinh(h*kappa_nmm); % Eq. 3.17
            G_nmm = Lambda2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                    omega_nmm,alpha_nmm,gamma_nmm,beta_nmm, g,h); % Eq. 3.19
            F_nmm = Gamma2(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                    omega_nmm,beta_nmm, g,h); % Eq. 3.21            
            
            % Transfer functions
            Ups_nm(cnm) = Upsilon_nm(omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), ...
                F_npm,F_nmm,G_npm,G_nmm, kappa_npm,kappa_nmm ,g,h);
            XI_nm(cnm) = Xi_nm(omega1(n),kappa(n), omega1(m), G_npm,F_npm,gamma_npm, G_nmm,F_nmm,gamma_nmm, h,g);
            Om_nm(cnm) = Omega_nm(omega1(n),kx(n),ky(n), omega1(m),kx(m),ky(m),kappa(m), F_npm,F_nmm,G_npm,G_nmm, kappa_npm,kappa_nmm, g,h);
            F13(n) = F13(n) + c(m)^2*Ups_nm(cnm); % Eq. 3.66 (second part)
            muStar(n) = muStar(n) + c(m)^2*XI_nm(cnm); % Eq. 3.84 last part (corrected, see Fuhrman et al., Eq. A3)
            omega3(n) = omega3(n) + c(m)^2*kappa(m)^2*Om_nm(cnm); % Eq. 3.73 (second part)
        end
    end
    omega = omega + omega3.*omega1; % Eq. 3.72 (adds additional contributions to Eq. 3.5a)
    muStar = muStar + F13.*cosh(h*kappa); % First part of Eq. 3.84 (corrected, see Fuhrman et al. 2023, Eq. A3)
end

% Second order,
if order >= 2
    % Self-self interactions
    A_2 = 1/(2*h)*(a.^2 - b.^2); B_2 = 1/h*(a.*b); % Eq. 3.11 (corrected here - See Fuhrman et al. 2023, Eq. A1)
    mu_2 = F_2.*cosh(h*kappa_2) - h*omega1; % Eq. 3.79
    
    % Mass flux coefficient
    M = c.^2.*omega1./(2*kappa).*coth(h*kappa); % Factor from Eq. 3.70
    
    % Sum and difference interactions
    cnm = 0; % Double-summation counter
    if progress_enabled
        total_cnm_o2 = N*(N-1);
        step_cnm_o2 = max(1, floor(total_cnm_o2*progress_pct_step/100));
        t_cnm_o2 = tic;
    end
    for n = 1:N
        for m = n+1:N
            for pm = [1 -1] % Loop over both sum & difference (pm = +/-)
                cnm = cnm + 1; % Update counter
                if progress_enabled && (mod(cnm, step_cnm_o2) == 0 || cnm == total_cnm_o2)
                    elapsed = toc(t_cnm_o2);
                    pct = 100*cnm/total_cnm_o2;
                    eta = elapsed*(100/pct - 1);
                    fprintf('MF12 order2 (n,m,pm): %d/%d (%.1f%%) elapsed %.1fs, ETA %.1fs\n', ...
                        cnm, total_cnm_o2, pct, elapsed, eta);
                end

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
                A_npm(cnm) = 1/h*(a(n)*a(m) - pm*b(n)*b(m)); % Eq. 3.10a
                B_npm(cnm) = 1/h*(a(m)*b(n) + pm*a(n)*b(m)); % Eq. 3.10b                
            end % Ends pm loop
        end % Ends m loop
    end % Ends n loop
end % End, if order>=2

% Third order
if order == 3
    % Single summations
    kappa_3 = 3*kappa; % Eq. 3.43
    gamma_2 = kappa_2.*sinh(h*kappa_2); % From Eq. 3.16    
    for n = 1:N 
        % Transfer functions
        A_3(n) = 1/2*ThetaA(a(n),b(n), a(n),b(n), a(n),b(n), h); % Eq. 3.38
        B_3(n) = 1/2*ThetaB(a(n),b(n), a(n),b(n), a(n),b(n), h); % Eq. 3.39
        F_3(n) = 1/32*(h^2*kappa(n)*h*omega1(n))/(sinh(h*kappa(n))^7)*(-11 + 2*cosh(2*h*kappa(n))); % Eq. 3.65
        G_3(n) = 3/128*h^2*kappa(n)^2/(sinh(h*kappa(n))^6)*(14 + 15*cosh(2*h*kappa(n)) + 6*cosh(4*h*kappa(n)) + cosh(6*h*kappa(n))); % Eq. 3.64
        mu_3(n) = F_3(n)*cosh(h*kappa_3(n)) - g*h^2/4*kappa(n)^2/omega1(n) + h/2*(F_2(n)*gamma_2(n) - omega1(n)*G_2(n)); % Eq. 3.80
    end
    
    % Double summations
    cnm = 0; % Double-summation counter
    if progress_enabled
        total_cnm_o3b = N*(N-1);
        step_cnm_o3b = max(1, floor(total_cnm_o3b*progress_pct_step/100));
        t_cnm_o3b = tic;
    end
    for n = 1:N 
        for m = n+1:N
            for pm = [1 -1] % +/- m
                cnm = cnm + 1; % Update counter
                if progress_enabled && (mod(cnm, step_cnm_o3b) == 0 || cnm == total_cnm_o3b)
                    elapsed = toc(t_cnm_o3b);
                    pct = 100*cnm/total_cnm_o3b;
                    eta = elapsed*(100/pct - 1);
                    fprintf('MF12 order3 (n,m,pm): %d/%d (%.1f%%) elapsed %.1fs, ETA %.1fs\n', ...
                        cnm, total_cnm_o3b, pct, elapsed, eta);
                end
                
                % Transfer functions, n+2m & n-2m
                omega_np2m(cnm) = omega1(n) + pm*2*omega1(m); % Eq. 3.44b
                kx_np2m(cnm) = kx(n) + pm*2*kx(m); ky_np2m(cnm) = ky(n) + pm*2*ky(m);
                kappa_np2m(cnm) = sqrt(kx_np2m(cnm)^2 + ky_np2m(cnm)^2); % Eq. 3.41
                alpha_np2m(cnm) = omega_np2m(cnm)*cosh(h*kappa_np2m(cnm)); % Eq. 3.51a
                gamma_np2m(cnm) = kappa_np2m(cnm)*sinh(h*kappa_np2m(cnm)); % Eq. 3.51b
                beta_np2m(cnm) = omega_np2m(cnm)^2*cosh(h*kappa_np2m(cnm)) - g*kappa_np2m(cnm)*sinh(h*kappa_np2m(cnm)); % Eq. 3.47
                A_np2m(cnm) = 1/2*ThetaA(a(n),b(n), a(m),pm*b(m), a(m),pm*b(m), h); % Eq. 3.36
                B_np2m(cnm) = 1/2*ThetaB(a(n),b(n), a(m),pm*b(m), a(m),pm*b(m), h); % Eq. 3.36
                G_np2m(cnm) = Lambda3(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                       kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), kappa_2(m),gamma_2(m),G_2(m),pm*F_2(m), ...
                       omega_np2m(cnm),alpha_np2m(cnm),gamma_np2m(cnm),beta_np2m(cnm), g,h); 
                F_np2m(cnm) = Gamma3(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                       kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), kappa_2(m),gamma_2(m),G_2(m),pm*F_2(m), ...
                       omega_np2m(cnm),beta_np2m(cnm), g,h); 
                mu_np2m(cnm) = Pi(omega1(n),kappa(n), pm*omega1(m),kappa(m), pm*omega1(m),kappa(m), ...
                        gamma_npm(cnm),G_npm(cnm),F_npm(cnm), gamma_npm(cnm),G_npm(cnm),F_npm(cnm), gamma_2(m),G_2(m),pm*F_2(m), F_np2m(cnm),kappa_np2m(cnm), g,h);
                    % Note the correction to the arguments, see Fuhrman et al. (2023), Eq. B2
                    
                % Transfer functions, 2n+m & 2n-m   
                omega_2npm(cnm) = 2*omega1(n) + pm*omega1(m); % Eq. 3.44c
                kx_2npm(cnm) = 2*kx(n) + pm*kx(m); ky_2npm(cnm) = 2*ky(n) + pm*ky(m);
                kappa_2npm(cnm) = sqrt(kx_2npm(cnm)^2 + ky_2npm(cnm)^2); % Eq. 3.42
                alpha_2npm(cnm) = omega_2npm(cnm)*cosh(h*kappa_2npm(cnm)); % Eq. 3.52a
                gamma_2npm(cnm) = kappa_2npm(cnm)*sinh(h*kappa_2npm(cnm)); % Eq. 3.52b
                A_2npm(cnm) = 1/2*ThetaA(a(n),b(n), a(n),b(n), a(m),pm*b(m), h); % From Eq. 3.32
                B_2npm(cnm) = 1/2*ThetaB(a(n),b(n), a(n),b(n), a(m),pm*b(m), h); % From Eq. 3.33
                beta_2npm(cnm) = omega_2npm(cnm)^2*cosh(h*kappa_2npm(cnm)) - g*kappa_2npm(cnm)*sinh(h*kappa_2npm(cnm)); % Eq. 3.48
                G_2npm(cnm) = Lambda3(omega1(n),kx(n),ky(n),kappa(n), omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                       kappa_2(n),gamma_2(n),G_2(n),F_2(n), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), ...
                       omega_2npm(cnm),alpha_2npm(cnm),gamma_2npm(cnm),beta_2npm(cnm), g,h); % Terms like: G_(2n+m), G_(2n-m), etc.   
                F_2npm(cnm) = Gamma3(omega1(n),kx(n),ky(n),kappa(n), omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                       kappa_2(n),gamma_2(n),G_2(n),F_2(n), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), kappa_npm(cnm),gamma_npm(cnm),G_npm(cnm),F_npm(cnm), ...
                       omega_2npm(cnm),beta_2npm(cnm), g,h);
                mu_2npm(cnm) = Pi(omega1(n),kappa(n), omega1(n),kappa(n), pm*omega1(m),kappa(m), ...
                        gamma_2(n),G_2(n),F_2(n), gamma_npm(cnm),G_npm(cnm),F_npm(cnm), gamma_npm(cnm),G_npm(cnm),F_npm(cnm), F_2npm(cnm),kappa_2npm(cnm), g,h);
            end
        end
    end % End of double summation
    
    % Build upper-triangular matrix of n+m indices (used in triple summation loop below)
    M_nm = zeros(N); % Initialize
    nm_indices = 1:N*(N-1)/2; % Vector of indices
    M_nm(triu(ones(N),1)==1) = nm_indices; % Matrix of indices (needs to be moved up)
    M_nm = 2*M_nm - 1; % Index for n+m (factor 2 is because a single vectors is used to store both n+m and n-m components)
    
    % Triple summations
    cnm = 0; c3 = 0; % Initialize counters
    if progress_enabled
        total_c3 = 4 * (N*(N-1)*(N-2) / 6);
        step_c3 = max(1, floor(total_c3*progress_pct_step/100));
        t_c3 = tic;
    end
    for n = 1:N 
        for m = n+1:N
            for pmm = [1 -1] % +/- m
                cnm = cnm + 1; % Counter for (m +/- n) quantities
                for p = m+1:N
                    for pmp = [1 -1] % +/- p
                        c3 = c3 + 1; % Update counter for (n +/- m +/- p) quantities
                        if progress_enabled && (mod(c3, step_c3) == 0 || c3 == total_c3)
                            elapsed = toc(t_c3);
                            pct = 100*c3/total_c3;
                            eta = elapsed*(100/pct - 1);
                            fprintf('MF12 order3 (n,m,p,pm): %d/%d (%.1f%%) elapsed %.1fs, ETA %.1fs\n', ...
                                c3, total_c3, pct, elapsed, eta);
                        end
                        
                        % (n +/- m +/- p) coefficients
                        omega_npmpp(c3) = omega1(n) + pmm*omega1(m) + pmp*omega1(p); % Eq. 3.44a
                        kx_npmpp(c3) = kx(n) + pmm*kx(m) + pmp*kx(p); ky_npmpp(c3) = ky(n) + pmm*ky(m) + pmp*ky(p);
                        kappa_npmpp(c3) = sqrt(kx_npmpp(c3)^2 + ky_npmpp(c3)^2); % Eq. 3.40
                        alpha_npmpp(c3) = omega_npmpp(c3)*cosh(h*kappa_npmpp(c3)); % Eq. 3.50a
                        beta_npmpp(c3) = omega_npmpp(c3)^2*cosh(h*kappa_npmpp(c3)) - g*kappa_npmpp(c3)*sinh(h*kappa_npmpp(c3)); % Eq. 3.46
                        gamma_npmpp(c3) = kappa_npmpp(c3)*sinh(h*kappa_npmpp(c3)); % Eq. 3.50b
                        A_npmpp(c3) = 1/2*ThetaA(a(n),b(n), a(m),pmm*b(m), a(p),pmp*b(p), h); % Eq. 3.32
                        B_npmpp(c3) = 1/2*ThetaB(a(n),b(n), a(m),pmm*b(m), a(p),pmp*b(p), h); % Eq. 3.33
                        
                        % Find cnp (n +/- p) and cmp (+/- m +/- p) indices for function input below
                        cnp = M_nm(n,p); cmp = M_nm(m,p); % Indices for (n+p) and (m+p)
                        if pmp == -1, cnp = cnp + 1; end % Adjust cnp (add one if n-p)
                        if pmm*pmp == -1, cmp = cmp + 1; end % Adjust cmp (add one if -m+p or m-p; If m+p or -m-p then no adjustment is needed!)
                        
                        % Transfer functions
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
                    end % End of pmp loop
                end % End of p loop
            end % End of pmm loop
        end % End of m loop
    end % End of triple summation
end % End of third order

%%% Save coefficients for output
% Linear coefficients
coeffs.g = g; coeffs.h = h; coeffs.N = N;
coeffs.a = a; coeffs.b = b; coeffs.kx = kx; coeffs.ky = ky;
coeffs.Ux = Ux; coeffs.Uy = Uy;
coeffs.kappa = kappa; coeffs.omega1 = omega1; coeffs.omega = omega;
coeffs.mu = mu; coeffs.muStar = muStar; coeffs.F = F; coeffs.c = c; coeffs.kappa_2 = kappa_2;
if order >= 2 % Second-order coefficients
    coeffs.A_2 = A_2; coeffs.B_2 = B_2;
    coeffs.F_2 = F_2; coeffs.G_2 = G_2; coeffs.mu_2 = mu_2;
    coeffs.F_npm = F_npm; coeffs.G_npm = G_npm; coeffs.A_npm = A_npm; coeffs.B_npm = B_npm;
    coeffs.mu_npm = mu_npm; coeffs.kappa_npm = kappa_npm; coeffs.omega_npm = omega_npm;
    coeffs.kx_2 = kx_2; coeffs.ky_2 = ky_2; coeffs.kx_npm = kx_npm; coeffs.ky_npm = ky_npm;
    coeffs.M = M;
end
if order >= 3 % Third-order coefficients
    coeffs.A_3 = A_3; coeffs.B_3 = B_3; coeffs.F_3 = F_3; coeffs.G_3 = G_3; 
    coeffs.mu_3 = mu_3; coeffs.kappa_3 = kappa_3;
    coeffs.F13 = F13;
    coeffs.A_np2m = A_np2m; coeffs.B_np2m = B_np2m; 
    coeffs.F_np2m = F_np2m; coeffs.G_np2m = G_np2m; coeffs.mu_np2m = mu_np2m;
    coeffs.kappa_np2m = kappa_np2m; coeffs.kappa_2npm = kappa_2npm;
    coeffs.A_2npm = A_2npm; coeffs.B_2npm = B_2npm; 
    coeffs.F_2npm = F_2npm; coeffs.G_2npm = G_2npm; coeffs.mu_2npm = mu_2npm;
    coeffs.omega_np2m = omega_np2m; coeffs.omega_2npm = omega_2npm;
    coeffs.A_npmpp = A_npmpp; coeffs.B_npmpp = B_npmpp;
    coeffs.F_npmpp = F_npmpp; coeffs.G_npmpp = G_npmpp; coeffs.mu_npmpp = mu_npmpp;
    coeffs.kappa_npmpp = kappa_npmpp; coeffs.omega_npmpp = omega_npmpp;
    coeffs.kx_np2m = kx_np2m; coeffs.ky_np2m = ky_np2m;
    coeffs.kx_2npm = kx_2npm; coeffs.ky_2npm = ky_2npm;
    coeffs.kx_npmpp = kx_npmpp; coeffs.ky_npmpp = ky_npmpp;
end

% Display higher-order output coefficients
if dispCoeffs == 1
    disp(' '); disp('Transfer function coefficients:')
    if order >= 2
        disp(['G_2n = ' num2str(G_2)]), disp(['F_2n = ' num2str(F_2)]), disp(['mu_2n = ' num2str(mu_2)])
        disp(['G_npm = ' num2str(G_npm)]), disp(['F_npm = ' num2str(F_npm)]), disp(['mu_npm = ' num2str(mu_npm)])        
    end
    if order == 3
        disp(['G_3n = ' num2str(G_3)]), disp(['F_3n = ' num2str(F_3)]), disp(['mu_3n = ' num2str(mu_3)])
        disp(['G_npmpp = ' num2str(G_npmpp)]), disp(['F_npmpp = ' num2str(F_npmpp)]), disp(['mu_npmpp = ' num2str(mu_npmpp)])
        disp(['G_np2m = ' num2str(G_np2m)]), disp(['F_np2m = ' num2str(F_np2m)]), disp(['mu_np2m = ' num2str(mu_np2m)])
        disp(['G_2npm = ' num2str(G_2npm)]), disp(['F_2npm = ' num2str(F_2npm)]), disp(['mu_2npm = ' num2str(mu_2npm)])
    end
    if order == 3
        disp(' '), disp('Other important coefficients:')
        disp(['Upsilon_nn = ' num2str(Upsilon)]), disp(['Upsilon_nm = ' num2str(Ups_nm)])
        disp(['Xi_nn = ' num2str(Xi)]), disp(['Xi_nm = ' num2str(XI_nm)])
        disp(['Omega_nn = ' num2str(Omega)]), disp(['Omega_nm = ' num2str(Om_nm)])
    end
end

end % End of function


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