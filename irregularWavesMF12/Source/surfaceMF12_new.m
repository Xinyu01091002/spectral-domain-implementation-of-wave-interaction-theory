function [eta,phiS,Mx,My] = surfaceMF12(order,coeffs,x,y,t)
% A function to evaluate the third-order multi-directional irregular wave 
% theory of Madsen & Fuhrman (2012, MF12).  The function returns the free 
% surface elevation (eta), the free surface velocity potential (phiS), as 
% well as the components of the resulting mass flux vector (Mx,My). 
% Equation numbers in comments below correspond to those of MF12.
%
% Reference:
%
% Madsen, P.A. & Fuhrman, D.R. (2012) Third-order theory for
% multi-directional irregular waves. J. Fluid Mech. 698, 304-334.
%
% Programmed by David R. Fuhrman, November 17-23, 2022

% Initialize
eta = 0.*t;
phiS = 0*t + coeffs.Ux*x + coeffs.Uy*y; % First part of Eq. 3.76
Mx = coeffs.h*coeffs.Ux; My = coeffs.h*coeffs.Uy; % First part of Eq. 3.70

% First order
for n = 1:coeffs.N
    theta_n = coeffs.omega(n).*t - coeffs.kx(n)*x - coeffs.ky(n)*y; % Eq. 3.1
    eta = eta + coeffs.a(n)*cos(theta_n) + coeffs.b(n)*sin(theta_n); % Eq. 3.2
    phiS = phiS + (coeffs.mu(n) + coeffs.muStar(n))* ...
        (coeffs.a(n)*sin(theta_n) - coeffs.b(n)*cos(theta_n)); % Eq. 3.76 (linear part)
end

% Second order
if order >= 2
    % Self-self interactions
    for n = 1:coeffs.N 
        % Superimpose results
        theta_n = coeffs.omega(n).*t - coeffs.kx(n)*x - coeffs.ky(n)*y; % Eq. 3.1
        eta = eta + coeffs.G_2(n)*(coeffs.A_2(n)*cos(2*theta_n) + coeffs.B_2(n)*sin(2*theta_n)); % Eq. 3.7a
        phiS = phiS + coeffs.mu_2(n)*(coeffs.A_2(n)*sin(2*theta_n) - coeffs.B_2(n)*cos(2*theta_n)); % Eq. 3.76 (second line)
        
        % Mass flux
        Mx = Mx + coeffs.M(n)*coeffs.kx(n); My = My + coeffs.M(n)*coeffs.ky(n); % Eq. 3.70 (vector components)
    end
    
    % Sum and difference interactions
    cnm = 0; % Double-summation counter
    for n = 1:coeffs.N
        for m = n+1:coeffs.N
            for pm = [1 -1] % Loop over both sum & difference (pm = +/-)
                cnm = cnm + 1; % Update counter

                % Superimpose results
                theta_n = coeffs.omega(n).*t - coeffs.kx(n)*x - coeffs.ky(n)*y;
                theta_m = coeffs.omega(m).*t - coeffs.kx(m)*x - coeffs.ky(m)*y;
                theta_npm = theta_n + pm*theta_m; % Eq. 3.9
                eta = eta + coeffs.G_npm(cnm)*(coeffs.A_npm(cnm)*cos(theta_npm) ...
                                             + coeffs.B_npm(cnm)*sin(theta_npm)); % Eq. 3.7 (second line)
                phiS = phiS + coeffs.mu_npm(cnm)*(coeffs.A_npm(cnm)*sin(theta_npm) ...
                                                - coeffs.B_npm(cnm)*cos(theta_npm)); % Eq. 3.76 (third line)
            end % Ends pm loop
        end % Ends m loop
    end % Ends n loop
end % End, if order>=2

% Third order
if order == 3
    % Single summations
    for n = 1:coeffs.N 
        % Superimpose results
        theta_n = coeffs.omega(n).*t - coeffs.kx(n)*x - coeffs.ky(n)*y; % Eq. 3.1
        %eta = eta + 0*coeffs.G13(n)*(coeffs.a(n)*cos(theta_n) + coeffs.b(n)*sin(theta_n)); % Eq. 3.28 (first line), see p. 316 (top), G13n=0 here
        eta = eta + coeffs.G_3(n)*(coeffs.A_3(n)*cos(3*theta_n) + coeffs.B_3(n)*sin(3*theta_n)); % Eq. 3.28 (second line)
        phiS = phiS + coeffs.mu_3(n)*(coeffs.A_3(n)*sin(3*theta_n) - coeffs.B_3(n)*cos(3*theta_n)); % Eq. 3.76 (fourth line)
    end
    
    % Double summations
    cnm = 0; % Double-summation counter
    for n = 1:coeffs.N 
        for m = n+1:coeffs.N
            for pm = [1 -1] % +/- m
                cnm = cnm + 1; % Update counter
                
                % Superimpose results 
                theta_n = coeffs.omega(n).*t - coeffs.kx(n)*x - coeffs.ky(n)*y; % Eq. 3.1
                theta_m = coeffs.omega(m).*t - coeffs.kx(m)*x - coeffs.ky(m)*y;
                theta_np2m = theta_n + pm*2*theta_m; % Eq. 3.30a
                if pm==1
                eta = eta + coeffs.G_np2m(cnm)*(coeffs.A_np2m(cnm)*cos(theta_np2m) ...
                                              + coeffs.B_np2m(cnm)*sin(theta_np2m)); % Eq. 3.28 (third line)
                
                phiS = phiS + coeffs.mu_np2m(cnm)*(coeffs.A_np2m(cnm)*sin(theta_np2m) ...
                                                 - coeffs.B_np2m(cnm)*cos(theta_np2m)); % Eq. 3.76 (fifth line)
                end
                theta_2npm = 2*theta_n + pm*theta_m; % Eq. 3.30b
                if pm==1
                eta = eta + coeffs.G_2npm(cnm)*(coeffs.A_2npm(cnm)*cos(theta_2npm) ...
                                              + coeffs.B_2npm(cnm)*sin(theta_2npm)); % Eq. 3.28 (fourth line)
                
                phiS = phiS + coeffs.mu_2npm(cnm)*(coeffs.A_2npm(cnm)*sin(theta_2npm) ...
                                                 - coeffs.B_2npm(cnm)*cos(theta_2npm)); % Eq. 3.76 (sixth line)
                end
            end
        end
    end % End of double summation
      
    % Triple summations
    c3 = 0; % Initialize counter
    for n = 1:coeffs.N 
        for m = n+1:coeffs.N
            for pmm = [1 -1] % +/- m
                for p = m+1:coeffs.N
                    for pmp = [1 -1] % +/- p
                        c3 = c3 + 1; % Update counter for (n +/- m +/- p) quantities
                        
                        % Superimpose results
                        theta_n = coeffs.omega(n).*t - coeffs.kx(n)*x - coeffs.ky(n)*y; % Eq. 3.1
                        theta_m = coeffs.omega(m).*t - coeffs.kx(m)*x - coeffs.ky(m)*y;
                        theta_p = coeffs.omega(p).*t - coeffs.kx(p)*x - coeffs.ky(p)*y;
                        theta_npmpp = theta_n + pmm*theta_m + pmp*theta_p; % Eq. 3.31
                        if pmm==1 && pmp==1
                        eta = eta + 2*coeffs.G_npmpp(c3)*(coeffs.A_npmpp(c3)*cos(theta_npmpp) ...
                                                      + coeffs.B_npmpp(c3)*sin(theta_npmpp)); % Eq. 3.28 (last line)
                        
                        phiS = phiS + 2*coeffs.mu_npmpp(c3)*(coeffs.A_npmpp(c3)*sin(theta_npmpp) ...
                                                         - coeffs.B_npmpp(c3)*cos(theta_npmpp)); % Eq. 3.76 (last line)
                        end     
                    end % End of pmp loop
                end % End of p loop
            end % End of pmm loop
        end % End of m loop
    end % End of triple summation
end % End of third order

end % End of function
