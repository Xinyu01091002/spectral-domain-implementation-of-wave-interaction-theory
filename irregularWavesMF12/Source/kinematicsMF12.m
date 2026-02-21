function [u,v,w,p,phi, uV,vV,a_x,a_y] = kinematicsMF12(order,coeffs,x,y,z,t)
% A function to evaluate velocity kinematics based on the third-order 
% multi-directional irregular wave theory of Madsen & Fuhrman (2012) (MF12).  
% The function returns distributions for the (x,y,z) velocity components 
% (u,v,w), the dynamic pressure p=p^+/rho, and the velocity potential (phi).  
% Additionally, the function returns uV=u*V and vV=v*V, where V=sqrt(u^2+v^2), 
% and the total acceleration a_x=du/dt+u*du/dx+v*du/dy+w*du/dz (similarly 
% for a_y) e.g. for potential use in force calculations using the Morison 
% equation. Equation numbers in comments below correspond to those of MF12.
%
% Reference:
%
% Madsen, P.A. & Fuhrman, D.R. (2012) Third-order theory for
% multi-directional irregular waves. J. Fluid Mech. 698, 304-334.
%
% Programmed by David R. Fuhrman, November 24, 2022

% Initialize
phi = coeffs.Ux*x + coeffs.Uy*y; % Velocity potential, first part of Eq. 3.76
Z = z + coeffs.h; % p. 308
u = coeffs.Ux; v = coeffs.Uy; w = 0; p = 0; % u, v, w, p=p^+/rho
V = sqrt(u^2 + v^2)*u; uV = u*V; vV = v*V; % (u,v)*V profiles for computing drag forces
a_x = 0; a_y = 0; % Acceleration profiles for computing inertial forces

% First order
for n = 1:coeffs.N
    % Determine common coefficients and superimpose results
    omega = coeffs.omega(n); kappa = coeffs.kappa(n);
    kx = coeffs.kx(n); ky = coeffs.ky(n);
    F = coeffs.F(n); A = coeffs.a(n); B = coeffs.b(n);
    [phi,u,v,w,p, uV,vV,a_x,a_y] = addKinematics(phi,u,v,w,p, uV,vV,a_x,a_y, omega,kx,ky,kappa,Z,F,A,B, x,y,t);
end

% Second order
if order >= 2
    % Self-self interactions
    for n = 1:coeffs.N 
        omega = 2*coeffs.omega(n); kappa = coeffs.kappa_2(n); % !!!
        kx = coeffs.kx_2(n); ky = coeffs.ky_2(n);         
        F = coeffs.F_2(n); A = coeffs.A_2(n); B = coeffs.B_2(n);
        [phi,u,v,w,p, uV,vV,a_x,a_y] = addKinematics(phi,u,v,w,p, uV,vV,a_x,a_y, omega,kx,ky,kappa,Z,F,A,B, x,y,t);
    end
    
    % Sum and difference interactions
    cnm = 0; % Double-summation counter
    for n = 1:coeffs.N
        for m = n+1:coeffs.N
            for pm = [1 -1] % Loop over both sum & difference (pm = +/-)
                cnm = cnm + 1; % Update counter
                omega = coeffs.omega_npm(cnm); kappa = coeffs.kappa_npm(cnm); 
                kx = coeffs.kx_npm(cnm); ky = coeffs.ky_npm(cnm);         
                F = coeffs.F_npm(cnm); A = coeffs.A_npm(cnm); B = coeffs.B_npm(cnm);
                [phi,u,v,w,p, uV,vV,a_x,a_y] = addKinematics(phi,u,v,w,p, uV,vV,a_x,a_y, omega,kx,ky,kappa,Z,F,A,B, x,y,t);                
            end % Ends pm loop
        end % Ends m loop
    end % Ends n loop
end % End, if order>=2

% Third order
if order == 3
    % Single summations
    for n = 1:coeffs.N 
        % Third-order corrections to first-order terms
        omega = coeffs.omega(n); kappa = coeffs.kappa(n); 
        kx = coeffs.kx(n); ky = coeffs.ky(n);         
        F = coeffs.F13(n); A = coeffs.a(n); B = coeffs.b(n);
        [phi,u,v,w,p, uV,vV,a_x,a_y] = addKinematics(phi,u,v,w,p, uV,vV,a_x,a_y, omega,kx,ky,kappa,Z,F,A,B, x,y,t);

        % Self-self-self interactions
        omega = 3*coeffs.omega(n); kappa = coeffs.kappa_3(n); % !!!
        kx = 3*coeffs.kx(n); ky = 3*coeffs.ky(n);         
        F = coeffs.F_3(n); A = coeffs.A_3(n); B = coeffs.B_3(n);
        [phi,u,v,w,p, uV,vV,a_x,a_y] = addKinematics(phi,u,v,w,p, uV,vV,a_x,a_y, omega,kx,ky,kappa,Z,F,A,B, x,y,t);
    end
    
    % Double summations
    cnm = 0; % Double-summation counter
    for n = 1:coeffs.N 
        for m = n+1:coeffs.N
            for pm = [1 -1] % +/- m
                cnm = cnm + 1; % Update counter
                
                % n +/- 2m quantities 
                omega = coeffs.omega_np2m(cnm); kappa = coeffs.kappa_np2m(cnm); 
                kx = coeffs.kx_np2m(cnm); ky = coeffs.ky_np2m(cnm);         
                F = coeffs.F_np2m(cnm); A = coeffs.A_np2m(cnm); B = coeffs.B_np2m(cnm);
                [phi,u,v,w,p, uV,vV,a_x,a_y] = addKinematics(phi,u,v,w,p, uV,vV,a_x,a_y, omega,kx,ky,kappa,Z,F,A,B, x,y,t);

                % 2n +/- m quantities 
                omega = coeffs.omega_2npm(cnm); kappa = coeffs.kappa_2npm(cnm); 
                kx = coeffs.kx_2npm(cnm); ky = coeffs.ky_2npm(cnm);         
                F = coeffs.F_2npm(cnm); A = coeffs.A_2npm(cnm); B = coeffs.B_2npm(cnm);
                [phi,u,v,w,p, uV,vV,a_x,a_y] = addKinematics(phi,u,v,w,p, uV,vV,a_x,a_y, omega,kx,ky,kappa,Z,F,A,B, x,y,t);
            end
        end
    end % End of double summation
      
    % Triple summations
    c3 = 0; % Initialize counter
    for n = 1:coeffs.N 
        for m = n+1:coeffs.N
            for pmm = [1 -1] % +/- m
                for q = m+1:coeffs.N
                    for pmp = [1 -1] % +/- p
                        c3 = c3 + 1; % Update counter for (n +/- m +/- p) quantities
                        omega = coeffs.omega_npmpp(c3); kappa = coeffs.kappa_npmpp(c3); 
                        kx = coeffs.kx_npmpp(c3); ky = coeffs.ky_npmpp(c3);         
                        F = coeffs.F_npmpp(c3); A = coeffs.A_npmpp(c3); B = coeffs.B_npmpp(c3);
                        [phi,u,v,w,p, uV,vV,a_x,a_y] = addKinematics(phi,u,v,w,p, uV,vV,a_x,a_y, omega,kx,ky,kappa,Z,F,A,B, x,y,t);
                   end % End of pmp loop
                end % End of p loop
            end % End of pmm loop
        end % End of m loop
    end % End of triple summation
end % End of third order

end % End of function


%%% Internal function
% A function to generally add new kinematic contributions
function [phi,u,v,w,p, uV,vV,a_x,a_y] = addKinematics(phi,u,v,w,p, uV,vV,a_x,a_y, omega,kx,ky,kappa,Z,F,A,B, x,y,t)
    theta = omega*t - kx*x - ky*y; % Phase function
    factorZ = F*cosh(kappa*Z); % Common factor
    phiAdd = factorZ.*(A*sin(theta) - B*cos(theta)); phi = phi + phiAdd; % Velocity potential
    uAdd = kx*factorZ.*(-A*cos(theta) - B*sin(theta)); u = u + uAdd; % x velocity, u = d(phi)/dx
    vAdd = ky*factorZ.*(-A*cos(theta) - B*sin(theta)); v = v + vAdd; % y velocity, v = d(phi)/dy
    wAdd = F.*sinh(kappa*Z).*kappa.*(A*sin(theta) - B*cos(theta)); w = w + wAdd;% z velocity, w = d(phi)/dz
    pAdd = -factorZ.*omega.*(A*cos(theta) + B*sin(theta)); p = p + pAdd; % Dynamic pressure p = p^+/rho = - d(phi)/dt
    
    % Drag forces
    V = sqrt(uAdd.^2 + vAdd.^2); uV = uV + uAdd.*V; vV = vV + vAdd.*V;
    
    % Inertial forces
    if omega == 0, omega = eps; end % Avoids divide by zero below
    T = tanh(kappa*Z).*pAdd;
    u_t = kx*omega*phiAdd; v_t = ky*omega*phiAdd; % Time derivatives
    u_x = -kx^2*phiAdd; v_x = -kx*ky*phiAdd; % x derivatives
    u_y = -kx*ky*phiAdd; v_y = -ky^2*phiAdd; % y derivatives
    u_z = kx*kappa/omega*T; v_z = ky*kappa/omega*T; % z derivatives
    a_x = a_x + u_t + u.*u_x + v.*u_y + w.*u_z; % Total x acceleration (including convective acceleration)
    a_y = a_y + v_t + u.*v_x + v.*v_y + w.*v_z; % Totay y acceleration
end