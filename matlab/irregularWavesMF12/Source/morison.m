function [F_x,F_y,M_x,M_y, F_Dx,F_Dy,F_Ix,F_Iy, M_Dx,M_Dy,M_Ix,M_Iy] = morison(z,h, uV,vV,a_x,a_y, K_D,K_I)
% Integrates kinematic quantities to obtain drag and inertial forces and 
% moments based on the vectorized Morison (Morison et al. 1950) equation.
%
% Reference:
%
% Morison, J.R., O'Brien, M.P., Johnson, J.W. & Schaaf, S.A. (1950) The
% force exerted by surface waves on pile. J. Petrol. Tech. 2, 149-154.

% Forces
F_Dx = trapz(z,K_D.*uV); F_Dy = trapz(z,K_D.*vV); % Drag force components F_D = (F_Dx,F_Dy)
F_Ix = trapz(z,K_I.*a_x); F_Iy = trapz(z,K_I.*a_y); % Inertial force components F_I = (F_Ix,F_Iy)
F_x = F_Dx + F_Ix; F_y = F_Dy + F_Iy; % Total force components F = (F_x,F_y)

% Moments (about bottom z=-h)
Z = z + h; % Shift coordinate system such that Z=0 is the sea bed
M_Dx = -trapz(z,K_D.*Z.*vV); M_Dy = trapz(z,K_D.*Z.*uV); % Drag moment components M_D = (M_Dx,M_Dy)
M_Ix = -trapz(z,K_I.*Z.*a_y); M_Iy = trapz(z,K_I.*Z.*a_x); % Inertia moment components M_I = (M_Ix,M_Iy)
M_x = M_Dx + M_Ix; M_y = M_Dy + M_Iy; % Total moment components M = (M_x,M_y)
end