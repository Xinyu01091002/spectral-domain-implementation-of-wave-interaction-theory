% An example comparing a first-, second-, and third-order standing wave
clear all, close all, addpath ../Source;

% Input
order = 3;
h = 2*pi; g = 9.81; % Water depth (m), gravitational acceleration (m/s^2)
a = [0.25 0.25 0]; b = [0 0 0]; % Harmonic amplitudes (requires 3, but can be zeros!)
%a = [0.25 0.1 0]; b = [0 0 0]; % Un-comment for a partial standing wave
kx = [1 -1 0.02]; ky = 0*kx; % Wave numbers (1/m)
Ux = 0; Uy = 0; % Current velocities (m/s)

% Create spatial mesh
y = 0; Lx = 2*pi/kx(1); X = [0:Lx/50:Lx];

% Determine coefficients
coeffs3 = coeffsMF12(3,g,h,a,b,kx,ky,Ux,Uy); % Third order

% Animate
omega = coeffs3.omega(1); T = 2*pi/omega; % Determine wave period
figure(1)
for t = 0:T/100:5*T
    [eta1] = surfaceMF12(1,coeffs3,X,y,t); % First order
    [eta2] = surfaceMF12(2,coeffs3,X,y,t); % Second order    
    [eta3] = surfaceMF12(3,coeffs3,X,y,t); % Third order
    plot(X,eta1,'k--'), xlabel('x (m)'), ylabel('\eta (m)')
    hold on, plot(X,eta2,'r--'), hold off
    hold on, plot(X,eta3,'b-'), hold off
    xlim([0 X(end)]), ylim(3*a(1)*[-1 1])
    title('First-, second- & third-order standing wave')
    %legend('First-order','Second-order','Third-order','Location','SouthWest') % Slows animation
    grid on, drawnow
end