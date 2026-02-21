% An example comparing a first-, second-, and third-order Stokes wave
clear all, close all, addpath ../Source;

% Input
h = 2*pi; g = 9.81; % Water depth (m), gravitational acceleration (m/s^2)
a = [0.45 0 0]; b = [0 0 0]; % Harmonic amplitudes (requires minimum of 3, but can be zero!)
kx = [1 0.01 0.02]; ky = 0*kx; % Wave numbers (1/m)
Ux = 0; Uy = 0; % Current velocities (m/s)
t = 0; % Initial time (s)

% Create spatial mesh
y = 0; Lx = 2*pi/kx(1); X = [-Lx:Lx/50:Lx];

% Second order (adjusting return current to give zero Stokes drift)
coeffs2 = coeffsMF12(2,g,h,a,b,kx,ky,Ux,Uy); % Determine second-order coefficients
[eta2,phiS2,Mx,My] = surfaceMF12(2,coeffs2,X,y,t); % Obtain the wave field
disp(' '), disp(['Mass flux without return current: (Mx,My) = (' num2str(Mx) ',' num2str(My) ') m^2/s'])
Ux = -Mx/h; Uy = -My/h; % Add return current
coeffs2 = coeffsMF12(2,g,h,a,b,kx,ky,Ux,Uy); % Re-run with return current added
[eta2,phiS,Mx,My] = surfaceMF12(2,coeffs2,X,y,t); % Obtain the wave field
disp(['Mass flux after adding return current: (Mx,My) = (' num2str(Mx) ',' num2str(My) ') m^2/s'])

% First order (for comparison)
coeffs1 = coeffsMF12(1,g,h,a,b,kx,ky,Ux,Uy); % Determine coefficients
[eta1] = surfaceMF12(1,coeffs1,X,y,t); % Obtain the wave field

% Third order
coeffs3 = coeffsMF12(3,g,h,a,b,kx,ky,Ux,Uy); % Determine coefficients
[eta3] = surfaceMF12(3,coeffs3,X,y,t); % Obtain the wave field

% Plot the initial (t=0) free surfaces
figure(1), plot(X,eta1,'k--'), xlabel('x (m)'), ylabel('\eta (m)')
hold on, plot(X,eta2,'r--'), plot(X,eta3,'b-'), hold off, grid on
xlim([X(1), X(end)]), title('First-, second- & third-order Stokes wave')
legend('First-order','Second-order','Third-order')

% Animate (demonstrates amplitude dispersion)
omega = coeffs3.omega(1); T = 2*pi/omega;
figure(2), disp(' '), disp('Animating ...')
for t = 0:T/80:3*T
    [eta1] = surfaceMF12(1,coeffs1,X,y,t); % First order
    [eta2] = surfaceMF12(2,coeffs2,X,y,t); % Second order
    [eta3] = surfaceMF12(3,coeffs3,X,y,t); % Third order
    plot(X,eta1,'k--'), xlabel('x (m)'), ylabel('\eta (m)')
    hold on, plot(X,eta2,'r--'), plot(X,eta3,'b-'), hold off, grid on
    xlim([X(1) X(end)]), title('First-, second- & third-order Stokes wave')
    %legend('First-order','Second-order','Third-order','Location','SouthWest') % Slows animation
    drawnow
end