% Example creating a spatial plot of the free surface and underlying 
% horizontal velocity profile using a basic three-wave input.  This case 
% corresponds to that given by Madsen & Fuhrman (2012) for checking their 
% irregular wave theory, based on their Fig. 1 with transfer function and 
% other coefficient values provided in their in Tables 2 and 3.
clear all, close all, addpath ../Source;

% Input
order = 3; % Order of theory
h = 1; g = 9.81; % Water depth (m), gravitational acceleration (m/s^2)
Ux = 0; Uy = 0; % Current velocities (m/s)
t = 0; % Time (s)
a = [0.025 0.025 0.05]; b = 0*a; % Harmonic amplitudes (m)
hkappa = [1.71666 1.85737 2.02010]; kappa = hkappa/h; % Exact kh values from MF12, p. 319
phi = [30 -20 0]*pi/180; % Angle
kx = kappa.*cos(phi); ky = kappa.*sin(phi); % Wave number vector components (1/m)

% Create a spatial mesh
x = [-10:0.1:10]; y = [-10:0.1:10];
[X,Y] = meshgrid(x,y);

% Run the theory with no current
coeffs1 = coeffsMF12(1,g,h,a,b,kx,ky,Ux,Uy); % Linear coefficients
coeffs = coeffsMF12(order,g,h,a,b,kx,ky,Ux,Uy); % Third-order coefficients
[eta3,PhiS,Mx,My] = surfaceMF12(order,coeffs,X,Y,t); % Third-order wave field (with Eulerian drift allowed)

% Add a return current, and re-run third-order theory
Ux = -Mx/h; Uy = -My/h; % Current to remove Eulerian drift
coeffs = coeffsMF12(order,g,h,a,b,kx,ky,Ux,Uy,1); % The final (optional) argument of 1 displays the coefficients
[eta3,PhiS,Mx,My] = surfaceMF12(order,coeffs,X,Y,t); % Third-order wave field 
[eta1] = surfaceMF12(1,coeffs1,X,Y,t); % First-order wave field

% Also obtain kinematics beneath the crest
x0 = 0; y0 = 0; t0 = 0;
i0 = find(x==x0); j0 = find(y==y0);
eta0 = eta3(i0,j0); z = linspace(-h,eta0,100);
[u3] = kinematicsMF12(order,coeffs,x0,y0,z,t0); % Third-order velocity profile
[u1] = kinematicsMF12(1,coeffs1,x0,y0,z,t0); % First-order velocity profile

% Plot the third-order wave field (Compare with Fig. 1a of MF12)
subplot(3,1,1), surf(X,Y,eta3), title('Third-order wave field')
xlabel('x (m)'), ylabel('y (m)'), zlabel('\eta (m)')
axis equal, shading interp

% Plot a line along y=0 (Compare with Fig. 1b of MF12)
subplot(3,1,2), plot(x,eta1(j0,:),'k--')
hold on, plot(x,eta3(j0,:),'b-'), hold off
xlabel('x (m)'), ylabel('\eta (m)'), ylim([-1 1].*0.1501)
set(gca,'YTick',[-0.15:0.05:0.15]), grid on
title('Surface elevation along y=0')

% Plot the horizontal velocity profile beneath the crest (compare with Fig. 1c of MF12)
subplot(3,1,3), plot(u3,z,'b-'), xlabel('u (m/s)'), ylabel('z (m)')
hold on, plot(u1,z,'k--'), hold off
xlim([0 0.6]), ylim([-h 0.2]), grid on
title('Velocity profile beneath the crest at (x,y)=(0,0)')
legend('Third-order','First-order','Location','SouthEast')