% An example generating a doubly-periodic short-crested wave field
clear all, close all, addpath ../Source;

% Input
order = 3; % Order of theory
h = 2*pi; g = 9.81; % Water depth (m), gravitational acceleration (m/s^2)
a = [0.2 0.2 0]; b = [0 0 0]; % Harmonic amplitudes
Ux = 0; Uy = 0; % Current vectory (m/s)
t = 0; % Time (s)
kappa = [1 1 1]; % Wave number magnitudes (1/m)
phi = [15 -15 0]*pi/180; % Angle
kx = kappa.*cos(phi); ky = kappa.*sin(phi); % Wave number vector components (1/m)

% Create a spatial mesh
Lx = 2*pi/kx(1); Ly = 2*pi/ky(1);
x = [-1.5*Lx:Lx/20:1.5*Lx]; y = [-Ly/2:Ly/20:Ly/2];
[X,Y] = meshgrid(x,y);

% Run the theory
coeffs = coeffsMF12(order,g,h,a,b,kx,ky,Ux,Uy); % Determine coefficients
[eta,phiS,Mx,My] = surfaceMF12(order,coeffs,X,Y,t); % Obtain the wave field

% Plot the free surface
surf(X,Y,eta)
xlabel('x (m)'), ylabel('y (m)'), zlabel('\eta (m)')
axis equal, shading interp, title('Third-order short-crested waves')
colormap jet;

% Demonstrate 3D velocity kinematics
z = [-h:h/20:0]; [X,Y,Z] = meshgrid(x,y,z); % Generate 3D grid
[u,v,w] = kinematicsMF12(order,coeffs,X,Y,Z,t); % Third-order velocity field
figure(), slice(X,Y,Z,u,[],[0],[-1 -2]) % Make a slice plot for u
axis tight, axis equal, shading interp, colormap jet
xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
hcb = colorbar; title(hcb,'u (m/s)')
title('3D velocity kinematics')