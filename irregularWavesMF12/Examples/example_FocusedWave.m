% Demonstrates how to create time series and spatial plots for
% directionally spread irregular waves, based on a JONSWAP spectrum.
clear all, close all, addpath ../Source;

% Basic input
A = 15.7; % Linear amplitude sum (m)
kph = 1.2; % Peak wave number times water depth
h = 70; % Water depth (m)
g = 9.81; % Gravitational acceleration (m/s^2)
Ux = 0; Uy = 0; % Current (m/s)
x = 0; y = 0; % Position for time series (m)

% JONSWAP spectrum input
gamma = 3.3; % Peak enhancement factor
kp = kph/h; % Peak wave number (1/m)
omegap = sqrt(g*kp*tanh(kph)); % Peak angular frequency
Tp = 2*pi/omegap; % Peak period
Nfreq = 301; fcut = 3/Tp;
%tdur = 100*Tp; % Duration of time series
%dt = Tp/40; N = round(tdur/dt); % Time step (s), and length of series
%fN = 1/(2*dt); df = fN/(N/2); fcut = 3/Tp; % Nyquist frequency (Hz), frequency increment (Hz), cut-off frequency (Hz)
%Nfreq = round(fcut/df+1); % Number of discrete frequencies
%t = [0:dt:tdur-dt]; % Time duration (s)
ND = 10; s = ND/2; % Directional spreading parameter
seed = 1963; % Seed for random number generator
%t = 0; -2*Tp; % Time (s)

% Generate harmonic amplitudes and wave number components for the JONSWAP spectrum
rng(seed) % Re-set the random number generator with the seed
[a,b,kx,ky, fJ, SJ] = jonswap(g,h,1,Tp,fcut,s,gamma,Nfreq);

% Adjust to create a focused wave
a = sqrt(a.^2 + b.^2); % Eliminate phase difference
a = [a a];  b = 0*a; kx = [kx kx]; ky = [ky -ky]; % Reflect to make symmetric
a = a/sum(a)*A; % Rescale to give desired linear amplitude sum

% Create a spatial mesh
Lp = 2*pi/kp; % Peak wave length
Lx = 4*Lp; Ly = 0.4*sqrt(1+ND)*Lx; 
dx = Lp/20; dy = dx*sqrt(1+ND);
x = [-Lx/2:dx:Lx/2]; y = [-Ly/2:dy:Ly/2];
[X,Y] = meshgrid(x,y);

% Calculate coefficients (to second order)
disp(' '), disp('Generating wave field ...')
coeffs = coeffsMF12(2,g,h,a,b,kx,ky,Ux,Uy); % Determine coefficients

% Animate the free surface (use first order here, it is much faster)
disp('Animating free surface to first order ...')
figure(), jc = (length(y)-1)/2+1; % y=0 grid point
for t = -2*Tp:Tp/20:0
    [eta1] = surfaceMF12(1,coeffs,X,Y,t); % Obtain the wave field
    
    % Plot the surface
    subplot(2,1,1), surf(X,Y,eta1), axis equal, shading interp
    zlim([-0.5 1]*A), set(gca,'DataAspectRatio',[1 1 0.05])
    xlabel('x (m)'), ylabel('y (m)'), zlabel('\eta (m)')
    colormap jet, title('Focused wave field, first order')
    hcb = colorbar; caxis([-0.5 1]*A), title(hcb,'\eta (m)')
    
    % Plot along the centerline, y=0
    subplot(2,1,2), plot(x,eta1(jc,:),'k-'), grid on
    xlabel('x (m)'), ylabel('\eta (m)'), title('Surface elevation along y=0')
    ylim([-1 1]*A), drawnow
end

% Also calculate the surface to second order at the end
disp('Computing focused wave to second order ...')
[eta2] = surfaceMF12(2,coeffs,X,Y,t);
figure(), subplot(3,1,1), surf(X,Y,eta2), axis equal, shading interp
set(gca,'DataAspectRatio',[1 1 0.05])
xlabel('x (m)'), ylabel('y (m)'), zlabel('\eta (m)')
colormap jet, title('Focused wave field, Second order')
hcb = colorbar; title(hcb,'\eta (m)'), drawnow

% Compare the results along the centerline
disp('Comparing results along the centerline y=0 ...')
subplot(3,1,2), plot(x,eta1(jc,:),'k--'), grid on
hold on, plot(x,eta2(jc,:),'b-'), hold off
legend('First order','Second order')
title('Surface elevation along y=0')

% Compute the velocity profiles
disp('Computing velocity profiles beneath the crest ...')
eta01 = surfaceMF12(1,coeffs,0,0,0); z1 = linspace(-h,eta01,100);
eta02 = surfaceMF12(2,coeffs,0,0,0); z2 = linspace(-h,eta02,100);
[u1] = kinematicsMF12(1,coeffs,0,0,z1,0);
[u2] = kinematicsMF12(2,coeffs,0,0,z2,0);
subplot(3,1,3), plot(u1,z1,'k--'), ylim([-h eta02])
hold on, plot(u2,z2,'b-'), hold off, grid on
xlabel('u (m/s)'), ylabel('z (m)'), title('Velocity profile')