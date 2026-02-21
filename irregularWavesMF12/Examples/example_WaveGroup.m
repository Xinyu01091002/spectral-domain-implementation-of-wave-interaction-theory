% An example animating a wave group
clear all, close all, addpath ../Source;

% Input
h = 2*pi; g = 9.81; % Water depth (m), gravitational acceleration (m/s^2)
a = [0.01 0.02 0.01 0.005 0.005]*0.1; b = 0.*a; % Harmonic amplitudes (requires 3, but can be zeros!)
Ux = 0; Uy = 0; % Current vectory (m/s)
t=0; % Time (s)
kx = [1 1.05 0.95 1.1 0.9]; ky = 0*kx; % Wave number magnitudes (1/m)

% Create spatial mesh
y = 0; Lx = 2*pi/kx(1); X = [0:Lx/50:2*Lx];

% Compute coefficients
order = 1; coeffs1 = coeffsMF12(order,g,h,a,b,kx,ky,Ux,Uy);
order = 2; coeffs2 = coeffsMF12(order,g,h,a,b,kx,ky,Ux,Uy);
order = 3; coeffs3 = coeffsMF12(order,g,h,a,b,kx,ky,Ux,Uy);
% Calculate group length and period, and define spatial domain
Lg = 4*pi/abs(kx(2) - kx(1)); % Group length
Tg = 4*pi/abs(coeffs2.omega(2) - coeffs2.omega(1)); % Group period
X = [0:Lg/4000:Lg]; y = 0; % Spatial domain
N = length(X);
% Lg 是 X 向量的总长度
% 波数向量（用于 FFT）
k_unshiffed = (2*pi/Lg) * [0:floor(N/2), -floor(N-1)/2:-1];
% 如果您希望波数从负到正排列（中心为零），可以使用 fftshift

% Animate surface elevation
tvec = [0:Tg/1000:Tg];
for t = tvec
    order = 2; [eta2] = surfaceMF12(order,coeffs2,X,y,t); % Third order
    plot(X,eta2,'r--'), xlabel('x (m)'), ylabel('\eta (m)')
    order = 1; [eta1] = surfaceMF12(order,coeffs1,X,y,t); % First order
    VWAeta33=real((ifft(fft(hilbert(eta1)).*(k_unshiffed))).^2.*hilbert(eta1)*3/8);
%     hold on 
%     plot(X,eta1,'k')
%     plot(X,VWAeta33,'r--');
    hold on, plot(X,eta1,'k'), hold off
    order = 3; [eta3] = surfaceMF12(order,coeffs3,X,y,t); % First order
    hold on, plot(X,eta3,'b--'), hold off
    xlim([X(1) X(end)]), title('First- and second-order wave groups')
    drawnow
    figure 
%     plot(eta1);
    hold on 
%     plot(eta2-eta1);
    plot(eta3-eta2,DisplayName='Fuhrman');
    plot(VWAeta33,'r--',DisplayName='VWA');
    legend
end