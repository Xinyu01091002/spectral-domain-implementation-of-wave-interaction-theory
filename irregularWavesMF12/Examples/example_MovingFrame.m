% Example comparing irregular, directionally-spread wave time series 
% for fixed and moving frames of reference.
clear all, close all, addpath ../Source;

% Basic input
order = 2; % Order of theory
steep = 0.10; % Characteristic steepness, kp*Hm0/2
kph = 1.2; % Peak wave number times water depth
h = 70; % Water depth (m)
g = 9.81; % Gravitational acceleration (m/s^2)
Ux = 0; Uy = 0; % Current (m/s)
x = 0; y = 0; % Position for time series (m)

% JONSWAP spectrum input
gamma = 3.3; % Peak enhancement factor
kp = kph/h; % Peak wave number (1/m)
Hm0 = 2*steep/kp; % Spectral significant wave height (m)
omegap = sqrt(g*kp*tanh(kph)); % Peak angular frequency
Tp = 2*pi/omegap; % Peak period
tdur = 100*Tp; % Duration of time series
dt = Tp/20; N = round(tdur/dt); % Time step (s), and length of series
fN = 1/(2*dt); df = fN/(N/2); fcut = 3/Tp; % Nyquist frequency (Hz), frequency increment (Hz), cut-off frequency (Hz)
Nfreq = round(fcut/df+1); % Number of discrete frequencies
t = [0:dt:tdur-dt]; % Time duration (s)
ND = 50; s = ND/2; % Directional spreading parameter
seed = 1963; % Seed for random number generator

% Velocity of moving frame of reference
uShip = 10; vShip = 5; % Ship velocity components (m/s)
x = uShip*t; y = vShip*t; % x and y positions for moving frame

% Generate harmonic amplitudes and wave number components for the JONSWAP spectrum
rng(seed) % Re-set the random number generator with the seed
[a,b,kx,ky] = jonswap(g,h,Hm0,Tp,fcut,s,gamma,Nfreq);

% Generate the time series, fixed frame
disp(' '), disp('Generating time series (fixed frame) ...'), tic
coeffs = coeffsMF12(order,g,h,a,b,kx,ky,Ux,Uy); % Determine coefficients
[eta0] = surfaceMF12(order,coeffs,0,0,t); % Obtain the time series
toc

% Generate the time series, moving frame
disp(' '), disp('Generating time series (moving frame) ...'), tic
coeffs = coeffsMF12(order,g,h,a,b,kx,ky,Ux,Uy); % Determine coefficients
[eta] = surfaceMF12(order,coeffs,x,y,t); % Obtain the time series
toc

% Plot the time series
figure(1), subplot(2,1,1), plot(t,eta0,'r-'), xlim([0 tdur]) % Fixed frame
ylabel('\eta (m)'), grid on
title('Fixed frame')
subplot(2,1,2), plot(t,eta,'b-'), xlim([0 tdur]) % Moving frame
xlabel('t (s)'), ylabel('\eta (m)'), grid on
title('Moving frame')
