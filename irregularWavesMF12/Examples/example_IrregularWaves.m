% Demonstrates how to create time series and spatial plots for
% directionally spread irregular waves, based on a JONSWAP spectrum.
clear all, close all, addpath ../Source;

% Basic input
order = 2; % Order of theory
steep = 0.15; % Characteristic steepness, kp*Hm0/2
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

% Generate harmonic amplitudes and wave number components for the JONSWAP spectrum
rng(seed) % Re-set the random number generator with the seed
[a,b,kx,ky, fJ, SJ] = jonswap(g,h,Hm0,Tp,fcut,s,gamma,Nfreq);

% Generate the time series
disp(' '), disp('Generating time series ...'), tic
coeffs = coeffsMF12(order,g,h,a,b,kx,ky,Ux,Uy); % Determine coefficients
[eta] = surfaceMF12(order,coeffs,x,y,t); % Obtain the time series
toc

% Plot the time series
figure(1), subplot(3,1,1), plot(t,eta), xlim([0 tdur])
xlabel('t (s)'), ylabel('\eta (m)'), grid on
title('Second-order irregular wave time series')

% Compute and plot the spectrum based on the time series
disp('Computing spectrum ...')
[f,S] = spectrum(eta,dt); 
subplot(3,1,2), plot(f,S,'b-'), title('Spectrum')
xlabel('f (Hz)'), ylabel('S(f) (m^2 s)'), grid on
areaS = trapz(f,S); % Integrate the spectrum ~(Hm0/4)^2
hold on, plot(fJ,SJ,'k--'), hold off % JONSWAP spectrum
legend('Computed','JONSWAP')

% Calculate basic statistics from the time series
disp(' '), disp('Statistics from time series:')
etaBar = mean(eta); disp(['mean(eta) = ' num2str(etaBar) ' m'])
sigma2 = var(eta); sigma = sqrt(sigma2); disp(['variance = ' num2str(sigma2) ' m^2'])
disp(['(Hm0/4)^2 = ' num2str((Hm0/4)^2) ' m^2'])
disp(['Integral of spectrum = ' num2str(areaS) ' m^2'])
Sk = skewness(eta); disp(['skewness = ' num2str(Sk)])
K = kurtosis(eta); disp(['kurtosis = ' num2str(K)]), disp(' ')

% Create a probability density function
disp('Computing probability density function ...')
bins = sigma*[-6:0.2:6]; % Bins for eta (m)
p = probDensFun(eta,bins); subplot(3,1,3) % Compute probability density (1/m)
semilogy(bins/sigma,p*sigma,'b.','MarkerSize',8)
xlabel('\eta/\sigma'), ylabel('p(\eta)\sigma'), grid on, ylim([1e-8 1])
pG = 1/sigma/sqrt(2*pi)*exp(-bins.^2/2/sigma2); % Gaussian distribution
hold on, plot(bins/sigma,pG*sigma,'k--'), hold off
legend('Computed (time series)','Gaussian','Location','South')
title('Probability density function'), drawnow


%%% Generate a spatial plot

% Generate a spatial map
t = 0; % Set time to a scalar

% Generate a spatial mesh
Lp = 2*pi/kp; % Peak wave length
Lx = 20*Lp; Ly = 0.1*sqrt(1+ND)*Lx; 
dx = Lp/20; dy = dx*sqrt(1+ND);
x = [0:dx:Lx]; y = [0:dy:Ly];
[X,Y] = meshgrid(x,y);

% Generate the free surface time series
disp(' '), disp('Generating wave field ...'), tic
%order = 1; % Uncomment to use order=1 (much faster, order=2 takes ~1 min)
[eta2] = surfaceMF12(order,coeffs,X,Y,t); % Obtain the wave field
toc

% Plot the free surface
figure(3), surf(X,Y,eta2), axis equal, shading interp
xlabel('x (m)'), ylabel('y (m)'), zlabel('\eta (m)')
colormap jet, title('Second-order irregular wave field')
hcb = colorbar; title(hcb,'\eta (m)')

% Add results to the PDF
p2 = probDensFun(eta2(:),bins); 
figure(1), subplot(3,1,3)
semilogy(bins/sigma,p*sigma,'b.','MarkerSize',8)
xlabel('\eta/\sigma'), ylabel('p(\eta)\sigma'), grid on, ylim([1e-8 1])
hold on, plot(bins/sigma,p2*sigma,'r.','MarkerSize',8)
pG = 1/sigma/sqrt(2*pi)*exp(-bins.^2/2/sigma2); % Gaussian distribution
plot(bins/sigma,pG*sigma,'b--'), hold off

% Add PDF from Fuhrman et al. (2023)
Sk = skewness(reshape(eta2,1,[])); % Skewness
pF = pdfAiry(bins,Sk,1); % PDF of Fuhrman et al. (2023), with modified negative tail
hold on, plot(bins,pF,'k-'), hold off
legend('Computed (time series)','Gaussian','Location','South')
title('Probability density function'), drawnow
legend('Computed (time series)','Computed (wave field)','Gaussian','Fuhrman et al. (2023)','Location','South')