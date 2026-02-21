% Example demonstrating velocity and acceleration kinematics, as well as
% resulting forces and moments based on the Morison equation.
clear all, close all, addpath ../Source;

% Basic input
order = 2; % Order of theory
steep = 0.10; % Characteristic steepness, kp*Hm0/2
kph = 1.2; % Peak period times water depth
h = 70; % Water depth (m)
g = 9.81; % Gravitational acceleration
Ux = 0; Uy = 0; % Current (m/s)
x = 0; y = 0; % Position (m)
rho = 1000; D = 2; cD = 1; cM = 1; % Density (kg/m^3), monopile diameter (m), drag and inertial coefficient

% JONSWAP spectrum input
gamma = 3.3; % Peak enhancement factor
kp = kph/h; % Peak wave number (1/m)
Hm0 = 2*steep/kp; % Spectral significant wave height (m)
omegap = sqrt(g*kp*tanh(kph)); % Peak angular frequency
Tp = 2*pi/omegap; % Peak period
tdur = 10*Tp; % Duration of time series
dt = Tp/40; N = round(tdur/dt); % Time step (s), and length of series
fN = 1/(2*dt); df = fN/(N/2); fcut = 3/Tp; % Nyquist frequency (Hz), frequency increment (Hz), cut-off frequency (Hz)
Nfreq = round(fcut/df+1); % Number of discrete frequencies
t = [0:dt:tdur-dt]; % Time duration (s)
ND = 2; s = ND/2; % Directional spreading parameter
seed = 1963; % Seed for random number generator

% Generate harmonic amplitudes and wave number components for the JONSWAP spectrum
rng(seed) % Re-set the random number generator with the seed
[a,b,kx,ky] = jonswap(g,h,Hm0,Tp,fcut,s,gamma,Nfreq);

% Generate the time series
disp(' '), disp('Generating time series ...'), tic
coeffs = coeffsMF12(order,g,h,a,b,kx,ky,Ux,Uy); % Determine coefficients
[eta] = surfaceMF12(order,coeffs,x,y,t); % Obtain the time series
toc

% Plot the time series
figure(1), subplot(3,4,[1:4]), plot(t,eta), xlim([0 tdur]), grid on
xlabel('t (s)'), ylabel('\eta (m)'), title('Second-order irregular wave time series')
t0 = t(1); eta0 = eta(1); 
hold on, p = plot(t0,eta0,'r.','MarkerSize',8); hold off, 
p.XDataSource = 't0'; p.YDataSource = 'eta0';

% Animate the velocity field and force and momentum vectors
disp('Animating kinematics ...')
scale = 1;
K_D = rho*cD*D/2/1000; K_I = pi/4*rho*D^2*cM/1000; % Coefficient for drag and inertial forces
for j = 1:length(eta)
    t0 = t(j); eta0 = eta(j); % Current time and surface elevation
    z = linspace(-h,eta0,100); % z vector from bottom to surface
    [u,v,w,p,phi, uV,vV,a_x,a_y] = kinematicsMF12(order,coeffs,x,y,z,t0);
    
    % Plot velocity components
    subplot(3,4,5:6), plot(u,z), title('Velocities')
    xlim(scale*[-1 1]*Hm0*omegap), ylim([-h Hm0]), grid on
    hold on, plot(v,z,'r-'), plot(w,z,'k-'), hold off
    xlabel('(u,v,w) (m/s)'), ylabel('z (m)')
    legend('u','v','w','Location','SouthWest')
    
    % Plot acceleration components
    subplot(3,4,7:8), plot(a_x,z,'b-'), xlabel('(a_x,a_y) (m/s^2)')
    xlim(scale*[-1 1]*Hm0*omegap^2), ylim([-h Hm0]), grid on
    hold on, plot(a_y,z,'r-'), hold off, title('Accelerations')
    legend('a_x','a_y','Location','SouthWest')
    
    % Plot resultant forces
    [F_x,F_y,M_x,M_y, F_Dx,F_Dy,F_Ix,F_Iy, M_Dx,M_Dy,M_Ix,M_Iy] = morison(z,h,uV,vV,a_x,a_y,K_D,K_I); % Compute forces and moments
    subplot(3,4,9:10)
    quiver(0,0,F_Dx,F_Dy,1,'g'), hold on % Drag component
    quiver(F_Dx,F_Dy,F_Ix,F_Iy,1,'m') % Inertial component
    quiver(0,0,F_x,F_y,1,'b'), hold off, grid on % Resultant force
    scaleF = (K_D*h*Hm0*omegap + K_I*h*Hm0*omegap^2); % Scaling for forces
    axis equal, xlim(0.25*scale*scaleF*[-1 1]), ylim(0.25*scale*scaleF*[-1 1])
    xlabel('F_x (kN)'), ylabel('F_y (kN)'), title('Forces')
    legend('Drag','Inertial','Total','Location','EastOutside')
    
    % Plot resultant moments
    subplot(3,4,11:12), quiver(0,0,M_Dx,M_Dy,1,'g') % Drag component
    hold on, quiver(M_Dx,M_Dy,M_Ix,M_Iy,1,'m') % Inertial component
    quiver(0,0,M_x,M_y,1,'b'), hold off, grid on % Resultant moment
    scaleM = (K_D*h^2*Hm0*omegap + K_I*h^2*Hm0*omegap^2); % Scaling for moments
    axis equal, xlim(0.25*scale*scaleM*[-1 1]), ylim(0.25*scale*scaleM*[-1 1])
    xlabel('M_x (kN m)'), ylabel('M_y (kN m)'), title('Moments')
    refreshdata, drawnow
end