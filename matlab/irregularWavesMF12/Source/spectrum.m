function [f,S] = spectrum(eta,dt)
% Computes a spectrum from a given time series of the free surface elevation 
% (eta) with fixed time step (dt). The function returns a vector of discrete
% frequencies (f, in Hz) and spectral density S (m^2*s), where SI units 
% have been assumed. The integral of the spectrum (trapz(f,S))should match 
% the variance of the time series, provided it has zero mean. This code
% essentially follows that of Sumer & Fuhrman (2020), p. 216.
%
% Reference:
%
% Sumer, B.M. & Fuhrman, D.R. (2020) Turbulence in Coastal and Civil
% Engineering. World Scientific.
%
N = length(eta); % Length of time series
fN = 1/(2*dt); % Nyquist frequency
df = fN/(N/2); % Discrete frequency step
f = [0:df:fN]; % Vector of discrete frequencies
Eta = fft(eta)/N; % Take the Fourier transform
A = abs([Eta(1) 2*Eta(2:N/2) Eta(N/2+1)]); % Harmonic amps (m)
S = 0.5*A.^2/df; % Spectral density (m^2*s)
end