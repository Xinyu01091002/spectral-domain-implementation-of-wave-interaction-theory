function [a,b,kx,ky, f,S] = jonswap(g,h,Hm0,Tp,fcut,s,gammaVal,Nfreq)
% Returns cos (a) and sin (a) amplitudes as well as wave number components 
% (kx,ky) % for a JONSWAP spectrum with cos^(2*s) directional spreading.
% Single summation is used for the directional spreading. Outputs include 
% the discrete frequency vector (f) and Spectrum function (S).

% Determine discrete frequencies, up to cut off
fp = 1/Tp; f = linspace(0,fcut,Nfreq); df = f(2); omega = 2*pi*f;

% Determine wave number moduli
for n = 1:Nfreq, k(n) = fzero(@(k) omega(n)^2-g*k*tanh(k*h),omega(n)^2/g); end

% Frequency spectrum
sigma = 0.09*ones(1,Nfreq); sigma(find(f<fp)) = 0.07;
S = 5*fp^4*f.^(-5).*exp(-5/4*fp^4*f.^(-4))*gammaVal^exp(-0.5*((f-fp)/(sigma*fp)).^2); S(1) = 0; 
areaS = trapz(f,S); S = S/areaS*(Hm0/4)^2; % Normalize
areaS = trapz(f,S); % =(Hm0/4)^2
psi = 2*pi*rand(1,Nfreq); % Random numbers between 0 and 2*pi
a = sqrt(2*S*df).*cos(psi); b = sqrt(2*S*df).*sin(psi); % cos/sin amps

% Directional spreading
if isinf(s) % Uni-directional
    theta = 0;
else
    thetavec = linspace(-pi,pi,Nfreq); 
    D = 1/(2*sqrt(pi))*gamma(1+s)/gamma(1/2+s)*cos(thetavec/2).^(2*s);
    %figure(), plot(thetavec/pi*180,D), xlim([-180 180]), xlabel('\theta (deg)'), ylabel('D(\theta)')
    areaD = trapz(thetavec,D); P = cumtrapz(thetavec,D); % AreaD=1.0
    i1 = min(find(P>=1)); % Index where P effectively reaches 1
    if isempty(i1), i1 = length(P); end
    Psi = rand(1,Nfreq); % Random numbers between 0 and 1
    theta = interp1(P(1:i1),thetavec(1:i1),Psi); % Vq = interp1(X,V,Xq)
end
kx = k.*cos(theta); ky = k.*sin(theta); % Wave number components

% Eliminate zero amplitudes
A = sqrt(a.^2 + b.^2); i0 = find(A>0); 
a = a(i0); b = b(i0); kx = kx(i0); ky = ky(i0);
end