function p = probDensFun(eta,bins)
% Creates a probability density function from a time or space series (eta) 
% with discrete bins, following the example from Sumer & Fuhrman (2020), 
% p. 169.
%
% Reference:
%
% Sumer, B.M. & Fuhrman, D.R. (2020) Turbulence in Coastal and Civil
% Engineering. World Scientific.
%
p = hist(eta,bins); % Creates a histogram
area = trapz(bins,p); % Integrate histogram
p = p/area; % Normalize to yield PDF with unit area
end