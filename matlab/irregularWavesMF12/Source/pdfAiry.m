function p = pdfAiry(zeta,Sk,modNegTail)
% This function returns the second-order probability density function (PDF) 
% for the free surface elevation an irregular wave field, as described 
% and derived in Fuhrman et al. (2023).  Here p=p(zeta) is the PDF for a 
% given vector zeta=eta/sigma and scalar skewness Sk, where eta is the 
% free surface elevation and sigma is its standard deviation.  The 
% theoretical PDF, Eq. 2.32, is based on the Airy Ai function.  This is 
% replaced with a semi-theoretical negative tail to avoid spurious, non-
% physical oscillations for cases having Sk<=0.2, as in Eq. 5.3.  The 
% (optional) final argument modNegTail=0 will force use of the theoretical 
% PDF, whereas modNegTail=1 will force use of the modified negative tail.
%
% Example usage:
% >> zeta = [-8:0.01:8]; Sk = 0.3; p = pdfAiry(zeta,Sk); % Compute the PDF
% >> semilogy(zeta,p), xlim(8.*[-1 1]), ylim(10.^[-8 0]) % Plot the PDF
% >> xlabel('\zeta=\eta/\sigma'), ylabel('p(\zeta)=p(\eta)\sigma') % Label axes
%
% Programmed by: David R. Fuhrman, Tech. Univ. of Denmark, June 27, 2023
%
% Reference:
%
% Fuhrman, D.R., Klahn, M. & Zhai, Y. (2023) A new probability density 
% function for the surface elevation in irregular seas. J. Fluid Mech.
% 970, A38. https://doi.org/10.1017/jfm.2023.669

% Determine whether to modify the negative tail
if nargin < 3 % Only TRUE if modNegTail is not provided as input
    if Sk <= 0.2 % Recommended skewness threshold
        modNegTail = 1; % Modify negative negative tail
    else
        modNegTail = 0; % Do not modify negative tail
    end  
end

% Create the PDF
argExp = vpa(1/(3*Sk^2) + zeta/Sk); % Argument to exp, high precision
chi = vpa((2/Sk)^(1/3)*(1/(2*Sk)+zeta)); % Argument to Ai, high precision
p = (2/Sk)^(1/3)*exp(argExp).*airy(chi); % PDF, Eq. 2.31
if modNegTail % Modify left tail if desired
    chi0 = -1.17371; % Argument where Ai(chi0)=Ci(chi0) and Ai'(chi0)=Ci'(chi0)
    i0 = find(chi<chi0); % Indices where chi<chi0
    Ci = sqrt(airy(chi(i0)).^2+airy(2,chi(i0)).^2); % Airy envelope function
    p(i0) = (2/Sk)^(1/3)*exp(argExp(i0)).*Ci; % Modify negative tail, Eq. 5.3
end

end