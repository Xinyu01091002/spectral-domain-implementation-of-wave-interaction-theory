% This is an an implementation of the third-order, directionally-spread,
% irregular wave theory of Madsen and Fuhrman (2012), with corrections as
% indicated in the appendix of Fuhrman et al. (2023).  Any publications
% in which results from this code are used are kindy asked to please cite 
% both of these works appropriately.
% 
% The code provided can compute both free surface elevations and kinematics
% for irregular waves to a desired order (1, 2 or 3).  The primary source
% code lives in the Source folder.
%
% In this repository, clearer source-level aliases have also been added in
% the Source folder:
%   mf12_direct_coefficients
%   mf12_spectral_coefficients
%   mf12_direct_surface
%   mf12_spectral_surface
% These aliases are intended to be easier to understand than some of the
% original development-era file names, while preserving the same behavior.
%
% Code was programmed by David R. Fuhrman, November 2022
%
% To view the data set and access the most recent version see:
% https://doi.org/10.11583/DTU.22060124
%
%
% References
%
% Fuhrman, D.R., Klahn, M. & Zhai, Y. (2023) A new probability density 
% function for the surface elevation in irregular seas. J. Fluid Mech. 970,
% A38. https://doi.org/10.1017/jfm.2023.669
%
% Madsen, P.A. & Fuhrman, D.R. (2012) Third-order theory for
% multi-directional irregular waves. J. Fluid Mech. 698, 304-334. 
% https://doi.org/10.1017/jfm.2012.87
