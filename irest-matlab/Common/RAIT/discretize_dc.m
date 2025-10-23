function t = discretize_dc(mpoles, eps)

% DISCRETIZE_DC - Computes the non-equidistant complex discretization 
%                 on the unit disc that refers to the given poles.
%
% Usage: 
%     t = discretize_dc(mpoles,eps)
%   
% Input parameters:
%     mpoles : poles of the Blaschke product
%     eps    : accuracy of the complex discretization on the unit disc 
%
% Output parameters:
%     t : non-equidistant complex discretization
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

if max(abs(mpoles)) >= 1
    disp('Bad poles!');
    return;
end
if nargin < 3
    eps = 1e-6;
end

m = length(mpoles);
z = linspace(-pi,pi,m+1);
t = arg_inv(mpoles, z, eps);
