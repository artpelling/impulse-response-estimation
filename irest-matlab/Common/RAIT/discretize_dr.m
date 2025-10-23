function t = discretize_dr(mpoles, eps)

% DISCRETIZE_DR - Computes the non-equidistant real discretization 
%                 on the unit disc that refers to the given poles.
%
% Usage: 
%     t = discretize_dr(mpoles,eps)
%   
% Input parameters:
%     mpoles : poles of the Blaschke product
%     eps    : accuracy of the real discretization on the unit disc 
%
% Output parameters:
%     t : non-equidistant real discretization
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

if max(abs(mpoles)) >= 1
    disp('Bad poles!');
    return;
end
if nargin < 3
    eps = 1e-6;
end

mpoles = [0 mpoles];
m=length(mpoles);
z=-(m-1)*pi:pi:(m-1)*pi;
z=z./m;
t=argdr_inv(mpoles,z,eps);
