function mts = mtdc_system(mpoles, eps)

% MTDC_SYSTEM - Generates the discrete complex MT system.
%
% Usage: 
%     mts = mtdc_system(mpoles,eps)
%
% Input parameters:
%     mpoles : poles of the discrete complex MT system
%     eps    : accuracy of the discretization on the unit disc 
%
% Output parameters:
%     mts : the elements of the discrete complex MT system 
%           at the uniform sampling points as row vectors
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??.

if max(abs(mpoles)) >= 1
    disp('Bad poles!');
    return;
end
if nargin < 3
    eps = 1e-6;
end

m=length(mpoles);
t = discretize_dc(mpoles, eps);
mts = zeros(m,m+1); 

for j = 1:m
    mts(j,:)=MT(j-1,mpoles,exp(1i.*t));
end

% -------------------------------------------------------------------------

function r=MT(n, mpoles, z)

% Compute the values of the nth Malmquist-Takenaka function at z.

r=1;
for k=1:1:n
   r=r.* (z-mpoles(k)) ./ (1-conj(mpoles(k)).*z);
end
r=r.* sqrt(1-abs(mpoles(n+1)).^2) ./ (1-conj(mpoles(n+1)).*z);