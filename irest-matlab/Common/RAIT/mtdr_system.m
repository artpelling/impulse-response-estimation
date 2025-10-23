function [mts_re, mts_im] = mtdr_system(poles, eps)

% MTDR_SYSTEM - Generates the discrete real MT system.
%
% Usage: 
%     [mts_re,mts_im] = mtdr_system(poles,eps)
%
% Input parameters:
%     poles : poles of the discrete real MT system
%     eps   : accuracy of the real discretization on the unit disc 
%
% Output parameters:
%     mts_re : the real part of the discrete complex MT system 
%              at the non-equidistant discretization on the unit disc
%     mts_im : the imaginary part of the discrete complex MT system 
%              at the non-equidistant discretization on the unit disc
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

if max(abs(poles)) >= 1
    disp('Bad poles!');
    return;
end
if nargin < 3
    eps = 1e-6;
end

mpoles=[0 poles];
m=length(mpoles);
t = discretize_dr(poles, eps);
mts_re = zeros(m,length(t)); 
mts_im = zeros(m,length(t));

for j = 1:m
    mts_re(j,:)=real(MT(j-1,mpoles,exp(1i.*t)));
    mts_im(j,:)=imag(MT(j-1,mpoles,exp(1i.*t)));    
end

% -------------------------------------------------------------------------

function r=MT(n, mpoles, z)

% Compute the values of the nth Malmquist-Takenaka function at z.

r=1;
for k=1:1:n
   r=r.* (z-mpoles(k)) ./ (1-conj(mpoles(k)).*z);
end
r=r.* sqrt(1-abs(mpoles(n+1)).^2) ./ (1-conj(mpoles(n+1)).*z);