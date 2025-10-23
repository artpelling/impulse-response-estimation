function [cUk,cVk,err]=mtdr_coeffs(v, mpoles, eps)

% MTDC_COEFFS - Calculates the mtdr-coefficients of 'v' with respect 
%               to the discrete real MT system given by 'mpoles'.
% 
% Usage:
%   [cUk,cVk,err]=mtdr_coeffs(v,mpoles,eps)
%
% Input Parameters:
%     v      : an arbitrary vector
%     mpoles : poles of the discrete real MT system
%     eps    : accuracy of the real discretization on the unit disc
% 
% Output parameters:
%     cUk : the Fourier coefficients of 'v' with respect to the real part 
%           of the discrete real MT system defined by 'mpoles'
%     cVk : the Fourier coefficients of 'v' with respect to the imaginary 
%           part of the discrete real MT system defined by 'mpoles'
%     err : L^2 norm of the approximation error 
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

np = size(mpoles,1);
nv = size(v,1);
if np ~= 1
    error('Wrong pole parameters!');
end
if nv ~= 1
    error('Wrong vector to calc!');
end
if max(abs(mpoles)) >= 1
    disp('Bad poles!');
    return;
end
if nargin<3
    eps=1e-6;
end

m=length(mpoles)+1;
t = discretize_dr(mpoles, eps);
samples = subsample(v, t);
cUk=zeros(1,m);
cVk=zeros(1,m);
[mts_re, mts_im] = mtdr_system(mpoles, eps);

for i=1:1:m
    cUk(i)=dotdr(samples,mts_re(i,:),mpoles,t);
    cVk(i)=dotdr(samples,mts_im(i,:),mpoles,t); 
end

SRf = mtdr_generate(length(v), mpoles, cUk, cVk);
err = norm(SRf-v);
