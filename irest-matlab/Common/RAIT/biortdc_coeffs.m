function [co,err] = biortdc_coeffs(v, mpoles, eps)

% BIORTDC_COEFFS - Calculates the discrete biorthogonal-coefficients of 'v'
%                  with respect to the biorthogonal system given by 'poles'.
%
% Usage: 
%     [co,err] = biortdc_coeffs(v,mpoles,eps)
%
% Input parameters:
%     v      : an arbitrary vector  
%     mpoles : poles of the biorthogonal system
%     eps    : accuracy of the discretization on the unit disc
%
% Output parameters:
%     co  : the Fourier coefficients of v with respect to the discrete 
%           biorthogonal system defined by 'mpoles'
%     err : L^2 norm of the approximation error 
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??


[np,mp] = size(mpoles);
[nv,mv] = size(v);
if np ~= 1
    error('Wrong pole parameters!');
    return;
end
if nv ~= 1
    error('Wrong vector to calc!');
    return;
end
if max(abs(mpoles)) >= 1
    disp('Bad poles!');
    return;
end
if nargin<3
    eps=1e-6;
end

m=length(mpoles);
t = discretize_dc(mpoles, eps);
samples = subsample(v, t);
co=zeros(1,m);
mlf = mlfdc_system(mpoles, eps);

for i=1:1:m 
    co(i)=dotdc(samples(1:end),mlf(i,1:end),mpoles,t);
end

len=length(v);
bts = biort_system(len, mpoles);
err = norm(co*bts - v);
