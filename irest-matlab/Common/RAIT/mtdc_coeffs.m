function [co, err]=mtdc_coeffs(signal, mpoles, eps)

% MTDC_COEFFS - Calculates the mtdc-coefficients of 'v' with respect to the
%               discrete complex MT system given by 'mpoles'.
%
% Usage: 
%     [co,err] = mtdc_coeffs(signal,mpoles,eps)
%
% Input parameters:
%     v      : an arbitrary vector  
%     mpoles : poles of the discrete complex MT system
%     eps    : accuracy of the complex discretization on the unit disc
%
% Output parameters:
%     co  : the Fourier coefficients of 'v' with respect to the discrete 
%           complex MT system defined by 'mpoles'
%     err : L^2 norm of the approximation error 
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

np = size(mpoles,1);
nv = size(signal,1);
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

m=length(mpoles);
t = discretize_dc(mpoles, eps);
samples = subsample(signal, t);
co=zeros(1,m);
mts = mtdc_system(mpoles, eps);

for i=1:1:m 
    co(i)=dotdc(samples,mts(i,:),mpoles,t);
end

len=length(signal);
mts = mt_system(len, mpoles);
err = norm(co*mts - signal);
