function [co,err] = mlfdc_coeffs(signal, mpoles, eps)

% MLFDC_COEFFS - Calculates the mlfdc-coefficients of 'v' with respect to the
%                discrete modified basic rational system given by 'poles'.
%
% Usage: 
%     [co,err] = mlfdc_coeffs(signal,mpoles,eps)
%
% Input parameters:
%     v      : an arbitrary vector  
%     mpoles : poles of the modified basic rational system
%     eps    : accuracy of the discretization on the unit disc
%
% Output parameters:
%     co  : the Fourier coefficients of 'v' with respect to the discrete 
%           modified basic rational system defined by 'mpoles'
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
z=linspace(-pi,pi,m+1);
t=arg_inv(mpoles,z, eps);
samples = subsample(signal, t);
co=zeros(1,m);
bts = biortdc_system(mpoles, eps);

for i=1:1:m 
    co(i)=dotdc(samples,bts(i,:),mpoles,t);
end

len=length(signal);
mlf = mlf_system(len, mpoles);
err = norm(co*mlf - signal);