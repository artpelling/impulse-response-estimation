function SRf=mtdr_generate(len, mpoles, cUk, cVk)

% MTDR_GENERATE - Generates a function in the space spanned by 
%                 the discrete real MT system.
%
% Usage: 
%     v = mtdr_generate(len,mpoles,cUk,cVk)
%
% Input parameters:
%     len    : number of points in case of uniform sampling 
%     mpoles : poles of the discrete real MT system (row vector)
%     cUk    : coefficients of the linear combination to form (row vector)
%              with respect to the real part of the discrete real MT system
%              defined by 'mpoles'
%     cVk    : coefficients of the linear combination to form (row vector)
%              with respect to the imaginary part of the discrete real MT 
%              system defined by 'mpoles'
%
% Output parameters:
%     SRf : the generated function at the uniform sampling points as a
%           row vector
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

mpoles=[0 mpoles];
[np,mp] = size(mpoles);
[nu,mu] = size(cUk);
[nv,mv] = size(cVk);

if np ~= 1 || nu ~= 1 || mp ~= mu || nv ~= 1 || mp ~= mv || len < 2
    error('Wrong parameters!');
end
if max(abs(mpoles)) >= 1
    error('Bad poles!');
end

% The v is the linear combination of the discrete real MT system elements.
%
%          mtdr1  mtdr1  mtdr1 ... mtdr1
%          mtdr2  mtdr2  mtdr2 ... mtdr2
% co1 co2  v      v      v         v

mts = mt_system(len, mpoles);
m=length(mpoles);
SRf=cUk(1).*mts(1,:);
for i=2:1:m
     SRf=SRf+2*cUk(i).*real(mts(i,:))+2*cVk(i).*imag(mts(i,:));
end