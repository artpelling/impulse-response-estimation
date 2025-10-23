function v = mlfdc_generate(len, mpoles, coeffs)

% MLFDC_GENERATE - Generates a function in the space spanned by 
%                  the discrete modified basic rational function system.
%
% Usage: 
%     v = mlfdc_generate(len,mpoles,coeffs)
%
% Input parameters:
%     len    : number of points in case of uniform sampling 
%     mpoles : poles of the discrete modified basic rational system (row vector)
%     coeffs : coefficients of the linear combination to form (row vector)
%
% Output parameters:
%     v : the generated function at the uniform sampling points as a
%         row vector
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

[np,mp] = size(mpoles);
[nu,mu] = size(coeffs);

if np ~= 1 || nu ~= 1 || mp ~= mu || len < 2
    error('Wrong parameters!');
end
if max(abs(mpoles)) >= 1
    error('Bad poles!');
end

% The v is the linear combination of the discrete modified rational system elements.
%
%          mlfdc1  mlfdc1  mlfdc1 ... mlfdc1
%          mlfdc2  mlfdc2  mlfdc2 ... mlfdc2
% co1 co2  v       v       v          v

mlf = mlf_system(len, mpoles);
v=coeffs*mlf;
