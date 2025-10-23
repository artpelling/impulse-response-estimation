function v = biortdc_generate(len, mpoles, coeffs)

% BIORTDC_GENERATE - Generates a function in the space spanned by 
%                    the discrete biorthogonal system.
%
% Usage: 
%     v = biortdc_generate(len,mpoles,coeffs)
%
% Input parameters:
%     len    : number of points in case of uniform sampling 
%     mpoles : poles of the biorthogonal system (row vector)
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

% The v is the linear combination of the discrete biorthogonal system elements.
%
%          btdc1  btdc1  btdc1 ... btdc1
%          btdc2  btdc2  btdc2 ... btdc2
% co1 co2  v      v      v         v

bts = biort_system(len, mpoles);
v=coeffs*bts;
