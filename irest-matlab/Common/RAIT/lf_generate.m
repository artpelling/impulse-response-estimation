function v = lf_generate(len,poles,coeffs)

% LF_GENERATE - Generates a function in the space spanned by the LF system.
%
% Usage: 
%     v = lf_generate(len,poles,coeffs)
%
% Input parameters:
%     len    : number of points in case of uniform sampling 
%     poles  : poles of the rational system (row vector)
%     coeffs : coefficients of the linear combination to form (row vector)
%
% Output parameters:
%     v : the generated function at the uniform sampling points as a
%         row vector
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??


[np,mp] = size(poles);
[nl,ml] = size(coeffs);
if np ~= 1 || nl ~= 1 || mp ~= ml || len < 2
    error('Wrong parameters!');
end
if max(abs(poles)) >= 1
    error('Bad poles!');
end

% The v is the linear combination of the LF system elements.
%
%          lf1  lf1  lf1 ... lf1
%          lf2  lf2  lf2 ... lf2
% co1 co2  v    v    v       v

v = coeffs * lf_system(len, poles);
