function v = biort_generate(len, poles, coeffs)

% BIORT_GENERATE - Generates a function in the space spanned by 
%                  the biorthogonal system.
%
% Usage: 
%     v = biort_generate(len,poles,coeffs)
%
% Input parameters:
%     len    : number of points in case of uniform sampling 
%     poles  : poles of the biorthogonal system (row vector)
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
    disp('Wrong parameters!');
    return;
end
if max(abs(poles)) >= 1
    disp('Bad poles!');
    return;
end

% The v is the linear combination of the biorthogonal system elements.
%
%          bt1  bt1  bt1 ... bt1
%          bt2  bt2  bt2 ... bt2
% co1 co2  v    v    v       v

v = coeffs * biort_system(len, poles);