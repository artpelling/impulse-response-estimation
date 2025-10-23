function [co,err] = biort_coeffs(v,poles)

% BIORT_COEFFS - Calculates the biorthogonal-coefficients of 'v' with
%                respect to the biorthogonal system given by 'poles'.
%
% Usage: 
%     [co,err] = biort_coeffs(v,poles)
%
% Input parameters:
%     v     : an arbitrary vector  
%     poles : poles of the biorthogonal system
%
% Output parameters:
%     co  : the Fourier coefficients of v with respect to the biorthogonal 
%           system defined by poles
%     err : L^2 norm of the approximation error 
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

[np,mp] = size(poles);
[nv,mv] = size(v);
if np ~= 1
    error('Wrong pole parameters!');
    return;
end
if nv ~= 1
    error('Wrong vector to calc!');
    return;
end
if max(abs(poles)) >= 1
    disp('Bad poles!');
    return;
end

mlfs = mlf_system(mv,poles);
bts = biort_system(mv,poles);
co = (mlfs * v' / mv)';
err = norm(co*bts - v);
