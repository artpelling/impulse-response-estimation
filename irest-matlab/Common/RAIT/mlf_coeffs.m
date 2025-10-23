function [co,err] = mlf_coeffs(v, poles)

% MLF_COEFFS - Calculates the mlf-coefficients of 'v' with respect to the
%              modified basic rational functions system given by 'poles'.
%
% Usage: 
%     [co,err] = mlf_coeffs(v,poles)
%
% Input parameters:
%     v     : an arbitrary vector  
%     poles : poles of the modified basic rational system
%
% Output parameters:
%     co  : the Fourier coefficients of v with respect to the modified 
%           basic rational system defined by 'poles'
%     err : L^2 norm of the approximation error 
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

np = size(poles,1);
[nv,mv] = size(v);
if np ~= 1
    disp('Wrong pole parameters!');
    return;
end
if nv ~= 1
    disp('Wrong vector to calc!');
    return;
end
if max(abs(poles)) >= 1
    disp('Bad poles!');
    return;
end

bts  = biort_system(mv,poles);
co = (bts * v' / mv)';
mlf  = mlf_system(mv,poles);
err = norm(co*mlf - v);
