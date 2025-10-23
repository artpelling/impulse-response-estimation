function co=coeffd_conv(poles, coeffs, base1, base2, eps)

% COEFF_CONV - Converts the coefficients coeffs between the discrete 
%              systems base1 and base2.
%
% Usage: 
%     co=coeff_conv(poles,coeffs,base1,base2,eps)
%
% Input parameters:
%     poles  : poles of the discrete systems
%     coeffs : coefficients with respect to the discrete system 'base1' 
%     base1  : type of the discrete system to be converted
%     base2  : type of the converted discrete system  
%     eps    : accuracy of the discretization on the unit disc 
%
% Output parameters:
%     co : converted coefficients with respect to the system 'base2' 
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

if max(abs(poles)) >= 1
    error('Bad poles!');
end
if size(coeffs,1)~=1
    error('Coeffs should be row vector!');
end    
if nargin < 5
    eps = 1e-6;
end

switch base1 
    case 'mlfdc'
        g1 = mlfdc_system(poles,eps);
    case 'biortdc'
        g1 = biortdc_system(poles,eps);
    case 'mtdc'
        g1 = mtdc_system(poles,eps);
    otherwise
        error('Bad conversation!');
end

switch base2 
    case 'mlfdc'
        g2 = mlfdc_system(poles,eps);
    case 'biortdc'
        g2 = biortdc_system(poles,eps);
    case 'mtdc'
        g2 = mtdc_system(poles,eps);
    otherwise
        error('Bad conversation!');
end

F = g1 * g1' / length(coeffs);
G = g1 * g2' / length(coeffs);
co = G \ F * coeffs';
co=co';