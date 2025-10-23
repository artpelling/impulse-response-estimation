function co=coeff_conv(len, poles, coeffs, base1, base2)

% COEFF_CONV - Converts the coefficients coeffs between the continuous 
%              systems base1 and base2.
%
% Usage: 
%     co=coeff_conv(len,poles,coeffs,base1,base2)
%
% Input parameters:
%     len    : number of points in case of uniform sampling 
%     poles  : poles of the continuous systems
%     coeffs : coefficients with respect to the continuous system 'base1' 
%     base1  : type of the continuous system to be converted
%     base2  : type of the converted continuous system  
%
% Output parameters:
%     co : converted coefficients with respect to the system 'base2' 
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

np = size(poles,1);
if np ~= 1 || len < 2
    error('Wrong parameters!');
end
if max(abs(poles)) >= 1
    error('Bad poles!');
end
if size(coeffs,1)~=1
    error('Coeffs should be row vector!');
end    

switch base1 
    case 'lf'
        g1 = lf_system(len,poles);
    case 'mlf'
        g1 = mlf_system(len,poles);
    case 'biort'
        g1 = biort_system(len,poles);
    case 'mt'
        g1 = mt_system(len,poles);
    otherwise
        error('Bad conversation!');
end

switch base2 
    case 'lf'
        g2 = lf_system(len,poles);
    case 'mlf'
        g2 = mlf_system(len,poles);
    case 'biort'
        g2 = biort_system(len,poles);
    case 'mt'
        g2 = mt_system(len,poles);
    otherwise
        error('Bad conversation!');
end

F = g1 * g1' / len;
G = g1 * g2' / len;
co = G \ F * coeffs';
co=co';