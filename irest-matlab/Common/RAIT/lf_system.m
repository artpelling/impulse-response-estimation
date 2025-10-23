function lfs = lf_system(len, poles)

% LF_SYSTEM - Generates the linearly independent system defined by poles.
%
% Usage: 
%     lfs = lf_system_(len,poles)
%
% Input parameters:
%     len   : number of points in case of uniform sampling 
%     poles : poles of the rational system
%
% Output parameters:
%     lfs : the elements of the linearly independent system at the uniform
%           sampling points as row vectors
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??


[np,mp] = size(poles);
if np ~= 1 || len < 2
    error('Wrong parameters!');
end
if max(abs(poles)) >= 1
    error('Bad poles!');
end

lfs = zeros(mp,len);
t = linspace(-pi,pi,len+1);
t = t(1:len);
z = exp(1i*t);

for j = 1:mp
    rec = 1 ./ (1 - conj(poles(j)).*z);
    lfs(j,:) = rec .^ multiplicity_local(j,poles);
    lfs(j,:) = lfs(j,:) ./ sqrt(lfs(j,:) * lfs(j,:)' / len);
end


% -------------------------------------------------------------------------

function m = multiplicity_local(n,v)

% Returns the multiplicity of the nth element of the vector v.

m = 0;
for k = 1:n
    if v(k) == v(n)
        m = m + 1;
    end
end
