function bd = arg_der(a,t)

% ARG_DER - Derivatives of the argument function of a Blaschke product.
%
% Usage: 
%     bd = arg_der(a,t)
%
% Input parameters:
%     a : parameters of the Blaschke product
%     t : values in [-pi,pi), where the function values are needed
%
% Output parameters:
%     bd : the derivatives of the argument function at the points in t
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??


na = size(a,1);
nt = size(t,1);
if na ~= 1 || nt ~= 1
    error('Parameters should be row vectors!');
end
if max(abs(a)) >= 1
    error('Bad poles!');
end


bd = zeros(1,length(t));
for i = 1:length(a);
    bd = bd + arg_der_one(a(i),t);
end
bd = bd / length(a);


% -------------------------------------------------------------------------

function bd = arg_der_one(a,t)

r = abs(a);
fi = angle(a);

bd = (1 - r^2)./(1 + r^2 - 2*r*cos(t-fi));
