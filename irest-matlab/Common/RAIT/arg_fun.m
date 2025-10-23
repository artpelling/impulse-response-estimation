function b = arg_fun(a,t)

% ARG_FUN - Values of the argument function of a Blaschke product.
%
% Usage: 
%     b = arg_fun(a,t)
%
% Input parameters:
%     a : parameters of the Blaschke product
%     t : values in [-pi,pi), where the function values are needed
%
% Output parameters:
%     b : the values of the argument function at the points in t
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


b = zeros(1,length(t));
for i = 1:length(a);
    b = b + arg_fun_one(a(i),t);
end
b = b / length(a);


% -------------------------------------------------------------------------

function b = arg_fun_one(a,t)

r = abs(a);
fi = angle(a);
mu = (1+r)/(1-r);

gamma = 2 * atan((1/mu)*tan(fi/2));

b = 2*atan(mu*tan((t-fi)/2)) + gamma;
b = mod(b+pi,2*pi)-pi; % move it in [-pi,pi)
