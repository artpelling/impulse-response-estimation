function b = blaschkes(len,poles)

% BLASCHKES - Gives a sampled Blaschke function.
%
% Usage: 
%     b = blaschkes(len,poles)
%
% Input parameters:
%     len   : number of points in case of uniform sampling 
%     poles : parameters of the Blaschke product
%
% Output parameters:
%     b : the Blaschke product with given parameters sampled at uniform
%         points on the torus
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

t = linspace(-pi,pi,len+1);
z = exp(1i*t);
b = ones(1,len+1);
for p = 1:length(poles)
    b = b .* (z - poles(p))./(1 - conj(poles(p)).*z);
end

