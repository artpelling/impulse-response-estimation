function y = subsample(sample,x)

% SUBSAMPLE - Interpolates values between uniform sampling points.
%
% Usage: 
%     y = subsample(sample,x)
%
% Input parameters:
%     sample : row vector of uniformly sampled values on [-pi,pi) 
%     x      : the value at x is to be found
%
% Output parameters:
%     y : the interpolated value at the point x according to sample
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

len = size(sample,2);
sx = linspace(-pi,pi,len);

y = interp1(sx,sample,x,'linear');