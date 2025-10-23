function s=dotdc(F,G,poles,t)

% DOTDC - Computes complex discrete dot product of two function in H^2(ID).  
% 
% Usage:
%     s=dotdc(F,G,poles,t)
%
% Input parameters:
%     F,G : ID-->IC, anallytic functions on the unit disk
%     poles : poles of the rational system
%     t : argument(s)
%
% Output parameters:
%     s : values of the complex dot product of 'F' and 'G' at 't' 
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

s=sum(F(1:end-1).*conj(G(1:end-1))./kernel(exp(1i*t(1:end-1)),exp(1i*t(1:end-1)),poles));

