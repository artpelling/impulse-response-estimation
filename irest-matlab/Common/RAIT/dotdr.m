function s=dotdr(F, G, mpoles,t)

% DOTDR - Computes discrete real dot product of two function in H^2(ID).  
% 
% Usage:
%     s=dotdr(F,G,poles,t)
%
% Input parameters:
%     F,G : ID-->IR, anallytic functions on the unit disk
%     poles : poles of the rational system
%     t : argument(s)
%
% Output parameters:
%     s : values of the real dot product of 'F' and 'G' at 't' 
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

mpoles=[0 mpoles];
s=sum(F.*conj(G)./(2.*real(kernel(exp(1i*t(1:end)),exp(1i*t(1:end)),mpoles))-1));