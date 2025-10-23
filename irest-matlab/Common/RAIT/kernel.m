function r=kernel(y, z, mpoles)

% KERNEL - Computes the weight function of discrete dot product in H^2(D).  
% 
% Usage:
%     r=kernel(y,z,poles)
%
% Input parameters:
%     y,z   : arguments
%     poles : poles of the rational system
%
% Output parameters:
%     r : values of the weight function at arguments 'y' and 'z'
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??

n=1;
r=0;
m=length(mpoles);
if y==z
   for k=1:1:m
       alpha=angle(mpoles(k));
       R=abs(mpoles(k));
       t=angle(z);
       r=r+poisson(R,t-alpha);
   end
else
    for i=n:1:m
        r=r+MT(i-1,mpoles,y).*conj(MT(i-1,mpoles,z));   
    end
end

% -------------------------------------------------------------------------


function z=poisson(r, t)

% Compute the values of the poisson function at (r,t).

z=(1-r.^2)./(1-2.*r.*cos(t)+r.^2);

% -------------------------------------------------------------------------

function r=MT(n, mpoles, z)

% Compute the values of the nth Malmquist-Takenaka function at z.

r=1;
for k=1:1:n
   r=r.* (z-mpoles(k)) ./ (1-conj(mpoles(k)).*z);
end
r=r.* sqrt(1-abs(mpoles(n+1)).^2) ./ (1-conj(mpoles(n+1)).*z);