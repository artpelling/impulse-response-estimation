function p = coords2params(k)

% COORDS2PARAMS - Maps coordinates in IR^2 to parameters in ID. One row.
%
% Usage: 
%     p = coords2params(k)
%
% Input parameters:
%     k : row vector of coordinate pairs in IR^2
%
% Output parameters:
%     p : row vector of corresponding parameters in ID
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??


parnum = size(k,2)/2;
p = zeros(1,parnum);

for j = 1:parnum    
	u = k(1,2*j-1);
	v = k(1,2*j);
    r = sqrt(u^2 + v^2);
	x = u / sqrt(1 + r^2);
	y = v / sqrt(1 + r^2);
	z = x + 1i*y;
    p(1,j) = z;
end
