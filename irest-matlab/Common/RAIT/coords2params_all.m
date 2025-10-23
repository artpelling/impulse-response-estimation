function p = coords2params_all(k)

% COORDS2PARAMS_ALL - Maps coordinates in IR^2 to parameters in ID. 
%
% Usage: 
%     p = coords2params_all(k)
% 
% Input parameters:
%     k : matrix of coordinate pairs in IR^2, rows are considered as
%         vertices of the simplex
%
% Output parameters:
%     p : row vector of corresponding parameters in ID
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??


vertnum = size(k,1);
parnum = size(k,2)/2;
p = zeros(vertnum,parnum);

for i = 1:vertnum
    p(i,:) = coords2params(k(i,:));
end