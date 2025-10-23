function bo = bisection_order(n)

% BISECTION_ORDER - Gives a better order for multiple bisection runs.
%
% Usage: 
%     bo = bisection_order(n)
%
% Input parameters:
%     n : number of points
%
% Output parameters:
%     bo : a 3-by-(n+1) matrix with the order of calculation
%
% Copyright: (C) ELTE IK NumAnal, GPL 1.1 ??


% initialization

bo = zeros(n+1,3); 
bo(1,:) = [0,-1,-1];
bo(2,:) = [n,-1,-1];
bo(3,:) = [floor((n-0)/2),0,n];

watch = 3; % which column is currently watched
fill = 4; % where to fill the new values

% fill the matrix with the ordering

while fill <= n+1
    % names
    ch = bo(watch,1); % child
    p1 = bo(watch,2); % parents
    p2 = bo(watch,3); % INVAR: p1 < ch < p2
    
    % if there is place for a nother element...
    % the child with parent 1
    if ch-p1 > 1 && fill <= n+1
        gch = floor((ch+p1)/2); % 'grand child'
        bo(fill,:) = [gch,p1,ch];
        fill = fill + 1;
    end
    % the child with parent 2
    if p2-ch > 1 && fill <= n+1
        gch = floor((ch+p2)/2); % 'grand child'
        bo(fill,:) = [gch,ch,p2];
        fill = fill + 1;
    end
    watch = watch + 1;
end
