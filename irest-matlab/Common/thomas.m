function x = thomas(a, b, c, d)
% THOMAS_ALGORITHM  Solve a tridiagonal system Ax = d using the Thomas method.
%   a, b, c : sub-, main-, and super-diagonals (column vectors)
%   d       : right-hand side
% All vectors must be of length n.

    n = length(b);

    % Forward sweep
    c_star = zeros(n,1);
    d_star = zeros(n,1);

    % First row
    denom = b(1);
    c_star(1) = c(1) / denom;
    d_star(1) = d(1) / denom;

    % Eliminate forward
    for i = 2:n
        denom = b(i) - a(i) * c_star(i-1);
        c_star(i) = c(i) / denom;
        d_star(i) = (d(i) - a(i) * d_star(i-1)) / denom;
    end

    % Back substitution
    x = zeros(n,1);
    x(n) = d_star(n);
    for i = n-1:-1:1
        x(i) = d_star(i) - c_star(i) * x(i+1);
    end
end