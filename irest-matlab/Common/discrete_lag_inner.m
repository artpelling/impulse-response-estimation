function C = discrete_lag_inner(A, B, a, z)
    % Generate the discrete measure
    p = length(a);
    sigma = zeros(1, length(z));
    for k=1:p
        sigma = sigma + (1 - abs(a(k))^2)./(abs(1 - conj(a(k)).*z).^2);
    end

    S = diag(1./sigma);

    C = B'*(S*A);
end