function R = generate_rational(c, alpha, z)
    R = zeros(1, length(z));
    M = length(c);
    for k=1:M
        R = R + c(k)*(1 ./ ( 1 - conj(alpha(k)).*z ));
    end
end