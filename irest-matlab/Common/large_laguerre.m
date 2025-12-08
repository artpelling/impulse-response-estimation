function L = large_laguerre(a,n,z)
    L = zeros(length(z), n);
    f = sqrt(1 - abs(a)^2)./(1 - conj(a).*z);
    B = (z - a)./(1 - conj(a).*z);
    for k=1:n
        L(:,k) = f.';
        f = f.*B;
    end
end