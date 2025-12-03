function A = construct_eq_for_disc_lag(ul, a)
    N = length(ul);
    A = zeros(N);
    for k=1:N
        c1 = [zeros(1, k-1), ul(1:N-k+1)];
        c2 = [zeros(1,k), conj(a)*ul(1:N-k)];
        A(:,k) = c1+c2;
    end
end