function Q = z_trans(q, w)
%ZTRANS Compute the Z transform of the sequence q, over the Torus points
%specified in w.
    N = length(w);
    n = 0:length(q)-1;
    Q = zeros(N, 1);
    for k=1:N
        z = w(k).^n;
        Q(k) = sum(q.*z);
    end
end