function y = causal_conv(h, u)
%CAUSAL_CONV causal convolution to compute the output of causual SISO-LTI
%systems given the impulse response and the input.
    N = length(u);
    M = length(h);
    y = zeros(1, N);

    for n=1:N
        for k=1:n
            y(n) = y(n) + h(k)*u(n-k+1);
        end
    end
end