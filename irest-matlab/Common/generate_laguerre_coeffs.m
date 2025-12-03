function cfs = generate_laguerre_coeffs(c, alpha, a, N, disc)
   M = length(alpha);
   cfs = zeros(1, N);
   for k=1:N
       for j=1:M
            s = conj(c(j)).*sqrt(1 - abs(a)^2)./(1 - a*conj(alpha(j))) * ...
                ( (conj(alpha(j)) - conj(a))/(1 - a*conj(alpha(j))))^(k-1); 
            if disc
                s = s/( 1 - (( (conj(alpha(j)) - conj(a))/(1 - a*conj(alpha(j)))))^N );
            end
            cfs(k) = cfs(k) + s;
       end
   end
end