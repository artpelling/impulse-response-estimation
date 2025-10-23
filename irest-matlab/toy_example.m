%% Clear & close all
clear all;
close all;

%% Import libs for rational functions
addpath("Common/RAIT/");

%% Constants
M = 5000; % Initial number of (equidistant) sampling points on the boundary
N = 512; % final index of impulse response -> even for 512, we get numerical errors if we are not careful with computation
n = 0:N-1;
a = 0.4 + 1i*0.2; % Blaschke parameter (this determines the Laguerre system)

%% Callbacks
B_inv_a = @(z) (z + a)./(1 + conj(a).*z); % For inverting the Blaschke-transformation

%% Define example system - 3 poles
alpha = [ 0.7 - 1i*0.5; 0.7 + 1i*0.5; 0.9]; % Mirror image poles
c = [1 + 1i*2; 1 - 1i*2; 1]; % Residues

% Compute actual impulse response up to index N
h_true = sum( c.*conj(alpha).^(n) );

%% Compute discretization set & generate Laguerre functions

% Initial discretization set -> this is a first example, so let's leave it
% with equidistant sampling right now
f = linspace(-pi, pi, M+1); f(end) = [];
z = exp(1i*f); % Sampling points on the Torus

%% Generate system input and output

% Generate input: geometric sequence -> we can compute the Z transform exactly
% on the grid by the frequencies f. We can also do this with other input
% choices, but I wanted to be sure that I did everything correctly here.
r = 0.9;
U = fftshift(ifft(r.^(0:M-1)))*M;

% OPTIONAL: We can check the input has been generated correctly witht the
% below formulas
% U_hat = 1./(1 - r.*z); 
% figure; plot(abs(U-U_hat))

% Generate output

% First generate the true transfer function sampled over z
H = zeros(size(z));
for k=1:length(c)
    H = H + c(k)./(1 - conj(alpha(k)).*z);
end

% Now generate true output
Y = U.*H;

%% Solve for discrete Laguerre coefficients

% Again, this is a toy example, this should not be done like this, but I
% want to show here that the pipeline works
b = periodize_poles(a, N);
c_lag = mt_coeffs(Y./U, b);

%% Use sequence of Laguerre coeffs as an impulse response and go to the Z domain
angs = linspace(-pi, pi, N+1); angs(end) = [];
zz = exp(1i*angs);
H_hat = fftshift(ifft(c_lag))*N;

%% Recover the impulse response
h = zeros(N, 1);
phi1 = (1 - conj(a).* B_inv_a(zz) )/sqrt(1 - abs(a)^2);
phi2 = B_inv_a(zz);
for k=1:N
    phik = phi1.*phi2.^(k-1);
    h(k) = (H_hat * phik')/N;
end

%% Plot results
figure;
subplot(211);
plot(h_true, "LineWidth", 2); hold on;
plot(real(h), "r");
axis tight;
grid on;
xlabel("Discrete timestamps");
ylabel("Impulse response values");
legend("True IR", "Recovered IR");
title("Impulse response comparison");
subplot(212);
plot(abs(h - h_true.'), "LineWidth", 2);
axis tight;
grid on;
xlabel("Discrete timestamps");
ylabel("$| h - \hat{h}|$", 'Interpreter', 'latex');
title("Absolute error");

saveas(gcf, 'figures/toy.png');
