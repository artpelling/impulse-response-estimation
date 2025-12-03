%% Clear & Close all
clear all;
close all;

%% Add needed libs
addpath("Common/RAIT/");

%% Load data
data = load("../benchmarks/data/fabian.mat");

%% Save input, output and true impulse response
h_true = data.h;
u = data.u;
y = data.y;

%% Constants
M = length(u); % Initial number of (equidistant) sampling points on the boundary
N = length(h_true); % final index of impulse response -> even for 512, we get numerical errors if we are not careful with computation
n = 0:N-1;
a = 0.2; % Blaschke parameter (this determines the Laguerre system)

%% Callbacks
B_inv_a = @(z) (z + a)./(1 + conj(a).*z); % For inverting the Blaschke-transformation

%% Generate system input and output in the Frequency domain

U = fftshift(ifft(u))*M;
Y = fftshift(ifft(y))*M;

%% Solve for discrete Laguerre coefficients

% Again, this is a toy example, this should not be done like this, but I
% want to show here that the pipeline works
b = periodize_poles(a, N);
L = mt_system(M, b).';
c_lag = ((diag(U)*L)\Y.').';

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