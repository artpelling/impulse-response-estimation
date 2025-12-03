%% Clear & close all
clear all;
close all;

%% Import libs for rational functions
addpath("Common");
addpath("Common/RAIT");

%% For reproducibility
rng('default');

%% Constants
N = 20000; % final index of impulse response -> even for 512, we get numerical errors if we are not careful with computation
M = 10000; % Number of poles that define the system -> even number to generate poles simply
n = 0:N-1;
a = 0.1 + 1i*0.1; % Blaschke parameter (this determines the Laguerre system)

%% Callbacks
B_inv_a = @(z) (z + a)./(1 + conj(a).*z); % For inverting the Blaschke-transformation
B_a = @(z) (z-a)./(1 - conj(a).*z);

%% Define example system - 5000 poles
alpha_r = rand(M/2, 1); % abs for mirror image poles
c_r = rand(M/2, 1); % abs for residues

alpha_fi = -pi + 2*pi*rand(M/2,1);
c_fi = -pi + 2*pi*rand(M/2,1);

alpha = alpha_r.*exp(1i*alpha_fi);
c = c_r.*exp(1i*c_fi);

alpha = [alpha; conj(alpha)];
c = [c; conj(c)];

% Compute actual impulse response up to index N
h_true = sum( c.*conj(alpha).^(n) );

%% Some more constants (discretization sets)
% Initial discretization set -> this is a first example, so let's leave it
% with equidistant sampling right now
f = linspace(-pi, pi, N+1); f(end) = [];
z = exp(1i*f); % Equidistant sampling points on the Torus

%% Plot poles and true impulse response
figure;
subplot(121);
plot(real(h_true)); 
grid on;
title("True impulse response");
subplot(122);
plot(real(z), imag(z), 'k', 'LineWidth', 2);
hold on;
plot(real(alpha), imag(alpha), 'ro');
title("Location of mirror image poles");
grid on;
axis square;

%% Generate system input and output

% Generate input -> we use the zeroth Laguerre function to obtain a
% tridiagonal system for the Laguerre coefficients of the transfer
% function. In the time domain this means u_k = C * sqrt(1 - |a|^2) *
% conj(a)^k, where C is chosen arbitrarily. Let's discuss if this kind of
% input is feasible for the acoustic application. Other types of inputs can
% also be used, but then the SLE later on won't be as pretty.
c0 = 1000;
u = c0*sqrt(1 - abs(a)^2)*conj(a).^(0:N-1);
y = causal_conv(h_true, u);

%% Solve for discrete Laguerre coefficients

% Processing pipeline starts here, start recording
tic;

% First compute discretization set on T, for which the Laguerre matrix is
% unitary.

w = (z + a)./(1 + conj(a).*z);

% Generate Z transforms of input and output over w
U = z_trans(u, w);
Y = z_trans(y, w);

% Generate the big matrix of Laguerre functions
L = large_laguerre(a, N, w);

% First we compute the Laguerre coefficients (up to N) of U and Y. This we
% can do in a numerically safe manner.
aa = ones(1, N).*a; % vector of MT parameters
clu = discrete_lag_inner(U, L, aa, w).';
cly = discrete_lag_inner(Y, L, aa, w);

% Now we solve a tridiagonal system using the Thomas algorithm to obtain
% the Laguerre coefficients of the transfer function
A = construct_eq_for_disc_lag(clu, a);

c_lag = A\(cly*sqrt(1 - abs(a)^2));

%% Check if coefficients can be used to reconstruct H_true at the sampling points (it should interpolate)
H_true = generate_rational(c, alpha, w);
H_est = L*c_lag;

figure;
subplot(211);
plot(real(H_true), 'LineWidth', 2); hold on;
plot(real(H_est));
grid on;
legend("True transfer function", "Discrete orthonormal Laguerre reconstruction");
title("Real part");
subplot(212);
plot(imag(H_true), 'LineWidth', 2); hold on;
plot(imag(H_est));
legend("True transfer function", "Discrete orthonormal Laguerre reconstruction");
title("Imaginary part");
grid on;

%% Use sequence of Laguerre coeffs as an impulse response and go to the Z domain 

% Note: from here on out, I use equidistant sampling on the Torus. I'm sure
% this can be improved (i.e., find a new non-Laguerre discretization
% scheme and use Gauss quadratures instead of the Newton-Cotes guys below).
angs = linspace(-pi, pi, N+1); angs(end) = [];
zz = exp(1i*angs);
H_hat = fftshift(ifft(c_lag.'))*N;
%H_hat = z_trans(c_lag.', w);

%% Recover the impulse response
h = zeros(N, 1);
phi1 = (1 - conj(a).* B_inv_a(zz) )/sqrt(1 - abs(a)^2);
phi2 = B_inv_a(zz);
for k=1:N
    phik = phi1.*phi2.^(k-1);
    h(k) = trapz(H_hat.*conj(phik))/N;
end

% Time ends here
time = toc;
disp(["Impulse response recorey took: ", num2str(time), "seconds"]);

%% Plot results
figure;
subplot(211);
plot(real(h_true), "LineWidth", 2); hold on;
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

inds=1:1000;
figure;
subplot(211);
plot(real(h_true(inds)), "LineWidth", 2); hold on;
plot(real(h(inds)), "r");
axis tight;
grid on;
xlabel("Discrete timestamps");
ylabel("Impulse response values");
legend("True IR", "Recovered IR");
title("Impulse response comparison");
subplot(212);
plot(abs(h(inds) - h_true(inds).'), "LineWidth", 2);
axis tight;
grid on;
xlabel("Discrete timestamps");
ylabel("$| h - \hat{h}|$", 'Interpreter', 'latex');
title("Absolute error");

