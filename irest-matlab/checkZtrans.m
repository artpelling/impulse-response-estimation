% We can compute the Z transform at any "z" point on the Torus T (boundary of the unit disk). If we
% choose "z" to coincide with the equidistant sampling on T, then we should
% get the same result as fft. This script checks my Ztrans function.

clear all;
close all;

%% Add needed libraries
addpath("Common");

%% Constants
N = 5000; % Number of points in the generating sequence
fi = linspace(-pi, pi, N+1); fi(end)=[];
z = exp(1i*fi);

%% Compute generating sequence
x = randn(1, N);

% Compute Z transform using own function
X = z_trans(x, z);

% Compute using FFT
X_fft = fftshift(ifft(x))*N; % I use positive exponents, hence ifft

%% Plot results
figure;
subplot(131);
plot(fi, real(X)); hold on;
plot(fi, real(X_fft));
subplot(132);
plot(fi, imag(X)); hold on;
plot(fi, imag(X_fft));
subplot(133);
plot(fi, abs(X - X_fft.'));
