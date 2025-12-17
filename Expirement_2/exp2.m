% Experiment 2
clc; clear; close all;

% Define N
N1 = 8;
N2 = 16;

% Define n vectors
n1 = 0:N1-1;
n2 = 0:N2-1;

% Define signals
% x4(n) = cos(pi/4 * n)
x4_8 = cos(pi/4 * n1);
x4_16 = cos(pi/4 * n2);

% x5(n) = cos(pi/4 * n) + cos(pi/8 * n)
x5_8 = cos(pi/4 * n1) + cos(pi/8 * n1);
x5_16 = cos(pi/4 * n2) + cos(pi/8 * n2);

% Compute FFT
X4_8 = fft(x4_8, N1);
X4_16 = fft(x4_16, N2);

X5_8 = fft(x5_8, N1);
X5_16 = fft(x5_16, N2);

% Plotting
figure; 

% x4
subplot(2, 2, 1);
stem(0:N1-1, abs(X4_8), 'filled');
title('|DFT[x_4]| (N=8)');
xlabel('k'); ylabel('Magnitude');
grid on;

subplot(2, 2, 2);
stem(0:N2-1, abs(X4_16), 'filled');
title('|DFT[x_4]| (N=16)');
xlabel('k'); ylabel('Magnitude');
grid on;

% x5
subplot(2, 2, 3);
stem(0:N1-1, abs(X5_8), 'filled');
title('|DFT[x_5]| (N=8)');
xlabel('k'); ylabel('Magnitude');
grid on;

subplot(2, 2, 4);
stem(0:N2-1, abs(X5_16), 'filled');
title('|DFT[x_5]| (N=16)');
xlabel('k'); ylabel('Magnitude');
grid on;
