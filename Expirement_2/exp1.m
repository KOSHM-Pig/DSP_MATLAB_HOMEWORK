% Experiment 1
clc; clear; close all;

% Define sequences
n_short = 0:3;
x1_base = ones(1, 4); % R4(n)

n_full = 0:7;
x2 = zeros(1, 8);
x3 = zeros(1, 8);

for n = 0:7
    if n <= 3
        x2(n+1) = n + 1;
        x3(n+1) = 4 - n;
    elseif n <= 7
        x2(n+1) = 8 - n;
        x3(n+1) = n - 3;
    end
end

% x1 needs to be defined for at least length 8 context if we want to be consistent, 
% but base info is enough. We will zero pad.
x1 = [x1_base, zeros(1, 4)];

% FFT lengths
N1 = 8;
N2 = 16;

% Compute FFTs
X1_8 = fft(x1, N1);
X1_16 = fft(x1, N2);

X2_8 = fft(x2, N1);
X2_16 = fft(x2, N2);

X3_8 = fft(x3, N1);
X3_16 = fft(x3, N2);

% Plotting
figure;

% x1
subplot(3, 2, 1);
stem(0:N1-1, abs(X1_8), 'filled');
title('|DFT[x_1]| (N=8)');
xlabel('k'); ylabel('Magnitude');
grid on;

subplot(3, 2, 2);
stem(0:N2-1, abs(X1_16), 'filled');
title('|DFT[x_1]| (N=16)');
xlabel('k'); ylabel('Magnitude');
grid on;

% x2
subplot(3, 2, 3);
stem(0:N1-1, abs(X2_8), 'filled');
title('|DFT[x_2]| (N=8)');
xlabel('k'); ylabel('Magnitude');
grid on;

subplot(3, 2, 4);
stem(0:N2-1, abs(X2_16), 'filled');
title('|DFT[x_2]| (N=16)');
xlabel('k'); ylabel('Magnitude');
grid on;

% x3
subplot(3, 2, 5);
stem(0:N1-1, abs(X3_8), 'filled');
title('|DFT[x_3]| (N=8)');
xlabel('k'); ylabel('Magnitude');
grid on;

subplot(3, 2, 6);
stem(0:N2-1, abs(X3_16), 'filled');
title('|DFT[x_3]| (N=16)');
xlabel('k'); ylabel('Magnitude');
grid on;

% Relationship analysis
% x3 is x2 circularly shifted by 4. |X3(k)| = |X2(k)|.
