% Experiment 3
clc; clear; close all;

% Parameters
Fs = 64; % Sampling frequency
T = 1/Fs;

% Frequencies of the analog signal: 
% x(t) = cos(8pi*t) + cos(16pi*t) + cos(20pi*t)
% f1 = 4 Hz, f2 = 8 Hz, f3 = 10 Hz

N_values = [16, 32, 64];

figure;

for i = 1:length(N_values)
    N = N_values(i);
    
    % Discrete time vector
    n = 0:N-1;
    
    % Discrete signal x(n) obtained by sampling x(t) at t = n*T
    % x(n) = cos(8*pi*n*T) + cos(16*pi*n*T) + cos(20*pi*n*T)
    x = cos(8*pi*n*T) + cos(16*pi*n*T) + cos(20*pi*n*T);
    
    % FFT
    X = fft(x, N);
    
    % Frequency axis for plotting (optional, but good for analysis)
    % f = (0:N-1) * Fs / N; 
    % We will stick to k index for consistency with previous exp, 
    % but maybe label x-axis with Hz is better for "Analog Signal" analysis?
    % The prompt asks for "幅频特性" (Amplitude-Frequency Characteristic).
    % Let's stick to k but mention resolution in title or just plot against k.
    % Given previous instructions, simple stem is preferred.
    
    subplot(3, 1, i);
    stem(0:N-1, abs(X), 'filled');
    title(['|DFT[x_8]| (N=' num2str(N) ', F_s=64Hz)']);
    xlabel('k'); ylabel('Magnitude');
    grid on;
end
