clear; close all; clc;

%% 1. 实验步骤 (2): 信号生成与频谱显示
% 题目要求：给定一个具有加性噪声的信号x(t)，显示x(t)及其频谱
% 采样频率
Fs = 1000;
T = 1/Fs;
t = 0:T:1; % 时间长度 1秒

% 构造单频调幅信号 (目标信号)
% 载波 fc = 60 Hz (在通带 120Hz 内)
% 调制信号 fm = 5 Hz
fc = 60;
fm = 5;
s_am = (1 + 0.5*cos(2*pi*fm*t)) .* cos(2*pi*fc*t);

% 构造高频噪声
% 噪声频率 > 阻带截止频率 150Hz
% 选取 250Hz 和 400Hz 的干扰，加上少量白噪声
noise = 0.3*cos(2*pi*250*t) + 0.2*cos(2*pi*400*t) + 0.05*randn(size(t));

% 合成输入信号
x = s_am + noise;

% 计算输入信号频谱
NFFT = 2048;
f_axis = (0:NFFT/2)*Fs/NFFT;
X = fft(x, NFFT);
Xmag = abs(X(1:NFFT/2+1));
Xmag = Xmag / max(Xmag); % 归一化

% 绘图：输入信号
figure('Name', '输入信号及其频谱');
subplot(2,1,1);
plot(t, x);
title('输入信号 x(t) (含噪声)');
xlabel('时间 (s)'); ylabel('幅值');
grid on; xlim([0 0.5]); % 只显示前0.5秒以便观察

subplot(2,1,2);
plot(f_axis, 20*log10(Xmag + eps));
title('输入信号频谱 (dB)');
xlabel('频率 (Hz)'); ylabel('幅值 (dB)');
grid on; xlim([0 500]);
hold on;
yl = ylim;
plot([120 120], yl, 'g--', 'LineWidth', 1, 'DisplayName', '通带截止 120Hz');
plot([150 150], yl, 'r--', 'LineWidth', 1, 'DisplayName', '阻带截止 150Hz');
legend('show');
hold off;

%% 2. 实验步骤 (3) & (4): FIR 滤波器设计
% 题目要求：设计低通滤波器，提取调幅信号；选择窗函数，调用fir1设计
% 指标参数
fp = 120;       % 通带截止频率 Hz
fs_stop = 150;  % 阻带截止频率 Hz
Rp = 0.1;       % 通带最大衰减 dB
Rs = 60;        % 阻带最小衰减 dB

% 归一化频率 (Nyquist = 1)
Wp = fp / (Fs/2);
Ws = fs_stop / (Fs/2);

% 过渡带宽度
delta_w_norm = Ws - Wp; % 归一化过渡带宽度 (0.3 - 0.24 = 0.06)

% 选择窗函数：Blackman
% Blackman 窗通常能提供约 74dB 的阻带衰减 (主瓣旁瓣比约为 -58dB, 但渐进衰减更快)
% 估算阶数 M (或窗口长度 N = M + 1)
% Blackman 过渡带近似公式: Delta_f_norm = 5.56 / M
% M = 5.56 / delta_w_norm
M_est = ceil(5.56 / delta_w_norm); 
% 确保 M 为偶数 (Type I FIR) 或根据 fir1 自动调整
M = M_est + mod(M_est, 2); % 凑偶数阶，窗口长度为奇数

fprintf('计算得到的滤波器阶数 M = %d (窗口长度 N = %d)\n', M, M+1);

% 截止频率 (归一化)
Wc = (Wp + Ws) / 2;

% 生成 Blackman 窗
win = blackman(M+1);

% 使用 fir1 设计滤波器
b = fir1(M, Wc, win);

% 绘图：滤波器频响
[H, w] = freqz(b, 1, 1024);
magH = 20*log10(abs(H));
freq_hz = w/pi * (Fs/2);

figure('Name', 'FIR 滤波器频响特性');
plot(freq_hz, magH, 'LineWidth', 1.5);
grid on;
title(['FIR 低通滤波器幅频响应 (Blackman 窗, 阶数 M=', num2str(M), ')']);
xlabel('频率 (Hz)'); ylabel('幅度 (dB)');
hold on;
yl = ylim;
plot([fp fp], yl, 'g--', 'LineWidth', 1, 'DisplayName', ['Fp=', num2str(fp), 'Hz']);
plot([fs_stop fs_stop], yl, 'r--', 'LineWidth', 1, 'DisplayName', ['Fs=', num2str(fs_stop), 'Hz']);
plot(xlim, [-Rs -Rs], 'k--', 'LineWidth', 1, 'DisplayName', ['Rs=', num2str(-Rs), 'dB']);
hold off;
ylim([-100 10]);

%% 3. 实验步骤 (4): 快速卷积滤波与结果绘图
% 题目要求：调用fftfilt实现滤波；绘图显示频响、输出幅频特性和时域波形
% fftfilt 利用重叠相加法实现 FFT 卷积
y = fftfilt(b, x);

% 由于滤波带来群延时 (对于线性相位 FIR，延时为 M/2 个样本)
% 为了时域波形对齐，可以平移输出信号 (可选，但在对比时很有用)
delay = M / 2;
t_out = t; % 输出时间轴与输入一致

% 计算输出信号频谱
Y = fft(y, NFFT);
Ymag = abs(Y(1:NFFT/2+1));
Ymag = Ymag / max(Ymag);

% 绘图：滤波结果
figure('Name', '滤波输出结果');
subplot(2,1,1);
plot(t, y);
title('滤波后输出信号 y(t)');
xlabel('时间 (s)'); ylabel('幅值');
grid on; xlim([0 0.5]);

subplot(2,1,2);
plot(f_axis, 20*log10(Ymag + eps));
title('输出信号频谱 (dB)');
xlabel('频率 (Hz)'); ylabel('幅值 (dB)');
grid on; xlim([0 500]);
hold on;
yl = ylim;
plot([120 120], yl, 'g--', 'LineWidth', 1, 'DisplayName', '120Hz');
hold off;

%% 4. 验证对比
figure('Name', '时域波形对比');
plot(t, x, 'color', [0.7 0.7 0.7], 'DisplayName', '输入信号 (含噪)');
hold on;
plot(t, y, 'b', 'LineWidth', 1.5, 'DisplayName', '输出信号 (滤波后)');
plot(t, s_am, 'r--', 'LineWidth', 1, 'DisplayName', '原始AM信号 (参考)');
hold off;
title('时域波形对比');
xlabel('时间 (s)'); ylabel('幅值');
legend;
grid on; xlim([0.2 0.4]); % 放大观察局部
