% 时域采样实验：连续时间信号及其离散采样的频谱分析

% 连续时间信号参数（模拟“模拟信号”）
A=444.128;             % 振幅 A
alpha=50*sqrt(2)*pi;   % 衰减系数 alpha
Omega_0=50*sqrt(2)*pi; % 角频率 Omega_0

% --------- 利用 DFT(FFT) 近似连续时间信号的幅频特性 ----------
% 思路：先用很高的采样频率 Fs_fft 对连续时间信号 xa(t) 采样，
%       得到一个“近似连续”的离散信号，再对它做 FFT 求频谱。

Fs_fft=5000;           % 为近似连续时间用的高采样频率 (Hz)，要远大于 500 Hz
T_fft=0.5;             % 观察时间长度 (s)，足够长可以看到指数衰减过程
N_fft=round(Fs_fft*T_fft); % 采样点数

t_fft=(0:N_fft-1)/Fs_fft;                    % 连续时间 t 的离散取样点
xa=A*exp(-alpha*t_fft).*sin(Omega_0*t_fft);  % xa(t) = A*e^{-alpha t}*sin(Omega_0 t)

X_fft=fft(xa);                 % 对 xa(t) 的采样序列做 N_fft 点 FFT
f_fft=(0:N_fft-1)*Fs_fft/N_fft;% 对应的频率轴 0 ~ Fs_fft-Δf

% 只取 0~500 Hz 这一段来画幅频特性
idx_f=(f_fft>=0)&(f_fft<=500);
f=f_fft(idx_f);                % 0~500 Hz 的频率
Xa_mag=abs(X_fft(idx_f));      % 对应的幅度谱近似 |X_a(f)|

% --------- 公共参数：采样频率、观测时间、FFT 点数 ----------
Fs1=1000;              % 采样频率 1 (Hz)
Fs2=300;               % 采样频率 2 (Hz)
Fs3=200;               % 采样频率 3 (Hz)
Tobs=0.064;            % 观测时间 (s)
M=64;                  % FFT 点数

% ===== Fs = 1000 Hz =====
Ts1=1/Fs1;                             % 采样间隔
N1=round(Tobs*Fs1);                    % 采样点数
n1=0:N1-1;                             % 离散时间序号
t1=n1*Ts1;                             % 对应的连续时间
x1=A*exp(-alpha*t1).*sin(Omega_0*t1);  % 采样得到的序列 x1(n)
X1=fft(x1,M);                          % 对 x1(n) 做 64 点 FFT
k=0:M-1;                               % 频率索引
f1=k*Fs1/M;                            % 对应的频率轴
amp1=abs(X1);                          % 幅度谱
subplot(3,1,1);
plot(f1,amp1);
xlabel('f/Hz');
ylabel('|X_1(f)|');
title('Fs = 1000 Hz 时的幅度谱');

% ===== Fs = 300 Hz =====
Ts2=1/Fs2;
N2=round(Tobs*Fs2);
n2=0:N2-1;
t2=n2*Ts2;
x2=A*exp(-alpha*t2).*sin(Omega_0*t2);  % 采样得到的序列 x2(n)
X2=fft(x2,M);                          % 对 x2(n) 做 64 点 FFT
f2=k*Fs2/M;
amp2=abs(X2);

% ===== Fs = 200 Hz =====
Ts3=1/Fs3;
N3=round(Tobs*Fs3);
n3=0:N3-1;
t3=n3*Ts3;
x3=A*exp(-alpha*t3).*sin(Omega_0*t3);  % 采样得到的序列 x3(n)
X3=fft(x3,M);                          % 对 x3(n) 做 64 点 FFT
f3=k*Fs3/M;
amp3=abs(X3);

% ----------------- 6. 将所有图放在一个 Figure 中 -----------------
% 布局：4 行 2 列
% 第 1 行：连续时间信号 xa(t) 及其幅频特性
% 第 2 行：Fs1=1000 Hz 下的时域采样序列及频谱
% 第 3 行：Fs2=300  Hz 下的时域采样序列及频谱
% 第 4 行：Fs3=200  Hz 下的时域采样序列及频谱

figure;

% 第 1 行左图：连续时间信号 xa(t) 的波形（近似连续）
subplot(4,2,1);
plot(t_fft,xa);
xlabel('t/s');
ylabel('xa(t)');
title('连续时间信号 xa(t)（高采样率近似）');

% 第 1 行右图：xa(t) 的幅频特性
subplot(4,2,2);
plot(f,Xa_mag);
xlabel('f/Hz');
ylabel('|X_a(f)|');
title('利用 FFT 近似得到的 |X_a(f)|');

% 第 2 行左图：Fs1=1000 Hz 时的离散采样序列 x1(n)
subplot(4,2,3);
stem(t1,x1,'filled');
xlabel('t/s');
ylabel('x_1(n)');
title('Fs = 1000 Hz 时的采样序列');

% 第 2 行右图：Fs1=1000 Hz 时的频谱
subplot(4,2,4);
plot(f1,amp1);
xlabel('f/Hz');
ylabel('|X_1(f)|');
title('Fs = 1000 Hz 时的幅度谱');

% 第 3 行左图：Fs2=300 Hz 时的离散采样序列 x2(n)
subplot(4,2,5);
stem(t2,x2,'filled');
xlabel('t/s');
ylabel('x_2(n)');
title('Fs = 300 Hz 时的采样序列');

% 第 3 行右图：Fs2=300 Hz 时的频谱
subplot(4,2,6);
plot(f2,amp2);
xlabel('f/Hz');
ylabel('|X_2(f)|');
title('Fs = 300 Hz 时的幅度谱');

% 第 4 行左图：Fs3=200 Hz 时的离散采样序列 x3(n)
subplot(4,2,7);
stem(t3,x3,'filled');
xlabel('t/s');
ylabel('x_3(n)');
title('Fs = 200 Hz 时的采样序列');

% 第 4 行右图：Fs3=200 Hz 时的频谱
subplot(4,2,8);
plot(f3,amp3);
xlabel('f/Hz');
ylabel('|X_3(f)|');
title('Fs = 200 Hz 时的幅度谱');
