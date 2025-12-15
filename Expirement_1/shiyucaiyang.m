% 时域采样实验：连续时间信号及其离散采样的频谱分析

% 连续时间信号参数
A=444.128;             % 振幅 A
alpha=50*sqrt(2)*pi;   % 衰减系数 alpha
Omega_0=50*sqrt(2)*pi; % 角频率 Omega_0

% --------- 连续时间信号的幅频特性 |Xa(f)| ----------
f=0:1:500;             % 频率轴 0~500 Hz
Omega=2*pi*f;          % 对应的角频率轴 Omega
Xa=A*Omega_0./sqrt((alpha^2-Omega.^2+Omega_0^2).^2+(2*alpha.*Omega).^2); % 幅频特性公式
figure;
plot(f,Xa);            % 绘制 |Xa(f)| 曲线
xlabel('f/Hz');
ylabel('|X_a(f)|');
title('连续时间信号 X_a(f) 的幅频特性');

% --------- 公共参数：采样频率、观测时间、FFT 点数 ----------
Fs1=1000;              % 采样频率 1 (Hz)
Fs2=300;               % 采样频率 2 (Hz)
Fs3=200;               % 采样频率 3 (Hz)
Tobs=0.064;            % 观测时间 (s)
M=64;                  % FFT 点数

figure;

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
subplot(3,1,2);
plot(f2,amp2);
xlabel('f/Hz');
ylabel('|X_2(f)|');
title('Fs = 300 Hz 时的幅度谱');

% ===== Fs = 200 Hz =====
Ts3=1/Fs3;
N3=round(Tobs*Fs3);
n3=0:N3-1;
t3=n3*Ts3;
x3=A*exp(-alpha*t3).*sin(Omega_0*t3);  % 采样得到的序列 x3(n)
X3=fft(x3,M);                          % 对 x3(n) 做 64 点 FFT
f3=k*Fs3/M;
amp3=abs(X3);
subplot(3,1,3);
plot(f3,amp3);
xlabel('f/Hz');
ylabel('|X_3(f)|');
title('Fs = 200 Hz 时的幅度谱');
