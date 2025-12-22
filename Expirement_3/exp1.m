clc; clear; close all;

Fs = 4000;
T = 0.125;
t = 0:1/Fs:T-1/Fs;

fm1 = 25;  fc1 = 225;
fm2 = 50;  fc2 = 500;
fm3 = 100; fc3 = 1000;

m1 = cos(2*pi*fm1*t);
m2 = cos(2*pi*fm2*t);
m3 = cos(2*pi*fm3*t);

y1 = m1 .* cos(2*pi*fc1*t);
y2 = m2 .* cos(2*pi*fc2*t);
y3 = m3 .* cos(2*pi*fc3*t);

s = y1 + y2 + y3;

Rp = 0.1;
Rs = 60;

Fp_lp = 280;
Fs_lp = 350;
[n_lp, Wn_lp] = ellipord(Fp_lp/(Fs/2), Fs_lp/(Fs/2), Rp, Rs);
[b_lp, a_lp] = ellip(n_lp, Rp, Rs, Wn_lp, 'low');
s1 = filtfilt(b_lp, a_lp, s);

Fs1_bp = 350;
Fp1_bp = 420;
Fp2_bp = 580;
Fs2_bp = 800;
[n_bp, Wn_bp] = ellipord([Fp1_bp Fp2_bp]/(Fs/2), [Fs1_bp Fs2_bp]/(Fs/2), Rp, Rs);
[b_bp, a_bp] = ellip(n_bp, Rp, Rs, Wn_bp, 'bandpass');
s2 = filtfilt(b_bp, a_bp, s);

Fs_hp = 800;
Fp_hp = 850;
[n_hp, Wn_hp] = ellipord(Fp_hp/(Fs/2), Fs_hp/(Fs/2), Rp, Rs);
[b_hp, a_hp] = ellip(n_hp, Rp, Rs, Wn_hp, 'high');
s3 = filtfilt(b_hp, a_hp, s);

disp('滤波器频率指标(Hz):');
disp(['低通:  Fp=', num2str(Fp_lp), '  Fs=', num2str(Fs_lp)]);
disp(['带通:  Fp=[', num2str(Fp1_bp), ' ', num2str(Fp2_bp), ']  Fs=[', num2str(Fs1_bp), ' ', num2str(Fs2_bp), ']']);
disp(['高通:  Fp=', num2str(Fp_hp), '  Fs=', num2str(Fs_hp)]);

NFFT = 4096;
S = fft(s, NFFT);
f = (0:NFFT/2) * Fs / NFFT;
Smag = abs(S(1:NFFT/2+1));
Smag = Smag / max(Smag);

figure;
subplot(2, 1, 1);
plot(t, s, 'LineWidth', 1);
xlabel('t/s');
ylabel('s(t)');
title('(a) s(t) 的波形');
grid on;
xlim([t(1), t(end)]);

subplot(2, 1, 2);
stem(f, Smag, 'filled');
xlabel('f/Hz');
ylabel('幅度');
title('(b) s(t) 的频谱');
grid on;
xlim([0, Fs/2]);

S1 = fft(s1, NFFT);
S2 = fft(s2, NFFT);
S3 = fft(s3, NFFT);
S1mag = abs(S1(1:NFFT/2+1)); S1mag = S1mag / max(S1mag);
S2mag = abs(S2(1:NFFT/2+1)); S2mag = S2mag / max(S2mag);
S3mag = abs(S3(1:NFFT/2+1)); S3mag = S3mag / max(S3mag);

figure;
subplot(3, 2, 1);
plot(t, s1, 'LineWidth', 1);
xlabel('t/s'); ylabel('s_1(t)'); title('低通分离后的时域波形'); grid on; xlim([t(1), t(end)]);

subplot(3, 2, 2);
stem(f, S1mag, 'filled');
xlabel('f/Hz'); ylabel('幅度'); title('低通分离后的频谱'); grid on; xlim([0, Fs/2]);

subplot(3, 2, 3);
plot(t, s2, 'LineWidth', 1);
xlabel('t/s'); ylabel('s_2(t)'); title('带通分离后的时域波形'); grid on; xlim([t(1), t(end)]);

subplot(3, 2, 4);
stem(f, S2mag, 'filled');
xlabel('f/Hz'); ylabel('幅度'); title('带通分离后的频谱'); grid on; xlim([0, Fs/2]);

subplot(3, 2, 5);
plot(t, s3, 'LineWidth', 1);
xlabel('t/s'); ylabel('s_3(t)'); title('高通分离后的时域波形'); grid on; xlim([t(1), t(end)]);

subplot(3, 2, 6);
stem(f, S3mag, 'filled');
xlabel('f/Hz'); ylabel('幅度'); title('高通分离后的频谱'); grid on; xlim([0, Fs/2]);

[n_lp_butt, Wn_lp_butt] = buttord(Fp_lp/(Fs/2), Fs_lp/(Fs/2), Rp, Rs);
[b_lp_butt, a_lp_butt] = butter(n_lp_butt, Wn_lp_butt, 'low');
[n_lp_cheb, Wn_lp_cheb] = cheb1ord(Fp_lp/(Fs/2), Fs_lp/(Fs/2), Rp, Rs);
[b_lp_cheb, a_lp_cheb] = cheby1(n_lp_cheb, Rp, Wn_lp_cheb, 'low');

[n_bp_butt, Wn_bp_butt] = buttord([Fp1_bp Fp2_bp]/(Fs/2), [Fs1_bp Fs2_bp]/(Fs/2), Rp, Rs);
[b_bp_butt, a_bp_butt] = butter(n_bp_butt, Wn_bp_butt, 'bandpass');
[n_bp_cheb, Wn_bp_cheb] = cheb1ord([Fp1_bp Fp2_bp]/(Fs/2), [Fs1_bp Fs2_bp]/(Fs/2), Rp, Rs);
[b_bp_cheb, a_bp_cheb] = cheby1(n_bp_cheb, Rp, Wn_bp_cheb, 'bandpass');

[n_hp_butt, Wn_hp_butt] = buttord(Fp_hp/(Fs/2), Fs_hp/(Fs/2), Rp, Rs);
[b_hp_butt, a_hp_butt] = butter(n_hp_butt, Wn_hp_butt, 'high');
[n_hp_cheb, Wn_hp_cheb] = cheb1ord(Fp_hp/(Fs/2), Fs_hp/(Fs/2), Rp, Rs);
[b_hp_cheb, a_hp_cheb] = cheby1(n_hp_cheb, Rp, Wn_hp_cheb, 'high');

[H_lp_butt, w] = freqz(b_lp_butt, a_lp_butt, NFFT, Fs);
H_lp_cheb = freqz(b_lp_cheb, a_lp_cheb, NFFT, Fs);
H_bp_butt = freqz(b_bp_butt, a_bp_butt, NFFT, Fs);
H_bp_cheb = freqz(b_bp_cheb, a_bp_cheb, NFFT, Fs);
H_hp_butt = freqz(b_hp_butt, a_hp_butt, NFFT, Fs);
H_hp_cheb = freqz(b_hp_cheb, a_hp_cheb, NFFT, Fs);

L_lp_butt = -20*log10(abs(H_lp_butt) + eps);
L_lp_cheb = -20*log10(abs(H_lp_cheb) + eps);
L_bp_butt = -20*log10(abs(H_bp_butt) + eps);
L_bp_cheb = -20*log10(abs(H_bp_cheb) + eps);
L_hp_butt = -20*log10(abs(H_hp_butt) + eps);
L_hp_cheb = -20*log10(abs(H_hp_cheb) + eps);

figure;
subplot(3, 1, 1);
plot(w, L_lp_butt, 'LineWidth', 1); hold on;
plot(w, L_lp_cheb, 'LineWidth', 1); hold off;
grid on; xlim([0, Fs/2]); ylim([0, 100]);
xlabel('f/Hz'); ylabel('损耗/dB'); title('低通滤波器损耗函数对比'); legend('butter', 'cheby1');

subplot(3, 1, 2);
plot(w, L_bp_butt, 'LineWidth', 1); hold on;
plot(w, L_bp_cheb, 'LineWidth', 1); hold off;
grid on; xlim([0, Fs/2]); ylim([0, 100]);
xlabel('f/Hz'); ylabel('损耗/dB'); title('带通滤波器损耗函数对比'); legend('butter', 'cheby1');

subplot(3, 1, 3);
plot(w, L_hp_butt, 'LineWidth', 1); hold on;
plot(w, L_hp_cheb, 'LineWidth', 1); hold off;
grid on; xlim([0, Fs/2]); ylim([0, 100]);
xlabel('f/Hz'); ylabel('损耗/dB'); title('高通滤波器损耗函数对比'); legend('butter', 'cheby1');

idx_lp_pass = (w <= Fp_lp);
idx_lp_stop = (w >= Fs_lp);

idx_bp_pass = (w >= Fp1_bp) & (w <= Fp2_bp);
idx_bp_stop = (w <= Fs1_bp) | (w >= Fs2_bp);

idx_hp_pass = (w >= Fp_hp);
idx_hp_stop = (w <= Fs_hp);

disp('butter/cheby1 设计的滤波器阶数对比:');
disp(['低通:  butter=', num2str(n_lp_butt), '  cheby1=', num2str(n_lp_cheb)]);
disp(['带通:  butter=', num2str(n_bp_butt), '  cheby1=', num2str(n_bp_cheb)]);
disp(['高通:  butter=', num2str(n_hp_butt), '  cheby1=', num2str(n_hp_cheb)]);

disp('损耗指标验证(单位:dB, 越大表示衰减越强):');
disp(['低通 butter:  通带最大损耗=', num2str(max(L_lp_butt(idx_lp_pass))), '  阻带最小损耗=', num2str(min(L_lp_butt(idx_lp_stop)))]);
disp(['低通 cheby1:  通带最大损耗=', num2str(max(L_lp_cheb(idx_lp_pass))), '  阻带最小损耗=', num2str(min(L_lp_cheb(idx_lp_stop)))]);
disp(['带通 butter:  通带最大损耗=', num2str(max(L_bp_butt(idx_bp_pass))), '  阻带最小损耗=', num2str(min(L_bp_butt(idx_bp_stop)))]);
disp(['带通 cheby1:  通带最大损耗=', num2str(max(L_bp_cheb(idx_bp_pass))), '  阻带最小损耗=', num2str(min(L_bp_cheb(idx_bp_stop)))]);
disp(['高通 butter:  通带最大损耗=', num2str(max(L_hp_butt(idx_hp_pass))), '  阻带最小损耗=', num2str(min(L_hp_butt(idx_hp_stop)))]);
disp(['高通 cheby1:  通带最大损耗=', num2str(max(L_hp_cheb(idx_hp_pass))), '  阻带最小损耗=', num2str(min(L_hp_cheb(idx_hp_stop)))]);

y1n = filter(b_lp_cheb, a_lp_cheb, s);
y2n = filter(b_bp_cheb, a_bp_cheb, s);
y3n = filter(b_hp_cheb, a_hp_cheb, s);

figure;
subplot(3, 1, 1);
plot(t, y1n, 'LineWidth', 1);
xlabel('t/s'); ylabel('y_1(n)'); title('低通滤波分离的 y_1(n) 时域波形'); grid on; xlim([t(1), t(end)]);

subplot(3, 1, 2);
plot(t, y2n, 'LineWidth', 1);
xlabel('t/s'); ylabel('y_2(n)'); title('带通滤波分离的 y_2(n) 时域波形'); grid on; xlim([t(1), t(end)]);

subplot(3, 1, 3);
plot(t, y3n, 'LineWidth', 1);
xlabel('t/s'); ylabel('y_3(n)'); title('高通滤波分离的 y_3(n) 时域波形'); grid on; xlim([t(1), t(end)]);
