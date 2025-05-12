% Load ECG Data
data = readtable('C:\AllData\Semester6\DSP\labs\dspProject\archive\100.csv');
time = data.time_ms;
Sectime = time / 1000;
ecg_MLII = data.MLII;      % Lead MLII
ecg_V1   = data{:,4};        % Lead V1
Fs = 360;
SamplesToPlot = 10 * Fs;

%% === TIME-DOMAIN PLOT ===
figure;
subplot(2,1,1);
plot(Sectime(1:SamplesToPlot), ecg_MLII(1:SamplesToPlot), 'g');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
title('Original ECG - MLII (Time Domain)');

subplot(2,1,2);
plot(Sectime(1:SamplesToPlot), ecg_V1(1:SamplesToPlot), 'm');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
title('Original ECG - V1 (Time Domain)');

%% === FREQUENCY SPECTRA ===
L = length(ecg_MLII);
f = Fs * (0:(L/2)) / L;

Y1 = fft(ecg_MLII);  Y2 = fft(ecg_V1);
P1 = abs(Y1 / L);    P2 = abs(Y2 / L);
P1 = P1(1:L/2+1);    P2 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P2(2:end-1) = 2*P2(2:end-1);

figure;
subplot(2,1,1);
plot(f, P1, 'g'); title('Frequency Spectrum - MLII');
xlabel('Frequency (Hz)'); ylabel('Magnitude');

subplot(2,1,2);
plot(f, P2, 'm'); title('Frequency Spectrum - V1');
xlabel('Frequency (Hz)'); ylabel('Magnitude');

%% === SPECTROGRAMS ===
figure;
subplot(2,1,1);
spectrogram(ecg_MLII, 256, 200, 512, Fs, 'yaxis');
title('Spectrogram - MLII');

subplot(2,1,2);
spectrogram(ecg_V1, 256, 200, 512, Fs, 'yaxis');
title('Spectrogram - V1');

%% === SIMULATE NOISES FOR BOTH LEADS ===
t = Sectime;
powerline = 0.2 * sin(2*pi*50*t);       % 50 Hz
baseline  = 0.4 * sin(2*pi*0.5*t);      % 0.5 Hz
emg       = 0.1 * randn(size(ecg_MLII));% same size for both leads

noisy_MLII = ecg_MLII + powerline + baseline + emg;
noisy_V1   = ecg_V1   + powerline + baseline + emg;

% Plot noisy signals
figure;
subplot(2,1,1);
plot(t(1:SamplesToPlot), noisy_MLII(1:SamplesToPlot), 'r');
title('Noisy ECG - MLII');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(2,1,2);
plot(t(1:SamplesToPlot), noisy_V1(1:SamplesToPlot), 'r');
title('Noisy ECG - V1');
xlabel('Time (s)'); ylabel('Amplitude');

%% === DESIGN FILTERS ===
wo = 50 / (Fs/2); bw = wo / 35;
[bn, an] = iirnotch(wo, bw);
fc_hp = 0.7 / (Fs/2); [bh, ah] = butter(4, fc_hp, 'high');
fc_lp = 40 / (Fs/2); [bl, al] = butter(4, fc_lp, 'low');

% Filter MLII
filt_MLII = filter(bn, an, noisy_MLII);
filt_MLII = filter(bh, ah, filt_MLII);
filt_MLII = filter(bl, al, filt_MLII);

% Filter V1
filt_V1 = filter(bn, an, noisy_V1);
filt_V1 = filter(bh, ah, filt_V1);
filt_V1 = filter(bl, al, filt_V1);

%% === COMPARISON TIME DOMAIN ===
figure;
subplot(2,1,1);
plot(t, ecg_MLII, 'g'); hold on;
plot(t, filt_MLII, 'b');
legend('Clean', 'Filtered'); title('MLII: Clean vs Filtered');

subplot(2,1,2);
plot(t, ecg_V1, 'm'); hold on;
plot(t, filt_V1, 'b');
legend('Clean', 'Filtered'); title('V1: Clean vs Filtered');

%% === FFT COMPARISON ===
Yf1 = fft(filt_MLII); Yf2 = fft(filt_V1);
Pf1 = abs(Yf1/L); Pf2 = abs(Yf2/L);
Pf1 = Pf1(1:L/2+1); Pf2 = Pf2(1:L/2+1);
Pf1(2:end-1) = 2*Pf1(2:end-1);
Pf2(2:end-1) = 2*Pf2(2:end-1);

figure;
subplot(2,1,1);
plot(f, P1, 'r--'); hold on;
plot(f, Pf1, 'b');
legend('Noisy', 'Filtered'); title('FFT MLII');

subplot(2,1,2);
plot(f, P2, 'r--'); hold on;
plot(f, Pf2, 'b');
legend('Noisy', 'Filtered'); title('FFT V1');

%% === SPECTROGRAM OF FILTERED ===
figure;
subplot(2,1,1);
spectrogram(filt_MLII, 256, 200, 512, Fs, 'yaxis');
title('Filtered Spectrogram - MLII');

subplot(2,1,2);
spectrogram(filt_V1, 256, 200, 512, Fs, 'yaxis');
title('Filtered Spectrogram - V1');

%% === METRICS (RMSE & SNR) ===
rmse_MLII = sqrt(mean((ecg_MLII - filt_MLII).^2));
snr_MLII = snr(filt_MLII, filt_MLII - ecg_MLII);

rmse_V1 = sqrt(mean((ecg_V1 - filt_V1).^2));
snr_V1 = snr(filt_V1, filt_V1 - ecg_V1);

disp(['MLII RMSE: ', num2str(rmse_MLII), ', SNR (dB): ', num2str(snr_MLII)]);
disp(['V1   RMSE: ', num2str(rmse_V1),   ', SNR (dB): ', num2str(snr_V1)]);
