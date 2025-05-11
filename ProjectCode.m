% Load ECG Data
data = readtable('C:\AllData\Semester6\DSP\labs\dspProject\archive\100.csv');
time = data.time_ms / 1000;
ecg_clean = data.MLII;
Fs = 360; % Sampling frequency

% Plot in Time Domain
figure;
plot(time, ecg_clean);
xlabel('Time (s)'); ylabel('Amplitude (mV)');
title('Original ECG Signal - Time Domain');

% Frequency Spectrum
L = length(ecg_clean);
f = Fs * (0:(L/2))/L;
Y = fft(ecg_clean);
P = abs(Y/L);
P1 = P(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure;
plot(f, P1);
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Spectrum');

% Spectrogram
figure;
spectrogram(ecg_clean, 256, 200, 512, Fs, 'yaxis');
title('ECG Spectrogram');
% Simulate Noises
t = time;
powerline = 0.2 * sin(2*pi*50*t);         % 50 Hz
baseline = 0.4 * sin(2*pi*0.5*t);         % 0.5 Hz
emg = 0.1 * randn(size(ecg_clean));       % High-freq noise

% Combine with ECG
noisy_ecg = ecg_clean + powerline + baseline + emg;

% Plot
figure;
plot(t, noisy_ecg);
xlabel('Time (s)'); ylabel('Amplitude');
title('Noisy ECG Signal');

% 1. Notch Filter (50 Hz)
wo = 50/(Fs/2);      % Normalized frequency
bw = wo/35;          % Bandwidth (narrow notch)
[bn, an] = iirnotch(wo, bw);

% 2. High-pass Filter (cutoff = 0.7 Hz)
fc_hp = 0.7 / (Fs/2);    % Normalized cutoff
[bh, ah] = butter(4, fc_hp, 'high');

% 3. Low-pass Filter (cutoff = 40 Hz)
fc_lp = 40 / (Fs/2);     % Normalized cutoff
[bl, al] = butter(4, fc_lp, 'low');

% Apply filters sequentially
filtered_ecg = filter(bn, an, noisy_ecg);     % Notch filter
filtered_ecg = filter(bh, ah, filtered_ecg);  % High-pass filter
filtered_ecg = filter(bl, al, filtered_ecg);  % Low-pass filter

% Time-domain Comparison
figure;
plot(t, ecg_clean, 'g'); hold on;
plot(t, filtered_ecg, 'b');
legend('Original Clean', 'Filtered ECG');
title('Comparison: Clean vs Filtered');
xlabel('Time (s)');

% FFT Comparison
Yf = fft(filtered_ecg);
Pf = abs(Yf/L); Pf1 = Pf(1:L/2+1); Pf1(2:end-1) = 2*Pf1(2:end-1);

figure;
plot(f, P1, 'r--'); hold on;
plot(f, Pf1, 'b');
legend('Noisy ECG', 'Filtered ECG');
title('Frequency Spectrum Comparison');

% Spectrogram of Filtered Signal
figure;
spectrogram(filtered_ecg, 256, 200, 512, Fs, 'yaxis');
title('Filtered ECG Spectrogram');

% RMSE
rmse = sqrt(mean((ecg_clean - filtered_ecg).^2));

% SNR
snr_val = snr(filtered_ecg, filtered_ecg - ecg_clean);

disp(['RMSE: ', num2str(rmse)]);
disp(['SNR (dB): ', num2str(snr_val)]);
