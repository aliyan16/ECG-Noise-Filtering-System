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
P = abs(Y/L); P1 = P(1:L/2+1); P1(2:end-1) = 2*P1(2:end-1);

figure;
plot(f, P1);
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Spectrum');

% Spectrogram
figure;
spectrogram(ecg_clean, 256, 200, 512, Fs, 'yaxis');
title('ECG Spectrogram');
