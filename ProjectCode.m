% Load ECG Data
data = readtable('C:\AllData\Semester6\DSP\labs\dspProject\archive\100.csv'); % Make sure this path is correct
time = data.time_ms;
Sectime = time / 1000;
ecg_MLII = data.MLII;      % Lead MLII - Keep as original clean signal
ecg_V1   = data{:,4};        % Lead V1 - Keep as original clean signal
Fs = 360;                   % Sampling Frequency for MIT-BIH
L = length(ecg_MLII);       % Signal Length
t = Sectime;                % Time vector

% Select shorter segment for plotting and some analysis if needed
SamplesToPlot = 10 * Fs; % Plot 10 seconds
PlotEndSample = min(SamplesToPlot, L); % Ensure we don't exceed signal length

%% === 1a. Initial Analysis: Clean Signal Plots ===
fprintf('=== 1a. Analyzing Clean Signal ===\n');
figure('Name', 'Clean Signal Analysis');
% TIME-DOMAIN PLOT
subplot(2,2,1);
plot(t(1:PlotEndSample), ecg_MLII(1:PlotEndSample), 'g');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
title('Original ECG - MLII (Time Domain)');
grid on;

subplot(2,2,2);
plot(t(1:PlotEndSample), ecg_V1(1:PlotEndSample), 'm');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
title('Original ECG - V1 (Time Domain)');
grid on;

% FREQUENCY SPECTRA
f = Fs * (0:(L/2)) / L;
Y1 = fft(ecg_MLII);  Y2 = fft(ecg_V1);
P1 = abs(Y1 / L);    P2 = abs(Y2 / L);
P1_clean = P1(1:L/2+1);    P2_clean = P2(1:L/2+1);
P1_clean(2:end-1) = 2*P1_clean(2:end-1);
P2_clean(2:end-1) = 2*P2_clean(2:end-1);

subplot(2,2,3);
plot(f, P1_clean, 'g'); title('Clean Frequency Spectrum - MLII');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
xlim([0 Fs/2]); grid on;

subplot(2,2,4);
plot(f, P2_clean, 'm'); title('Clean Frequency Spectrum - V1');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
xlim([0 Fs/2]); grid on;

% SPECTROGRAMS
figure('Name', 'Clean Signal Spectrograms');
subplot(2,1,1);
spectrogram(ecg_MLII, hamming(256), 200, 512, Fs, 'yaxis');
title('Clean Spectrogram - MLII');

subplot(2,1,2);
spectrogram(ecg_V1, hamming(256), 200, 512, Fs, 'yaxis');
title('Clean Spectrogram - V1');

%% === 1b. Simulate Noises and Analyze Noisy Signal ===
fprintf('=== 1b. Simulating Noise and Analyzing Noisy Signal ===\n');
% Simulate Noises
powerline_freq = 50; % Set to 50 or 60 Hz depending on region
powerline = 0.2 * sin(2*pi*powerline_freq*t); % Powerline interference
baseline  = 0.4 * sin(2*pi*0.5*t);         % Baseline wander (0.5 Hz)
emg_noise_power = 0.1;                     % Adjust noise power if needed
emg       = sqrt(emg_noise_power) * randn(size(ecg_MLII)); % EMG-like noise

% Create Noisy Signals
noisy_MLII = ecg_MLII + powerline + baseline + emg;
noisy_V1   = ecg_V1   + powerline + baseline + emg;

% Plot noisy signals (Time Domain)
figure('Name', 'Noisy Signal Analysis');
subplot(2,2,1);
plot(t(1:PlotEndSample), noisy_MLII(1:PlotEndSample), 'r');
title('Noisy ECG - MLII (Time Domain)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

subplot(2,2,2);
plot(t(1:PlotEndSample), noisy_V1(1:PlotEndSample), 'r');
title('Noisy ECG - V1 (Time Domain)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

% Analyze Noisy Signals (Frequency Domain)
Yn1 = fft(noisy_MLII);  Yn2 = fft(noisy_V1);
Pn1 = abs(Yn1 / L);    Pn2 = abs(Yn2 / L);
Pn1_noisy = Pn1(1:L/2+1);    Pn2_noisy = Pn2(1:L/2+1);
Pn1_noisy(2:end-1) = 2*Pn1_noisy(2:end-1);
Pn2_noisy(2:end-1) = 2*Pn2_noisy(2:end-1);

subplot(2,2,3);
plot(f, Pn1_noisy, 'r'); title('Noisy Frequency Spectrum - MLII');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
xlim([0 Fs/2]); grid on;

subplot(2,2,4);
plot(f, Pn2_noisy, 'r'); title('Noisy Frequency Spectrum - V1');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
xlim([0 Fs/2]); grid on;

% Spectrograms of Noisy Signals
figure('Name', 'Noisy Signal Spectrograms');
subplot(2,1,1);
spectrogram(noisy_MLII, hamming(256), 200, 512, Fs, 'yaxis');
title('Noisy Spectrogram - MLII');

subplot(2,1,2);
spectrogram(noisy_V1, hamming(256), 200, 512, Fs, 'yaxis');
title('Noisy Spectrogram - V1');

%% === 2. Filter Design & Application (Static Filters) ===
% Using the designs from your original code (assumed from FDA Tool)
fprintf('=== 2. Designing and Applying Static Filters ===\n');

% Notch Filter for Powerline
wo_notch = powerline_freq / (Fs/2); % Normalize frequency
bw_notch = wo_notch / 35;      % Bandwidth (Q factor = 35)
[bn, an] = iirnotch(wo_notch, bw_notch);
fprintf('Notch Filter (IIR) Coefficients (Order %d):\n', max(length(an), length(bn))-1);
% disp(' b_notch:'); disp(bn); disp(' a_notch:'); disp(an);

% High-pass Filter for Baseline Wander
fc_hp = 0.7; % Cut-off frequency in Hz (adjust as needed)
order_hp = 4;
[bh, ah] = butter(order_hp, fc_hp / (Fs/2), 'high');
fprintf('High-pass Filter (Butterworth Order %d, fc=%.1f Hz) Coefficients:\n', order_hp, fc_hp);
% disp(' b_hp:'); disp(bh); disp(' a_hp:'); disp(ah);

% Low-pass Filter for High-Frequency Noise
fc_lp = 40; % Cut-off frequency in Hz (adjust as needed, e.g., 40-50 Hz)
order_lp = 4;
[bl, al] = butter(order_lp, fc_lp / (Fs/2), 'low');
fprintf('Low-pass Filter (Butterworth Order %d, fc=%.1f Hz) Coefficients:\n', order_lp, fc_lp);
% disp(' b_lp:'); disp(bl); disp(' a_lp:'); disp(al);

% Apply filters sequentially
fprintf('Applying static filters...\n');
filt_MLII_static = filter(bn, an, noisy_MLII);      % Remove powerline
filt_MLII_static = filter(bh, ah, filt_MLII_static); % Remove baseline wander
filt_MLII_static = filter(bl, al, filt_MLII_static); % Remove high-freq noise

filt_V1_static = filter(bn, an, noisy_V1);          % Remove powerline
filt_V1_static = filter(bh, ah, filt_V1_static);     % Remove baseline wander
filt_V1_static = filter(bl, al, filt_V1_static);     % Remove high-freq noise
fprintf('Static filtering complete.\n');

%% === 3. Adaptive Filtering (LMS for Powerline Noise Cancellation) ===
fprintf('=== 3. Applying Adaptive Filter (LMS) ===\n');
% We'll use LMS to remove the powerline component adaptively.
% Requires a reference input correlated with the noise to be removed.

% Create reference input: a sine wave at powerline frequency
ref_powerline = 0.2 * sin(2*pi*powerline_freq*t); % Use the same amplitude for reference? Or just sin? Let's try unit amplitude first.
ref_powerline_norm = sin(2*pi*powerline_freq*t); % Normalized reference

% LMS parameters
lms_order = 32;     % Filter order (number of taps); Adjust as needed
lms_mu = 0.001;     % Step size (learning rate); ** CRITICAL - Tune this! **
                    % Start small. If slow convergence, increase slightly. If unstable, decrease.

% Initialize LMS filter object
lms = dsp.LMSFilter('Length', lms_order, 'StepSize', lms_mu);

% Process MLII lead with LMS
% Input: noisy signal (primary input)
% Reference input: signal correlated with noise (powerline ref)
% Output 'y': estimated noise (estimated powerline interference)
% Output 'e': error signal (primary input - estimated noise) = hopefully cleaner signal
fprintf('Running LMS filter for MLII lead...\n');
[y_lms_mlii, e_lms_mlii] = lms(noisy_MLII, ref_powerline_norm);
fprintf('LMS filtering for MLII complete.\n');

% For V1 lead, reset the filter state and run again (or use a separate filter)
reset(lms); % Reset states if using the same object
fprintf('Running LMS filter for V1 lead...\n');
[y_lms_v1, e_lms_v1] = lms(noisy_V1, ref_powerline_norm);
fprintf('LMS filtering for V1 complete.\n');

% Optional: Apply the static HPF and LPF to the LMS output signal
% This creates a combined adaptive + static approach
filt_MLII_lms_combined = filter(bh, ah, e_lms_mlii); % Remove baseline wander
filt_MLII_lms_combined = filter(bl, al, filt_MLII_lms_combined); % Remove high-freq noise

filt_V1_lms_combined = filter(bh, ah, e_lms_v1);     % Remove baseline wander
filt_V1_lms_combined = filter(bl, al, filt_V1_lms_combined);     % Remove high-freq noise

%% === 4. Evaluation & Comparison ===
fprintf('=== 4. Evaluating Filter Performance ===\n');

% --- 4a. Visual Comparison (Time Domain) ---
figure('Name', 'Time Domain Comparison: Clean vs Static vs LMS');
% MLII Lead
subplot(2,1,1);
plot(t(1:PlotEndSample), ecg_MLII(1:PlotEndSample), 'g', 'LineWidth', 1); hold on;
plot(t(1:PlotEndSample), noisy_MLII(1:PlotEndSample), 'r:', 'LineWidth', 0.5);
plot(t(1:PlotEndSample), filt_MLII_static(1:PlotEndSample), 'b', 'LineWidth', 1);
plot(t(1:PlotEndSample), filt_MLII_lms_combined(1:PlotEndSample), 'c--', 'LineWidth', 1); % LMS+Static combined
legend('Clean', 'Noisy', 'Static Filtered', 'LMS+Static Filtered');
title('MLII: Clean vs Noisy vs Filtered Signals');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
grid on;

% V1 Lead
subplot(2,1,2);
plot(t(1:PlotEndSample), ecg_V1(1:PlotEndSample), 'm', 'LineWidth', 1); hold on;
plot(t(1:PlotEndSample), noisy_V1(1:PlotEndSample), 'r:', 'LineWidth', 0.5);
plot(t(1:PlotEndSample), filt_V1_static(1:PlotEndSample), 'b', 'LineWidth', 1);
plot(t(1:PlotEndSample), filt_V1_lms_combined(1:PlotEndSample), 'c--', 'LineWidth', 1); % LMS+Static combined
legend('Clean', 'Noisy', 'Static Filtered', 'LMS+Static Filtered');
title('V1: Clean vs Noisy vs Filtered Signals');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
grid on;

% --- 4b. Frequency Domain Comparison (FFT) ---
% Calculate FFT for filtered signals
Yf_s1 = fft(filt_MLII_static); Yf_s2 = fft(filt_V1_static);
Pf_s1 = abs(Yf_s1/L); Pf_s2 = abs(Yf_s2/L);
Pf_s1 = Pf_s1(1:L/2+1); Pf_s2 = Pf_s2(1:L/2+1);
Pf_s1(2:end-1) = 2*Pf_s1(2:end-1);
Pf_s2(2:end-1) = 2*Pf_s2(2:end-1);

Yf_l1 = fft(filt_MLII_lms_combined); Yf_l2 = fft(filt_V1_lms_combined);
Pf_l1 = abs(Yf_l1/L); Pf_l2 = abs(Yf_l2/L);
Pf_l1 = Pf_l1(1:L/2+1); Pf_l2 = Pf_l2(1:L/2+1);
Pf_l1(2:end-1) = 2*Pf_l1(2:end-1);
Pf_l2(2:end-1) = 2*Pf_l2(2:end-1);

figure('Name', 'Frequency Domain Comparison (FFT)');
subplot(2,1,1);
plot(f, Pn1_noisy, 'r--', 'LineWidth', 0.5); hold on; % Noisy FFT
plot(f, Pf_s1, 'b', 'LineWidth', 1);             % Static Filtered FFT
plot(f, Pf_l1, 'c', 'LineWidth', 1);             % LMS+Static Filtered FFT
plot(f, P1_clean, 'g:', 'LineWidth', 1.5);        % Clean FFT for reference
legend('Noisy', 'Static Filtered', 'LMS+Static Filtered', 'Clean (Reference)');
title('FFT Comparison - MLII');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
xlim([0 Fs/2]); grid on;

subplot(2,1,2);
plot(f, Pn2_noisy, 'r--', 'LineWidth', 0.5); hold on; % Noisy FFT
plot(f, Pf_s2, 'b', 'LineWidth', 1);             % Static Filtered FFT
plot(f, Pf_l2, 'c', 'LineWidth', 1);             % LMS+Static Filtered FFT
plot(f, P2_clean, 'm:', 'LineWidth', 1.5);        % Clean FFT for reference
legend('Noisy', 'Static Filtered', 'LMS+Static Filtered', 'Clean (Reference)');
title('FFT Comparison - V1');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
xlim([0 Fs/2]); grid on;


% --- 4c. Spectrogram Comparison ---
figure('Name', 'Filtered Signal Spectrograms');
subplot(2,2,1);
spectrogram(filt_MLII_static, hamming(256), 200, 512, Fs, 'yaxis');
title('Static Filtered Spectrogram - MLII');

subplot(2,2,2);
spectrogram(filt_V1_static, hamming(256), 200, 512, Fs, 'yaxis');
title('Static Filtered Spectrogram - V1');

subplot(2,2,3);
spectrogram(filt_MLII_lms_combined, hamming(256), 200, 512, Fs, 'yaxis');
title('LMS+Static Filtered Spectrogram - MLII');

subplot(2,2,4);
spectrogram(filt_V1_lms_combined, hamming(256), 200, 512, Fs, 'yaxis');
title('LMS+Static Filtered Spectrogram - V1');


% --- 4d. Quantitative Metrics (SNR, RMSE) ---
fprintf('--- Quantitative Metrics ---\n');

% Calculate noise power for SNR calculation
noise_power_in_mlii = mean((noisy_MLII - ecg_MLII).^2);
noise_power_in_v1   = mean((noisy_V1 - ecg_V1).^2);

signal_power_mlii = mean(ecg_MLII.^2);
signal_power_v1   = mean(ecg_V1.^2);

% Input SNR
snr_in_mlii = 10 * log10(signal_power_mlii / noise_power_in_mlii);
snr_in_v1   = 10 * log10(signal_power_v1 / noise_power_in_v1);

% Output SNR (Static Filter)
noise_power_static_mlii = mean((filt_MLII_static - ecg_MLII).^2);
noise_power_static_v1   = mean((filt_V1_static - ecg_V1).^2);
snr_out_static_mlii = 10 * log10(signal_power_mlii / noise_power_static_mlii);
snr_out_static_v1   = 10 * log10(signal_power_v1 / noise_power_static_v1);
snr_improvement_static_mlii = snr_out_static_mlii - snr_in_mlii;
snr_improvement_static_v1   = snr_out_static_v1 - snr_in_v1;

% Output SNR (LMS+Static Filter)
noise_power_lms_mlii = mean((filt_MLII_lms_combined - ecg_MLII).^2);
noise_power_lms_v1   = mean((filt_V1_lms_combined - ecg_V1).^2);
snr_out_lms_mlii = 10 * log10(signal_power_mlii / noise_power_lms_mlii);
snr_out_lms_v1   = 10 * log10(signal_power_v1 / noise_power_lms_v1);
snr_improvement_lms_mlii = snr_out_lms_mlii - snr_in_mlii;
snr_improvement_lms_v1   = snr_out_lms_v1 - snr_in_v1;

% RMSE
rmse_static_mlii = sqrt(noise_power_static_mlii);
rmse_static_v1   = sqrt(noise_power_static_v1);
rmse_lms_mlii = sqrt(noise_power_lms_mlii);
rmse_lms_v1   = sqrt(noise_power_lms_v1);

% Display Metrics
fprintf('Lead MLII:\n');
fprintf('  Input SNR: %.2f dB\n', snr_in_mlii);
fprintf('  Static Filter Output SNR: %.2f dB (Improvement: %.2f dB), RMSE: %.4f\n', snr_out_static_mlii, snr_improvement_static_mlii, rmse_static_mlii);
fprintf('  LMS+Static Filter Output SNR: %.2f dB (Improvement: %.2f dB), RMSE: %.4f\n', snr_out_lms_mlii, snr_improvement_lms_mlii, rmse_lms_mlii);

fprintf('Lead V1:\n');
fprintf('  Input SNR: %.2f dB\n', snr_in_v1);
fprintf('  Static Filter Output SNR: %.2f dB (Improvement: %.2f dB), RMSE: %.4f\n', snr_out_static_v1, snr_improvement_static_v1, rmse_static_v1);
fprintf('  LMS+Static Filter Output SNR: %.2f dB (Improvement: %.2f dB), RMSE: %.4f\n', snr_out_lms_v1, snr_improvement_lms_v1, rmse_lms_v1);


% --- 4e. R-Peak Detection & Morphological Preservation ---
fprintf('--- R-Peak Detection Metrics ---\n');

% Simple R-peak detection function using findpeaks
% NOTE: This is a *very* basic detector. Robust detection often needs more steps
% (differentiation, squaring, moving window integration). Adjust parameters as needed.
min_peak_height_factor = 0.4; % Adjust based on signal amplitude
min_peak_dist_sec = 0.25; % Minimum distance between peaks in seconds (physiological refractory period)
min_peak_dist_samples = round(min_peak_dist_sec * Fs);

% Detect peaks on Clean signal (Ground Truth)
[pks_clean_mlii, locs_clean_mlii] = findpeaks(ecg_MLII, 'MinPeakHeight', min_peak_height_factor*max(ecg_MLII), 'MinPeakDistance', min_peak_dist_samples);
[pks_clean_v1, locs_clean_v1] = findpeaks(ecg_V1, 'MinPeakHeight', min_peak_height_factor*max(ecg_V1), 'MinPeakDistance', min_peak_dist_samples);
fprintf('Detected %d R-peaks in clean MLII\n', length(locs_clean_mlii));
fprintf('Detected %d R-peaks in clean V1\n', length(locs_clean_v1));

% Detect peaks on Statically Filtered signal
% Use filtered signal's max for threshold, or potentially clean signal's max
[pks_static_mlii, locs_static_mlii] = findpeaks(filt_MLII_static, 'MinPeakHeight', min_peak_height_factor*max(filt_MLII_static), 'MinPeakDistance', min_peak_dist_samples);
[pks_static_v1, locs_static_v1] = findpeaks(filt_V1_static, 'MinPeakHeight', min_peak_height_factor*max(filt_V1_static), 'MinPeakDistance', min_peak_dist_samples);

% Detect peaks on LMS+Static Filtered signal
[pks_lms_mlii, locs_lms_mlii] = findpeaks(filt_MLII_lms_combined, 'MinPeakHeight', min_peak_height_factor*max(filt_MLII_lms_combined), 'MinPeakDistance', min_peak_dist_samples);
[pks_lms_v1, locs_lms_v1] = findpeaks(filt_V1_lms_combined, 'MinPeakHeight', min_peak_height_factor*max(filt_V1_lms_combined), 'MinPeakDistance', min_peak_dist_samples);


% Function to evaluate detection performance
function [Se, P_plus, TP, FP, FN] = evaluate_rpeaks(detected_locs, clean_locs, tolerance_samples)
    TP = 0;
    FP = 0;
    FN = 0;
    detected_matched = false(size(detected_locs));

    for i = 1:length(clean_locs)
        clean_loc = clean_locs(i);
        % Find detected peaks within the tolerance window
        found_match = false;
        for j = 1:length(detected_locs)
            if ~detected_matched(j) && abs(detected_locs(j) - clean_loc) <= tolerance_samples
                TP = TP + 1;
                detected_matched(j) = true; % Mark as matched
                found_match = true;
                break; % Match found for this clean peak, move to next clean peak
            end
        end
        if ~found_match
            FN = FN + 1; % Clean peak missed
        end
    end
    % Any remaining unmatched detected peaks are false positives
    FP = sum(~detected_matched);

    % Calculate Se and P+
    if (TP + FN) > 0
        Se = TP / (TP + FN);
    else
        Se = NaN; % Avoid division by zero if no clean peaks
    end
    if (TP + FP) > 0
        P_plus = TP / (TP + FP);
    else
        P_plus = NaN; % Avoid division by zero if no detected peaks
    end
end

% Evaluate performance
tolerance_ms = 100; % Tolerance window in milliseconds (e.g., +/- 100ms around true peak)
tolerance_samples = round(tolerance_ms / 1000 * Fs);

fprintf('R-Peak Detection Tolerance: +/- %d samples (%.0f ms)\n', tolerance_samples, tolerance_ms);

% Static Filter Evaluation
[Se_static_mlii, P_plus_static_mlii, TP_s_mlii, FP_s_mlii, FN_s_mlii] = evaluate_rpeaks(locs_static_mlii, locs_clean_mlii, tolerance_samples);
[Se_static_v1, P_plus_static_v1, TP_s_v1, FP_s_v1, FN_s_v1] = evaluate_rpeaks(locs_static_v1, locs_clean_v1, tolerance_samples);

% LMS+Static Filter Evaluation
[Se_lms_mlii, P_plus_lms_mlii, TP_l_mlii, FP_l_mlii, FN_l_mlii] = evaluate_rpeaks(locs_lms_mlii, locs_clean_mlii, tolerance_samples);
[Se_lms_v1, P_plus_lms_v1, TP_l_v1, FP_l_v1, FN_l_v1] = evaluate_rpeaks(locs_lms_v1, locs_clean_v1, tolerance_samples);

fprintf('MLII - Static Filter R-Peaks: Se=%.2f%%, P+=%.2f%% (TP=%d, FP=%d, FN=%d / %d total clean)\n', Se_static_mlii*100, P_plus_static_mlii*100, TP_s_mlii, FP_s_mlii, FN_s_mlii, length(locs_clean_mlii));
fprintf('MLII - LMS+Static Filter R-Peaks: Se=%.2f%%, P+=%.2f%% (TP=%d, FP=%d, FN=%d / %d total clean)\n', Se_lms_mlii*100, P_plus_lms_mlii*100, TP_l_mlii, FP_l_mlii, FN_l_mlii, length(locs_clean_mlii));
fprintf('V1   - Static Filter R-Peaks: Se=%.2f%%, P+=%.2f%% (TP=%d, FP=%d, FN=%d / %d total clean)\n', Se_static_v1*100, P_plus_static_v1*100, TP_s_v1, FP_s_v1, FN_s_v1, length(locs_clean_v1));
fprintf('V1   - LMS+Static Filter R-Peaks: Se=%.2f%%, P+=%.2f%% (TP=%d, FP=%d, FN=%d / %d total clean)\n', Se_lms_v1*100, P_plus_lms_v1*100, TP_l_v1, FP_l_v1, FN_l_v1, length(locs_clean_v1));

% Optional: Plot detected peaks for visual verification (on a shorter segment)
figure('Name', 'R-Peak Detection Verification (MLII)');
plot(t(1:PlotEndSample), ecg_MLII(1:PlotEndSample), 'g'); hold on;
plot(t(locs_clean_mlii(locs_clean_mlii<=PlotEndSample)), pks_clean_mlii(locs_clean_mlii<=PlotEndSample), 'k^', 'MarkerFaceColor','k', 'MarkerSize', 8);
plot(t(1:PlotEndSample), filt_MLII_static(1:PlotEndSample), 'b');
plot(t(locs_static_mlii(locs_static_mlii<=PlotEndSample)), pks_static_mlii(locs_static_mlii<=PlotEndSample), 'bo', 'MarkerFaceColor','b');
plot(t(1:PlotEndSample), filt_MLII_lms_combined(1:PlotEndSample), 'c--');
plot(t(locs_lms_mlii(locs_lms_mlii<=PlotEndSample)), pks_lms_mlii(locs_lms_mlii<=PlotEndSample), 'co', 'MarkerFaceColor','c');
legend('Clean ECG', 'Clean R-Peaks', 'Static Filtered', 'Static Detected R-Peaks', 'LMS+Static Filtered', 'LMS Detected R-Peaks');
title('R-Peak Detection Example (MLII)');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
xlim([0 t(PlotEndSample)]);

fprintf('=== Analysis Complete ===\n');