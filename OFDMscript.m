%% Project Script for ECE3141 Friday 12pm D- OFDM 
% Hamish Cranston & Reece Harrison


%% Setup
clc;clear all;close all;


%% Parameters
N = 64; % [] number of subcarriers
T = 1e-3; % [s] symbol period
M = 16; %[] QAM in use
fs = N/T; %[Hz] sampling frequency
df = 1/T; %[Hz] subcarrier spacing
t = (0:N-1)/fs; %[s] time vector 


%% Generating N random subcarriers (Xn) (our datastream)
% This is simulated QAM-4, so 2 bits and 4 possible symbols. Adjust to
% center. No noise at the moment.

%This will only generate an individual qam signal. Could simulate
%parallelising serial data.

input_data = randi([0 M-1],N,1);

% Encode QAM
qamEncoded = qammod(input_data,M);



figure;
grid on
plot(qamEncoded,'*')
xlim([-5,5]);
ylim([-5,5])
xlabel("Real Axis A");
title("QAM encoded input")
ylabel("Imaginary Axis B")
grid on;

%% Inverse Forward Fourier Transformation (convert to time domain)
signal = ifft(qamEncoded);

% Amplitudes
signalAmplitudes = abs(signal);

%plot the actual signal (could add real+imag parts?)
figure;
plot(t*1e3,signalAmplitudes)
title("Magnitude of Time domain OFDM encoded signal");
xlabel("Time (ms)");
ylabel("Amplitude (arb)");


% Signal Amplitude Histogram
figure;
histogram(signalAmplitudes, 'Normalization', 'pdf');
xlabel('Amplitude');
ylabel('Probability Density');
title('Histogram of Signal Amplitudes');

%% Frequency Domain Spectrum (PSD)

% We use a larger FFT size
N_fft = N * 2; 
signal_freq = fft(signal, N_fft);

% Shift the spectrum so 0 Hz is in the center
signal_shifted = fftshift(signal_freq);

% Db Power spectrum we use 20*log10 for amplitude-to-dB conversion
PSD = 20*log10(abs(signal_shifted));
PSD_norm = PSD - max(PSD); % Normalize so peak is at 0 dB

% The range is from -fs/2 to fs/2
f_axis = linspace(-fs/2, fs/2, N_fft);

figure;
plot(f_axis/1e6, PSD_norm, 'LineWidth', 1.5);
grid on;
title('Frequency Domain Spectrum of Sent OFDM Signal');
xlabel('Frequency (MHz) (centered relative to carrier)'); % carrier frequency is at 0
ylabel('Relative Power (dB)');
ylim([-60 5]); % Adjust to see the noise floor/sidelobes

%% Calculating Peak-to-Average Power Ratio and Performance Statistics
PAPR = 10*log10(max(signalAmplitudes.^2)./mean(signalAmplitudes.^2));

fprintf("Peak-to-Average power ratio found to be %.2f dB\n",PAPR)

fprintf("\nFrequency Domain Performance\n")
fprintf("Subcarrier spacing %.0f Hz\n",df)
fprintf("Total Bandwidth %.2f MHz\n",N*df/1e6);
fprintf("Number of subcarriers %d\n",N)

fprintf("\nTime Domain Performance\n")
fprintf("Symbol duration (1 symbol) %.0f ms\n",T*1e3)
fprintf("Symbol rate (1 symbol) %.0f symbols/s\n",1/T)
fprintf("Max Data Rate without OFDM using %.0f-QAM (send one QAM symbol in one symbol period) %.2f kbps ((QAM) log2(M) * 1/T)\n",M,log2(M)/(T*1e3))
fprintf("Max Data Rate with OFDM using %.0f-QAM %.2f Mbps ((QAM) log2(M) * N/T)\n",M,N*log2(M)/(T*1e6))


%% Add noise to the transmitted signal.

noisy_signal = awgn(signal,60);

%% Forward Fourier Transform recovers Xn (our datastream)
recoveredbitPattern = fft(noisy_signal);

figure()
hold on
plot(recoveredbitPattern,'*')
plot(qamEncoded,'^')
legend("Recovered Pattern","Original Pattern")
title("Recovered Bit Pattern")
xlabel("Real Axis");
ylabel("Imaginary Axis")
grid on;

%This is now a relative tolerance
check=isapprox(qamEncoded,recoveredbitPattern,RelativeTolerance=0.2);
equalcheck = all(check(:)==1);
if equalcheck==1
    fprintf('\nThe original and recovered bit patterns are equal.\n');
else
    fprintf('\nThe original and recovered bit patterns are NOT equal.\n');
end

%% Constructing a signal across multiple transmission periods without cyclic prefix:

num_periods = 10;
pause_length = N/4;

pause_vec = zeros([1,pause_length]);

overall_signal = [];

signal_zero = zeros([1, N]);

cyclic_prefixes = [];

time_overall = (0:(N+pause_length)*num_periods-1)/fs; %[s] time vector 


for i = 1:num_periods

    input_data = randi([0 M-1],N,1);

% Encode QAM
    qamEncoded = qammod(input_data,M);

    signal = ifft(qamEncoded);

    prefix = signal((N-pause_length):N-1);

    cyclic_prefixes = [cyclic_prefixes,prefix',signal_zero];

    overall_signal = [overall_signal, pause_vec,signal'];


end

overall_abs = abs(overall_signal);

prefix_abs = abs(cyclic_prefixes);

combo_abs = overall_abs + prefix_abs;

PAPR = 10*log10(max(combo_abs.^2)./mean(combo_abs.^2));
fprintf("Peak-to-Average power ratio found to be %.2f dB with cyclic prefix\n",PAPR)

figure();
hold on
plot(time_overall*1e3,overall_abs)
plot(time_overall*1e3,prefix_abs)

legend("OFDM Time Domain Signal","Cyclic Prefix");

title("A ten sequence OFDM transmission with cyclic prefixes");
xlabel("Time (ms)")
ylabel("Amplitude (arb)")

%% Plot orthogonal carriers for an example.

fs = 1024;
N_fft = 8192;
time = (0:fs-1)/fs;
multipliers = (1:10) + 300 ;
multipliers = multipliers';

orthogonal_signals = cos(2*pi*multipliers*time);
orthogonal_signals = [zeros([10,fs*2]),orthogonal_signals,zeros([10,fs*2])];

fft_signals = fft(orthogonal_signals',N_fft);

freqs = linspace(0,fs,N_fft);

half_point = floor(N_fft/2);

figure()
plot(freqs(1:half_point),abs(fft_signals(1:half_point,:)))
xlim([295,316])
title("FFT of Orthogonal Frequencies")
xlabel("Frequency (Hz)");
ylabel("Amplitude (arb)");
