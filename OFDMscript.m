%% Project Script for ECE3141 Friday 12pm D- OFDM 
% Hamish Cranston & Reece Harrison


%% Setup
clc;clear all;close all;


%% Parameters
N = 8096; % [] number of subcarriers
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