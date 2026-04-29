%% Project Script for ECE3141 Friday 12pm D- OFDM 
% Hamish Cranston & Reece Harrison


%% Setup
clc;clear all;close all;


%% Parameters
N = 200; % [] number of subcarriers
T = 1; % [ms] symbol period


%% Generating N random subcarriers (Xn) (our datastream)
bitPattern=randi([0 1], 1, N)+ j*randi([0 1], 1, N);

%% Inverse Forward Fourier Transformation (convert to time domain)
signal = ifft(bitPattern);

% Amplitudes
signalAmplitudes = abs(signal);

% Signal Amplitude Histogram
figure;
histogram(signalAmplitudes, 'Normalization', 'pdf');
xlabel('Amplitude');
ylabel('Probability Density');
title('Histogram of Signal Amplitudes');

%% Calculating Peak-to-Average Power Ratio
PAPR = 10*log10(max(signalAmplitudes.^2)./mean(signalAmplitudes.^2));

fprintf("Peak-to-Average power ratio found to be %.2f dB\n",PAPR)

%% Forward Fourier Transform recovers Xn (our datastream)
recoveredbitPattern = fft(signal);

check=isapprox(bitPattern,recoveredbitPattern,"verytight");
equalcheck = all(check(:)==1);
if equalcheck==1
    fprintf('The original and recovered bit patterns are equal.\n');
else
    fprintf('The original and recovered bit patterns are NOT equal.\n');
end