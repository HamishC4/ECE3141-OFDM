%% Project Script for ECE3141 Friday 12pm D- OFDM 
% Hamish Cranston & Reece Harrison


%% Setup
clc;clear all;close all;


%% Parameters
N = 128; % [] number of subcarriers
T = 1e-3; % [s] symbol period
M = 16; %[] QAM in use
fs = N/T; %[Hz] sampling frequency
df = 1/T; %[Hz] subcarrier spacing
t = (0:N-1)/fs; %[s] time vector 
fs2 = 100* N/T;
t2 = 0:1/fs2 : T; %[s] time vector 
snr = 10; % SNR for AWGN, in dB
%% Generating N random subcarriers (Xn) (our datastream)
% This is simulated QAM-4, so 2 bits and 4 possible symbols. Adjust to
% center. No noise at the moment.

%This will only generate an individual qam signal. Could simulate
%parallelising serial data.

input_data = randi([0 M-1],N,1);

% Use the bit encoding later:
input_data_bin = int2bit(input_data,log2(M));
% Encode QAM

qamEncoded = qammod(input_data,M);

%% This section visualises. Only for small dimensions of data, irrelevant for larger simulations.
if (N <= 32)
    subcarriers = zeros(N,length(t2));
    
    %Plot each subcarrier
    figure();
    hold on;
    title("Visualisation of Individual Carrier Frequencies ");
    xlabel("Time (s)");
    ylabel("Amplitude (arb)")
    for k = 1:N
        f_k = (k-1-N/2)*df;
        subcarriers(k,:) = qamEncoded(k)* exp(1i * 2 * pi * f_k * t2);
        plot(t2,real(subcarriers(k,:))  + k*5,LineWidth=4)    
    end
    % Plot the sum of Signals
    sum_subcarriers = sum(subcarriers);
    figure()
    plot(t2,sum_subcarriers,LineWidth=5)
    title("Sum of Subcarriers (IFFT)");
    xlabel("Time (s)");
    ylabel("Amplitude (arb)");
    
    
    % Plot the cyclic prefix: Truncate.
    t_cyclic = [t2,T:1/fs2 : 5*T/4];
    sig_cyclic = sum_subcarriers((end - ceil(length(t2)/4)) : end -1);
    
    
    figure();
    hold on
    plot(t_cyclic(length(sig_cyclic)+1 : end),sum_subcarriers,LineWidth=5);
    plot(t_cyclic(1:length(sig_cyclic)),sig_cyclic,LineWidth=5);
    title("Sum of Subcarriers with cyclic prefix in time domain");
    xlabel("Time (s)");
    ylabel("Amplitude (arb)");
    legend("OFDM symbol","Cyclic Prefix");
end




%%

figure;
grid on
plot(qamEncoded,'^',MarkerSize=10,LineWidth=10,Color=[1 0.5 0])
xlim([-5,5]);
ylim([-5,5])
xlabel("Real Axis A");
title("QAM encoded input")
ylabel("Imaginary Axis B")
axis square;

%% Inverse Forward Fourier Transformation (convert to time domain)
signal = ifft(ifftshift(qamEncoded));

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
stem(f_axis/1e6, PSD_norm, 'LineWidth', 1.5);
axis square
% grid on;
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
noisy_signal = awgn(signal,snr,'measured');
%% Forward Fourier Transform recovers Xn (our datastream)
recoveredbitPattern = fftshift(fft(noisy_signal));

figure()
hold on
axis equal
plot(recoveredbitPattern,'*',LineWidth=10,MarkerSize=10)
plot(qamEncoded,'^',MarkerSize=10,LineWidth=1)
xlim([-2,2]);
ylim([-2,2]);
legend("Recovered Pattern","Original Pattern")
title("Recovered Bit Pattern")
xlabel("Real Axis");
ylabel("Imaginary Axis")


%This is now a relative tolerance
check=isapprox(qamEncoded,recoveredbitPattern,RelativeTolerance=0.2);
equalcheck = all(check(:)==1);
if equalcheck==1
    fprintf('\nThe original and recovered bit patterns are equal.\n');
else
    fprintf('\nThe original and recovered bit patterns are NOT equal.\n');
end

%% Demodulating QAM signal to determine BER

dataOutput = qamdemod(recoveredbitPattern,M);

dataOutput_bin = int2bit(dataOutput,log2(M));

check  = dataOutput_bin ~= input_data_bin;

check_symbol = dataOutput ~= input_data;

ber = sum(check)/length(dataOutput_bin);

qer = sum(check_symbol)/length(dataOutput);

fprintf("\nThe bit error rate with %.0f dB SNR is %f \n",snr,ber);

fprintf("\nThe number of mismatched QAM-%.0f symbols out of %.0f is %.0f\n",M,N,sum(check_symbol));

fprintf("\nThe error rate of QAM-%.0f symbols is %f \n",M,qer);

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

%% Plot the PAPR of a large range of subcarriers.

N = unique(round(logspace(1,log10(1000000),1000)));
M =256; %[] QAM in use

snr = 1000; % SNR for AWGN, in dB

paprs_array = zeros([1,length(N)]);

ser_array = zeros([1,length(N)]);

ber_array = zeros([1,length(N)]);

for i = 1:length(N)

input_data = randi([0 M-1],N(i),1);
input_data_bin = int2bit(input_data,log2(M));
qamEncoded = qammod(input_data,M);
signal = ifft(ifftshift(qamEncoded));
noisy_signal = signal;

paprs_array(i) = 10*log10(max(abs(noisy_signal).^2)./mean(abs(noisy_signal).^2));

recoveredbitPattern = fftshift(fft(noisy_signal));


dataOutput = qamdemod(recoveredbitPattern,M);


dataOutput_bin = int2bit(dataOutput,log2(M));

check  = dataOutput_bin ~= input_data_bin;

check_symbol = dataOutput ~= input_data;

ber_array(i) = sum(check)/length(dataOutput_bin);

ser_array(i) = sum(check_symbol)/length(dataOutput);



end

%% 

figure();
semilogx(N,paprs_array)
hold on
grid on
title("PAPR versus number of subcarriers 256-QAM")
xlabel("Number of subcarriers N");
ylabel("PAPR (dB)");



%% Timing Jitter Simulation - find symbol error rate due to jitter in timing mismatch between clocks.

% Matlab example code used here for understanding how to make OFDM signal, from mathworks website.

M = 16; % Modulation order for 16QAM
nfft  = 64;
cplen = 16;
nSym  = 100;
nullIdx  = [1:6 33 64-4:64]'; % Pilot and guard carriers used to help prevent leakage between channels/track offsets.
numDataCarrs = nfft-length(nullIdx);

jitter_stds = linspace(0.001,7,100); % Standard deviation list

jittered_sers = zeros([10,length(jitter_stds)]); % Returned SERs.

%Vector of jitters.
for i  = 1:length(jitter_stds)

    % Repeat 10 times so we can average data.
    for j = 1:10
        inSym = randi([0 M-1],numDataCarrs,nSym);
        
        %Modulate QAM/OFDM with varying cyclic prefix
        qamSig = qammod(inSym,M,UnitAveragePower=true);
        outSig = ofdmmod(qamSig,nfft,cplen,nullIdx);
        
        % I did not use one of matlabs inbuilt channels due to complexity. Instead
        % I only simulated the timing jitter that might occur over the channel 
        t_ideal = 1:length(outSig);
        
        jitterSTD = jitter_stds(i); % Flatten this out
        
        %Jitter with Standard deviation. This simulates timing jitters in
        %the receiver clock relative to the transmitter clock.
        timing_noise = jitterSTD*randn(size(t_ideal));

        
        % Add jitter
        t_jittered = t_ideal + timing_noise;

        % Ensure null indices from jitter are filtered.
        t_jittered = max(min(t_jittered,length(outSig)),1);
        
        % Interpolate the signal to get jitter.
        outSig_with_jitter = interp1(t_ideal,outSig,t_jittered,'spline');
        
        outSig_with_jitter = outSig_with_jitter(:);
        
        %Demodulate the signal.
        rxSig = ofdmdemod(outSig_with_jitter,nfft,cplen,8,nullIdx);
        rxData = qamdemod(rxSig,M,UnitAveragePower=true);
        %Find SER
        num_se = sum(rxData(:)~=inSym(:));
        jittered_sers(j,i) = num_se / numel(inSym);
    end
    

end

% Compute average from matrix
average_jitter_ser = mean(jittered_sers,1);

%%
% Plotting
figure();

plot(jitter_stds,average_jitter_ser)
title("Symbol Error Rate for Randomly Distributed Timing Jitter")
xlabel("Standard Deviation of Jitter from Timing Period")
ylabel("Symbol Error Rate (SER)")
grid on





