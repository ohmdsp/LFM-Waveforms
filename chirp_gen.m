% Script:           chirp_gen.m
% Author:           drohm   
%
% LFM waveform generator for simulation use.

% B = bandwidth of chirp
% T = Time of chirp duration
%  
%==========================================================================
%==========================================================================
close all;clear all
B = 30.0e6;         % Set the bandwidth of the chirp
T = 42e-6;          % Set the time duration of the chirp

time_B_product = B * T;
if(time_B_product < 5 )
    fprintf('************ Time Bandwidth product is TOO SMALL ***************')
    fprintf('\n Change B and or T')
    return
else
end
% Compute alpha
mu = 2. * pi * B / T;
npow = nextpow2(5 * B * T + 1);
npoints = 1*(2^(npow));

% Determine sampling times
delt = linspace(-T/2., T/2., npoints); 
M = length(delt);

% Compute the complex LFM representation
Ichannal = cos(mu .* delt.^2 / 2.); % Real part
Qchannal = sin(mu .* delt.^2 / 2.); % Imaginary Part
LFM = Ichannal + sqrt(-1) .* Qchannal;  % complex signal

SNR = 15;
w_n = randn(1,M) + i*randn(1,M);
w = (10^(-SNR/10))*w_n;        % Additive white noise
LFM = w+LFM;

% Compute the FFT of the LFM waveform
LFMFFT = fftshift(fft(LFM));

% Plot the real and Inginary parts and the spectrum
sampling_interval = T / npoints;
freqlimit = 0.5 / sampling_interval;
freq = linspace(-freqlimit,freqlimit,npoints);

figure(1)
subplot(4,1,1)
plot(delt,Ichannal,'k');
axis([-T/2 T/2 -1.5 1.5])
grid
xlabel('Time - seconds')
ylabel('Units of Waveform')
title('Real part of an LFM waveform')

subplot(4,1,2)
plot(delt,Qchannal,'k');
axis([-T/2 T/2 -1.5 1.5])
grid
xlabel('Time - seconds')
ylabel('Units of Waveform')
title('Imaginary part of LFM waveform')

subplot(4,1,3)
plot(freq, abs(LFMFFT),'k');
grid
xlabel('Frequency - Hz')
ylabel('Amplitude spectrum')
title('Spectrum for an LFM waveform')

% Chirp detection method
% Computer delay and multiply of chirp with itself 
delay = 2000;
chrp_d = real(LFM(delay:end)) ;
chrp = real(LFM(1:length(chrp_d))); 
output_ch = abs(fftshift(fft(chrp_d.*chrp))).^2;
y = output_ch/max(output_ch);
fs = 1/sampling_interval;
f = (fs/length(y))*(-(length(y)/2):(length(y)/2)-1);

% Plot delay and multiple result
subplot(4,1,4)
plot(10e-6*f,10*log10(y),'k');
grid
xlabel('Frequency (MHz)')
ylabel('Normalized Amplitude (dB)')
title('Chirp Detection Tone Plot')
ylim([-75 0])

