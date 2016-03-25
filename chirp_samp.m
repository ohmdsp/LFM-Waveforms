% Script:       chrip_samp.m 
% Author:       drohm
%
% Simulate aliasing from bandpass sampling a chirp with specified setup
% parameters.
%
% fc = carrier frequency
% B = bandwidth of chirp
% T = time duration of chirp
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
clear all;close all

fc = 130e6;     % carrier frequency
phi = 0;        % initial phase
fs = 105e6;     % sample rate (subsampling)
Ts = 1/fs;
T = 42e-6;      % chirp time duration
npoints = T*fs;
tp = linspace(-T/2, T/2, npoints);
B = 32e6;       % 32 MHz chirp bandwidth
mu = 2. * pi * B / T;
M = length(tp);

% Compute the complex LFM representation
Ichannal = cos(2*pi*fc.*tp + mu .* tp.^2 / 2.); % Real part
Qchannal = sin(2*pi*fc.*tp + mu .* tp.^2 / 2.); % Imaginary Part
LFM = Ichannal + sqrt(-1) .* Qchannal;  % complex signal

% Add Noise
SNR = 10;
w_n = randn(1,M) + i*randn(1,M);
w = (10^(-SNR/10))*w_n;        % Additive white noise for defined SNR
LFM = LFM+w;
sig = real(LFM);

% Computer FFT
pow = nextpow2(npoints);
fpoints = 2^pow;
sigfft = fftshift(fft(sig,fpoints));
freqlimit = .5*fs;
freq = linspace(-freqlimit,freqlimit,fpoints);
spec = abs(sigfft);
subplot(2,1,1)
plot(freq,spec,'k')
xlabel('Frequency (Hz)');ylabel('Magnitude')
axis tight

sigfft_long = (fft(sig,fpoints));
specl = abs(sigfft_long);
spec_long = fftshift([specl specl specl specl]);
freq_long = linspace(-2*freqlimit,2*freqlimit,4*fpoints);
subplot(2,1,2)
plot(freq_long,spec_long,'k')
xlabel('Frequency (Hz)');ylabel('Magnitude')
axis tight
