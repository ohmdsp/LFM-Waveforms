% Script:           chirp_samp.m
% Author:           D.R.Ohm   
% Rev:              1.2
% Date:             April 4, 2005
%
% Explores bandpas sampling for a LFM waveform (chirp signal).

% B = bandwidth of chirp
% T = Time of chirp duration
%
%==========================================================================
%==========================================================================
clear all;close all

fc = 130e6;
phi = 0;

fs = 105e6;
%fs = 75e6;
Ts = 1/fs;
T = 42e-6;
npoints = T*fs;
tp = linspace(-T/2, T/2, npoints);
B = 32e6;   % 32 MHz bandwidth
mu = 2. * pi * B / T;
M = length(tp);

%-Compute the complex LFM representation
Ichannal = cos(2*pi*fc.*tp + mu .* tp.^2 / 2.); % Real part
Qchannal = sin(2*pi*fc.*tp + mu .* tp.^2 / 2.); % Imaginary Part
LFM = Ichannal + sqrt(-1) .* Qchannal;  % complex valued signal

SNR = 10;
w_n = randn(1,M) + 1i*randn(1,M);
w = (10^(-SNR/10))*w_n;        % Additive white noise for defined SNR
LFM = LFM+w;
%LFMFFT = fftshift(fft(LFM));
sig = real(LFM);

%-Compute and plot spectrum of LFM waveform sampled at > BW
pow = nextpow2(npoints);
fpoints = 2^pow;
sigfft = fftshift(fft(sig,fpoints));
freqlimit = .5*fs;
freq = linspace(-freqlimit,freqlimit,fpoints);
spec = abs(sigfft);
subplot(2,1,1)
plot(freq,spec,'k')
axis tight

%-Show spectrum replicas due to sampling in time-domain
sigfft_long = (fft(sig,fpoints));
specl = abs(sigfft_long);
spec_long = fftshift([specl specl specl specl]);
freq_long = linspace(-2*freqlimit,2*freqlimit,4*fpoints);
subplot(2,1,2)
plot(freq_long,spec_long,'k')
axis tight
