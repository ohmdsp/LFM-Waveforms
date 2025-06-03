% Script:           chirp_gen.m
% Author:           D.R.Ohm   
% Rev:              1.2
% Date:             April 4, 2005
%
% Chirp generator and pulse compression for range detection.

% B = bandwidth of chirp
% T = Time of chirp duration
%
%==========================================================================
%==========================================================================
close all;clear all
B = 30.0e6;
T = 42e-6;

time_B_product = B * T;
if(time_B_product < 5 )
    fprintf('************ Time Bandwidth product is TOO SMALL ***************')
    fprintf('\n Change B and or T')
    return
else
end

%-Compute alpha (chirp rate)
mu = 2. * pi * B / T;
npow = nextpow2(5 * B * T + 1);
npoints = 1*(2^(npow));

%-Determine sampling times
delt = linspace(0, T, npoints); % 
M = length(delt);

%-Compute the complex LFM representation
Ichannel = cos(mu .* delt.^2 / 2.); % Real part
Qchannel = sin(mu .* delt.^2 / 2.); % Imaginary Part
LFM = Ichannel + sqrt(-1) .* Qchannel;  % complex signal

%-Add WGN Noise
SNR = 15;
w_n = randn(1,M) + 1i*randn(1,M);
w = (10^(-SNR/10))*w_n;        % Additive white noise
LFM = w+LFM;

%-Compute the FFT of the LFM waveform
LFMFFT = fftshift(fft(LFM));
sampling_interval = T / npoints;
freqlimit = 0.5 / sampling_interval;
freq = linspace(-freqlimit,freqlimit,npoints);

%-Plot the real, imaginary parts and the spectrum and magnitude
figure(1)
subplot(5,1,1)
plot(delt,Ichannel,'k');
axis([0 T -1.5 1.5])
grid
xlabel('Time - seconds')
ylabel('Units of Waveform')
title('Real part of an LFM waveform')

subplot(5,1,2)
plot(delt,Qchannel,'k');
axis([0 T -1.5 1.5])
grid
xlabel('Time - seconds')
ylabel('Units of Waveform')
title('Imaginary part of LFM waveform')

%-Plot the reversed real part
subplot(5,1,3)
plot(delt,flip(Ichannel),'k');
axis([0 T -1.5 1.5]); grid
xlabel('Time - seconds')
ylabel('Units of Waveform')
title('Real part of an LFM waveform (Reversed)')

subplot(5,1,4)
plot(freq.*1/1e6, abs(LFMFFT),'k'); grid
xlabel('Frequency - Mz')
ylabel('Magnitude spectrum')
title('Spectrum for an LFM waveform')

%-Add delay (in samples) to return chirp
delay = 3000;
tau = delay*sampling_interval;
disp(['Delay of return is = ' num2str(tau*1/1e-6) ' microseconds']);
chirp_rx = [zeros(1,delay) LFM(1:end)] ;  % create delayed (truncated) chirp
chirp_tx = [LFM(1:end) zeros(1,delay)];   % truncate tx ref chirp

%-Compute pulse compression (aka matched filter) and then plot spectrum again 
Nfft = length(chirp_tx);
%out_comp = fft(conj(flip(chirp_tx)), Nfft) .* fft(chirp_rx, Nfft);
out_comp = fft(conj(flip(chirp_tx)), Nfft) .* fft(chirp_rx, Nfft);
out_comp_ifft = ifft(out_comp, Nfft);
out_comp_mag = abs(out_comp_ifft).^2;
y = out_comp_mag/max(out_comp_mag);
delty = sampling_interval.*(1:length(y));

%-Plot output of matched filter (pulse compression)
subplot(5,1,5)
plot(delty.*1/1e-6,10*log10(y),'k');
xlabel('Delay (us)')
ylabel('Normalized Amplitude (dB)')
title('LFM Range Delay')
axis([0 T*1/1e-6 -75 10]); grid
