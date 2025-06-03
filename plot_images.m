function plot_images(bw, fc, fs)

for k=-200:2:200,
    SpectrumLocations(200+k+1) = ((k/2)*fs) - fc;
    SpectrumLocations(200+k+2) = ((k/2)*fs) + fc;
end

% Want about 1000 points across BW
incr = bw/1000;
faxis = (-2*fc):incr:(2*fc);
NewSpectrumNeg = zeros(1, length(faxis));
NewSpectrumPos = zeros(1, length(faxis));
s = (-bw/2):incr:(bw/2);
OrigSpectrum = tripuls(s, bw);
% Generate all the images from the original negative spectral image
for k=1:2:length(SpectrumLocations),
    Location = SpectrumLocations(k) - (bw/2);
    index = floor((2*fc + Location)/incr) + 1;
    if ((index > 0)&&(index <= (length(NewSpectrumNeg)-length(OrigSpectrum)))),
        NewSpectrumNeg(index:(index+length(OrigSpectrum))-1) = OrigSpectrum;
    end
end
% Generate all the images from the original positive spectral image
for k=2:2:length(SpectrumLocations),
    Location = SpectrumLocations(k) - (bw/2);
    index = floor((2*fc + Location)/incr) + 1;
    if ((index > 0)&&(index <= (length(NewSpectrumPos)-length(OrigSpectrum)))),
        NewSpectrumPos(index:(index+length(OrigSpectrum))-1) = OrigSpectrum;
    end
end
figure
plot(faxis*1e-6, NewSpectrumNeg, 'b-', faxis*1e-6, NewSpectrumPos, 'r-');
PlotTitle = sprintf('Red=PosImages, Blue=NegImages, Fs=%1.2fMHz, BW=%dMHz, Carrier=%dMHz)', fs*1e-6, bw*1e-6, fc*1e-6);
title (PlotTitle);
xlabel('Frequency (MHz)');
