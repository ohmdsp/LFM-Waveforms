clear

Bandwidth = 32e+6;
Carrier = 130e+6;

m = 1:1:200;
fs = ((Carrier - (Bandwidth/2))./m) + ((Carrier + (Bandwidth/2))./(m+1));

% Use this for spectrum centered at 1/4 the sampling rate
%m = 1:1:200;
%fs = (4*Carrier)./m;

%Use this to assure no negative images between DC and first positive image
%m = 2:2:200;
%fs = (2*Carrier - Bandwidth)./m;

index = 1;
NumInvalid = 0;
for k = 1:length(m)
   if (fs(k) > 2*Bandwidth)
      f(index) = fs(k);
      index = index + 1;
   else
      NumInvalid = NumInvalid + 1;
   end
end

% Number of sampling rates in our range of m that fall below 2*BW
if (index == 1)
    error('Could not compute any valid bandpass sampling rates for this signal')
    return;
end

[Fmin, MinIndex] = min(f);
sort(f);
plot_images(Bandwidth, Carrier, Fmin);
