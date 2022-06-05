%Exercise 3A

%The Initial Function
%x(t) = 10*cos(2*pi*20*t)-4*sin(2*pi*40*t+5);

sampleNum  = 128;

%Sampling at 240Hz>2*40 No Aliasing  
f = 40;
T_s = 1/f;

%Setup time axis
t = 0 : 1/10000 : (sampleNum-1) * T_s; %low step approximates continues signal for comparison
ts = 0 : T_s : (sampleNum-1) * T_s; %Vector ends at (sampleNum-1) * T_s in order to take sampleNum samples

xt = 10*cos(2*pi*20*t)-4*sin(2*pi*40*t+5);
xst = 10*cos(2*pi*20*ts)-4*sin(2*pi*40*ts+5);
%Copy at at T_s
xst2 = 10*cos(2*pi*(20-1/f)*ts)-4*sin(2*pi*(40-1/f)*ts+5);

%DTF to plot the frequency spectrum
FFTSamplingFreq  = 1024;
Xs = fftshift(fft(xst,FFTSamplingFreq));
f_axis = -1/2 : 1/FFTSamplingFreq: 1/2 - 1/FFTSamplingFreq;

Xs2 = (fft(xst2,FFTSamplingFreq));
Xsum = Xs + Xs2;

% %Plot signals in time domain
% figure();
% hold on
% plot(ts,xst,'DisplayName','Sampled signal');
% plot(t,xt,'DisplayName','Original Signal');
% 
% title('Comparison of Original and Sampled Signals');
% xlabel('Time');
% ylabel('Amplitude');
% legend show
% hold off;


%Plot Signals in Frequency Domain
%plot copies to confirm no aliasing at frequency domain
figure();
plot(f_axis,Xs);
title('Spectrum of Sampled Signal');
xlabel('Frequency in Hz');
ylabel('Amplitude');

figure();
plot(f_axis,Xs2);
title('Spectrum Copy T_s');
xlabel('Frequency in Hz');
ylabel('Amplitude');

%plot sum of copies
figure();
plot(f_axis,Xsum);
title('Summed Spectrum of Sampled Signal and Copy of it at T_s');
xlabel('Frequency in Hz');
ylabel('Amplitude');

