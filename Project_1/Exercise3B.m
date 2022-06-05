 %Second Part
fs = 8000;
n = 0 : 100;
phi = 3*pi/4;

f0 = 100 : 125 : 475;
xn1 = {};

for i = 1 : length(f0)
    xn1{i} = sin(2*pi*(f0(i)/fs)*n + phi);
end    

Xn1 = {};
f_axis = -1/2 : 1/1024: 1/2 - 1/1024;

% figure();
% hold on
for i = 1 : length(f0)
    figure();
    Xn1{i} =(fft(xn1{i},1024));
    plot(f_axis, Xn1{i},'DisplayName',[num2str(f0(i)), ' Hz']);
    xlabel('Frequency in Hz');
    ylabel('Amplitude');
    title('Spectrum of 100-475Hz Signals Sampled at 8000Hz');
    legend show
    
end
% hold off

f0 = 7525 : 125 : 7900;
xn2 = {};
for i = 1 : length(f0)
    xn2{i} = sin(2*pi*(f0(i)/fs)*n + phi);
end    

Xn2 = {};
f_axis = -1/2 : 1/1024: 1/2 - 1/1024;

% figure();
% hold on
for i = 1 : length(f0)
    figure();
    Xn2{i} = (fft(xn2{i},1024));
    plot(f_axis, Xn2{i},'DisplayName',[num2str(f0(i)), ' Hz']);
    xlabel('Frequency in Hz');
    ylabel('Amplitude');
    title('Spectrum of 7525-7900Hz Signals Sampled at 8000Hz');
    legend show
end
% hold off
