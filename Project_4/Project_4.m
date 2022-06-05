close all;
clear;

%% 1
Wc=0.4*pi;
Fs=100;
N=21;
n=N+1;
Fc=Wc/(2*pi);
Fn=Fc/Fs;
Wn=2*pi*Fn;

B1=fir1(N,Wn,rectwin(n));
[H1, w1] = freqz(B1);
B2=fir1(N,Wn,hamming(n));
[H2, w2] = freqz(B2);
figure(1);
plot(w1, abs(H1),'r');
hold on;
plot(w2, abs(H2),'b');
ylabel('frequency response amplitude');
xlabel('Frequency');
legend('rectangular window','hamming window');
title('Rectangular and Hammming window with same length');
hold off;

%% 2.1
%%%%%%%%%%%%%%%%%%%%%%%%%%
 %close all;

%initialize values and normalize frequency
Wc=0.5*pi;
Fs=100;
Fc=Wc/(2*pi);
Fn=Fc/Fs;
Wn=2*pi*Fn;

%initializing the lengths of the windows
N1=21;
N2=41;
n1=N1+1;
n2=N2+1;

%hamming windows made and designed
Bham1=fir1(N1,Wn,hamming(n1));
[H3, w3] = freqz(Bham1);
Bham2=fir1(N2,Wn,hamming(n2));
[H4, w4] = freqz(Bham2);
figure(2);
plot(w3,abs(H3),'k');
hold on;
plot(w4,abs(H4),'g');
legend('N=21','N=41');
ylabel('frequency response amplitude');
xlabel('Frequency');
title('Hamming windows');
hold off;

%hanning windows made and designed
figure(3);
Bhan1=fir1(N1,Wn,hanning(n1));
[H5, w5] = freqz(Bhan1);
Bhan2=fir1(N2,Wn,hanning(n2));
[H6, w6] = freqz(Bhan2);
plot(w5,abs(H5),'k');
hold on;
plot(w6,abs(H6),'g');
legend('N=21','N=41');
ylabel('frequency response amplitude');
xlabel('Frequency');
title('Hanning windows');
hold off;

%% 2.2
%making the axis of signal x(t)
% close all;
t = [-4:1/Fs:4];
N = 2048;
xt = sin(15*t) + 0.25*sin(200*t);
xf = fftshift(fft(xt,N)*(1/Fs));
F_axis = [-Fs/2:Fs/N:Fs/2 - Fs/N];


%%%%%%%%%%%%%%%%%%%%%%%%%
%filter the signal by hamming window with N=21

Ytham1 = filter2(Bham1,xt); 
Yfham1 = fftshift(fft(Ytham1,N)*(1/Fs));

figure(4);
subplot(2,1,1);
plot(F_axis, abs(xf));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('Fourier transform of the signal');
subplot(2,1,2);
plot(F_axis, abs(Yfham1));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('filtered signal by hamming 21');

%%%%%%%%%%%%%%%%%%%%%%%%%
%filter the signal by hamming window with N=41

Ytham2 = filter2(Bham2,xt); 
Yfham2 = fftshift(fft(Ytham2,N)*(1/Fs));

figure(5);
subplot(2,1,1);
plot(F_axis, abs(xf));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('Fourier transform of the signal');
subplot(2,1,2);
plot(F_axis, abs(Yfham2));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('filtered signal by hamming 41');

%%%%%%%%%%%%%%%%%%%%%%%%%
%filter the signal by hanning window with N=21

Ythan1 = filter2(Bhan1,xt); 
Yfhan1 = fftshift(fft(Ythan1,N)*(1/Fs));

figure(6);
subplot(2,1,1);
plot(F_axis, abs(xf));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('Fourier transform of the signal');
subplot(2,1,2);
plot(F_axis, abs(Yfhan1));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('filtered signal by hanning 21');

%%%%%%%%%%%%%%%%%%%%%%%%%
%filter the signal by hanning window with N=41


Ythan2 = filter2(Bhan2,xt);
Yfhan2 = fftshift(fft(Ythan2,N)*(1/Fs));

figure(7);
subplot(2,1,1);
plot(F_axis, abs(xf));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('Fourier transform of the signal');
subplot(2,1,2);
plot(F_axis, abs(Yfhan2));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('filtered signal by hanning 41');


%% 2.3

%initialize values and normalize frequency
Wc=0.5*pi;
Fs=50;
Fc=Wc/(2*pi);
Fn=Fc/Fs;
Wn=2*pi*Fn;

%initializing the lengths of the windows
N1=21;
N2=41;
n1=N1+1;
n2=N2+1;

%hamming and hanning windows made and designed
Bham1=fir1(N1,Wn,hamming(n1));
Bham2=fir1(N2,Wn,hamming(n2));

Bhan1=fir1(N1,Wn,hanning(n1));
Bhan2=fir1(N2,Wn,hanning(n2));

%making the axis of signal x(t)
t = [-4:1/Fs:4];
N = 2048;
xt = sin(15*t) + 0.25*sin(200*t);
xf = fftshift(fft(xt,N)*(1/Fs));
F_axis = [-Fs/2:Fs/N:Fs/2 - Fs/N];


%%%%%%%%%%%%%%%%%%%%%%%%%
%filter the signal by hamming window with N=21

Ytham1 = filter2(Bham1,xt); 
Yfham1 = fftshift(fft(Ytham1,N)*(1/Fs));

figure(8);
subplot(2,1,1);
plot(F_axis, abs(xf));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('Fourier transform of the signal');
subplot(2,1,2);
plot(F_axis, abs(Yfham1));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('filtered signal by hamming 21');

%%%%%%%%%%%%%%%%%%%%%%%%%
%filter the signal by hamming window with N=41

Ytham2 = filter2(Bham2,xt); 
Yfham2 = fftshift(fft(Ytham2,N)*(1/Fs));

figure(9);
subplot(2,1,1);
plot(F_axis, abs(xf));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('Fourier transform of the signal');
subplot(2,1,2);
plot(F_axis, abs(Yfham2));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('filtered signal by hamming 41');

%%%%%%%%%%%%%%%%%%%%%%%%%
%filter the signal by hanning window with N=21

Ythan1 = filter2(Bhan1,xt); 
Yfhan1 = fftshift(fft(Ythan1,N)*(1/Fs));

figure(10);
subplot(2,1,1);
plot(F_axis, abs(xf));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('Fourier transform of the signal');
subplot(2,1,2);
plot(F_axis, abs(Yfhan1));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('filtered signal by hanning 21');

%%%%%%%%%%%%%%%%%%%%%%%%%
%filter the signal by hanning window with N=41


Ythan2 = filter2(Bhan2,xt);
Yfhan2 = fftshift(fft(Ythan2,N)*(1/Fs));

figure(11);
subplot(2,1,1);
plot(F_axis, abs(xf));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('Fourier transform of the signal');
subplot(2,1,2);
plot(F_axis, abs(Yfhan2));
ylabel('frequency response amplitude');
xlabel('Frequency');
title('filtered signal by hanning 41');