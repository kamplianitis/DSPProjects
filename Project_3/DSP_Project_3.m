%askisi 1
close all;
clear all;

%sampling frequency 
fs=10*10^3;

%declare the values needed to create the butterworth filter

Wp = 2*pi*3*10^3; % we use the 2pi multi cause the command works in radians
Ws = 2*pi*4*10^3; %same here
Rp = 3; %passband riple
Rs = 50; %attenuation

[n,Wn] = buttord(Wp,Ws,Rp,Rs,'s');  % Gives minimum order of filter
%we use the s here in order to specify that our filter is analog.

%print the n to get the class of the system
disp(n)
%we will use the Wn afterwards cause it has the cut-off frequency of the
%analog signal

% we create Butterworth analog lowpass filter prototype.
[Z,P,K] = buttap(n); % this command returns the zeroes, the poles and the gain
%a lowpass filter prototype depending on its order (we will have n poles cause
% the system's order defines to us the number of the poles)

%we are still in a prototype filter
[NUM,DEN] = zp2tf(Z,P,K); % creation of a transfer func based on the poles the zeroes and the gain

% we create an analog lowpass filter that has a cutoff freq to Wn which
% is returned from the buttord that we used previously
[NUMT,DENT] = lp2lp(NUM,DEN,Wn); 


Nf=2048; %samples in the frequency zone
%the frequency-zone that we will use
f_zone =  [0:fs/(2*Nf):fs/2-fs/(2*Nf)];
%display the laplace transform given the num and the dnum we got from the
%function lp2lp that gives an analog filter. The frequency zone is
%multiplied with 2pi cause the frequency vector that is asked is in rad/s
hs = freqs(NUMT,DENT, 2*pi*f_zone);

freqs(NUMT,DENT, 2*pi*f_zone)

figure


[numd,dend] = bilinear(NUMT,DENT,fs);

[hz,w2] = freqz(numd,dend,Nf);
freqz(numd,dend,Nf)
figure


hold on
plot(f_zone, abs(hs))
plot(f_zone,abs(hz));
hold off;
title('Butterworth filters (digital and analog)');
xlabel('Frequency Zone (Hz)');
ylabel('Fequency Response of the filter');
legend('analog', 'digital');


figure 

plot(f_zone,abs(hs));
hold on;
semilogy(f_zone,abs(hz));
hold off;

% %askisi 2
 close all;

%the two orders that are asked
order = [2,16];
Ts = 0.2;
fs1 = 1/Ts;
omegac = 2;
fc= omegac/(2*pi);
ftel= fc/fs1;
omegatel=2*pi*ftel;
Rp = 3;


%samples
samples=256; 
w = [0:1/(samples): 1 - 1/samples];

%with the 1/samples we create the 256 samples that are asked

%creation of the chebyshev type 1 filter
[b,a]=cheby1(order(1),Rp,omegatel,'high');
cheb1 = freqz(b,a,samples);

[b1,a1]=cheby1(order(2),Rp,omegatel,'high');
cheb2 = freqz(b1,a1,samples);


freqz(b,a,samples)
figure

freqz(b1,a1,samples)
figure
hold on;
plot(w,mag2db(abs(cheb1)));
plot(w,mag2db(abs(cheb2)));
hold off;


    %askisi 3
    close all;

    samples = 500;
    Ts = 1/fs;
    t_zone = [0:samples-1]; %creating the zone for the xt
    f_zone2=[-fs/2 : fs/samples :fs/2-fs/samples];

    figure
    suptitle('Butterworth filter (attenuation 50)');
    subplot(4,1,1);
    xt=1+cos(1000*t_zone*Ts)+cos(16000*t_zone*Ts)+cos(30000*t_zone*Ts);
    plot(t_zone*Ts,xt);
    title('signal x(t)');

    subplot(4,1,2);
    xft = fftshift(fft(xt));
    plot(f_zone2,abs(xft));
    title('Fourier transform of signal x(t)');

    subplot(4,1,3);
    Ff=filter(numd,dend,abs(xft));
    plot(f_zone2,Ff);
    title('The filtered Fourier signal');

    subplot(4,1,4);
    Ff=filter(numd,dend,xt);
    plot(f_zone2,abs(Ff));
    title('The filtered signal x(t)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    Ts= 0.2;
    fs1 = 1/Ts;
    t_zone = [0:samples-1]; %creating the zone for the xt
    f_zone3=[-fs1/2 : fs1/samples :fs1/2-fs1/samples];


    figure
    suptitle('Chebyshev 16th order');
    
    subplot(4,1,1);
    xt2=1+cos(1.5*t_zone*Ts)+cos(5*t_zone*Ts);
    plot(t_zone*Ts,xt2);    
    title('signal x(t)');

    
    
    subplot(4,1,2);
    xft2 = fftshift(fft(xt2));
    plot(f_zone3,abs(xft2));
    title('Fourier transform of signal x(t)');

    subplot(4,1,3);
    Ff2=filter(b1,a1,abs(xft2));
    plot(f_zone3,Ff2);
    title('The filtered Fourier signal');

    subplot(4,1,4);
    Ff2=filter(b1,a1,xt);
    plot(f_zone3,abs(Ff2));
    title('The filtered signal x(t)');
