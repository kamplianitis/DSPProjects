%Askisi 2
close all;
%clear all;

% create the several timezones
Ts= 1/1000;
Ts1 = 1/48;
Ts2= 1/24;
Ts3= 1/12;

% create the signal that is asked
tz=  0:Ts:0.5;
xt = 5*cos(24*pi*tz) - 2*sin(1.5*pi*tz);

%for the fft we declare the Nf( how many sample will take) 
% we transform the timezone
Nf= 1024;
fs= 1/Ts;
f_axis = -fs/2:fs/Nf: (fs/2)- (1/Nf);


%fast Fourier Transform
X_f = fftshift(fft(xt,Nf));
X_F = X_f*Ts;
F_axis = fs*f_axis;

plot(F_axis, abs(X_F))
title('Fourier Transform of the given signal');
xlabel('F');
ylabel('Xf[F]');

%timezones depending on the Ts
tz1= [0: Ts1: 0.5];
tz2= [0: Ts2: 0.5];
tz3= [0: Ts3: 0.5];

%our signals depending on the Ts

xt1 = 5*cos(24*pi*tz1) - 2*sin(1.5*pi*tz1);
xt2 = 5*cos(24*pi*tz2) - 2*sin(1.5*pi*tz2);
xt3 = 5*cos(24*pi*tz3) - 2*sin(1.5*pi*tz3);

% now we create the plots of the signals
figure
hold all;
plot(tz, xt, 'r' )
plot(tz1,xt1,'b')
legend('Ts = 1/1000', 'Ts = 1/48')
title('Simultaneous depict of the initial signal with the singal that we take by sampling with Ts = 1/48');
xlabel('t(second)');
ylabel('X(t)');

hold off;

figure
hold all;
plot(tz, xt, 'r' )
plot(tz2,xt2,'b')
legend('Ts = 1/1000', 'Ts = 1/24')
title('Simultaneous depict of the initial signal with the singal that we take by sampling with Ts = 1/24');
xlabel('t(second)');
ylabel('X(t)');

hold off;

figure
hold all;
plot(tz, xt, 'r' )
plot(tz3,xt3,'b')
legend('Ts = 1/1000', 'Ts = 1/12')
title('Simultaneous depict of the initial signal with the singal that we take by sampling with Ts = 1/12');
xlabel('t(second)');
ylabel('X(t)');
hold off;