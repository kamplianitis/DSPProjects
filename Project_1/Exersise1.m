% the signals that we choose is the u[n] and the δ[n]
close all;
%clear all;

%sampling zone
N=100;
n = -5: N;

% now i create the signals

% δ[n]
delta=(n==0);
stem(n,delta)
title('Dirac pulse');
xlabel('n');
ylabel('δ[n]');


% now we create the u[n]
u1 = (n>=0); 
figure
stem(n, u1)
title('Unit step');
xlabel('n');
ylabel('u[n]');



%convolution without the conv function
%first things first we create the sampling axis of the convolution
nconv = 2*min(n): 2*max(n);
nconvl = length(nconv);

nl=length(n);
%zeropadding
X=[u1,zeros(1,nl)];
H=[delta,zeros(1,nl)];
%convolution
for i= 1:nl+nl-1
    Y(i)=0;
    for j=1:nl
        if(i-j+1>0)
            Y(i)=Y(i)+X(j)*H(i-j+1);
        end
    end
end

figure
stem(nconv, Y)
title('Convolution without the conv function');
xlabel('n');
ylabel('δ[n]*u[n]');

%convolution using the conv function
assurance=conv(double(u1),double(delta));
figure
stem(min(n)+min(n):max(n)+max(n),assurance)
title('Convolution with the conv function');
xlabel('n');
ylabel('δ[n]*u[n]');

close all;

Nf1=1024; % declares the number of the samples 
f_axis = -1/2:1/Nf1:1/2 - (1/Nf1); % axis in the frequencies

% fast fourier transform
X_fu1 = fftshift(fft(u1,Nf1));
X_fdelta = fftshift(fft(delta,Nf1));
figure 
plot(f_axis, abs(X_fu1))
title('Fourier Transform of u[n]');
xlabel('F');
ylabel('Xf1[F]');

figure 
plot(f_axis, abs(X_fdelta))
title('Fourier Transform of δ[n]');
xlabel('F');
ylabel('Xf2[F]');

Xf_final = X_fu1.*X_fdelta;
figure
plot(f_axis, abs(Xf_final))
title('Fourier Transform of the multiplication of the signals');
xlabel('F');
ylabel('Xf1[F]xXf2[F]');


% how we make sure the result we have is right
Xf= fftshift(fft(assurance,Nf1));
figure
plot(f_axis, abs(Xf));
title('Fourier Transform of the singal that we have from the convolution in the time var');
xlabel('F');
ylabel('Xf3[F]');


