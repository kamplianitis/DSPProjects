clear all; 
close all;

%sampling frequency ratio
fs = 1;
Ts = 1/fs;

%we use the tf function to simulate the tranfer function
num=[0, 0.2, 0];
dnum = [1,-0.7,-0.18];

H = tf(num, dnum, Ts)

%to create the zeroes- poles diagramm i have to find the roots of the num and dnum
rnum=roots(num);
rdnum=roots(dnum);
% poles and zeroes diagramm

zplane(num, dnum);

%the frequency zone that we asked for
fzone = [-pi:pi/128:pi];

figure
freqz(num, dnum)
figure
freqz(num, dnum, fzone)
%we use the tf function to simulate the tranfer function
%%%

num = [0,0,0.2,0];
dnum = [1,-1.7,0.52,0.18];

tf(num, dnum, Ts)
figure
zplane(num, dnum)
figure
freqz(num, dnum)
figure
freqz(num, dnum, fzone)

%% 2 %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define Symbolic Variable z
syms z

%Print transfer function in command line
transfFunc = (4 -3.5*z^-1)/(1 - 2.5*z^-1 + z^-2);
display('H(z)=');
pretty(transfFunc);

%Partial fraction decomposition
transfFuncNum = [4 -3.5 0];
transfFuncDenum = [1 -2.5 1];
[residues, poles] = residuez(transfFuncNum, transfFuncDenum);


%as soon as we know that we have 2 poles max we will have H1 + H2 that gives us the initial Hz
H1 = residues(1)/ (1-poles(1)*z^(-1));
H2 = residues(2)/ (1-poles(2)*z^(-1));

Hz = H1+H2;
display('=');
pretty(Hz)

figure
% Print poles and zeros
hold on
zplane(transfFuncNum,transfFuncDenum);
hold off
