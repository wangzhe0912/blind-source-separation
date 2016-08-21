clc
m                         =2;%number of mixtures
n                         =2;%number of sources
 [S1,FS]                   =wavread('m1_1.5_m2_1.5_s1_002.wav'); %采样频率Fs;
 [S2,FS]                   =wavread('m1_1.5_m2_1.5_s2_002.wav'); %采样频率Fs;
    S1                    =S1( 110:84473);
    S2                    =S2( 129:84473+19);
 X                       =[S1'; S2'];
 X                       =X';
NFFT                      =1024;
OVERLAP                   =512;
N                         =500;
[Y1,Y2]                   = ica_f(X,NFFT,FS,OVERLAP,N);
% subplot(2,1,1);
% plot(Y1,'.');
% subplot(2,1,2);
% plot(Y2,'r');