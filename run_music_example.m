clear
clc
close all;
%% read speech signals
[vib_flute, samplingFreq] = audioread(['vib_flute.wav']);vib_flute=vib_flute(:,1);
%% generate noisy data
SNR=0;
noise=randn(size(vib_flute)); % White Gaussian noise
scaled_noise=addnoise(vib_flute,noise,SNR);
NoisySignal=vib_flute+scaled_noise;
%% process the data 
% the third argument is the pre-whitening flag,
% when prew_flag=0,  pre-whitening will be disabled, and 
% when prew_flag=1,  pre-whitening will be enabled
plot_flag=1;
prew_flag=0;
Par_set.f0Bounds = [100, 1000]; % Pitch boundaries for music signal (prior information)
tic
F0_result=BF0NLS(NoisySignal,samplingFreq,plot_flag,prew_flag,Par_set);
toc



