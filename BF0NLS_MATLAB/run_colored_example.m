clear
clc
close all;
%% read speech signals
[cleanspeech, samplingFreq] = audioread(['CleanSpeech.wav']);
%% generate noisy data
SNR=0;
load factoryNoise.mat
noise=resample(factoryNoise,samplingFreq,16000);
randstart=randi(length(noise)-length(cleanspeech)-100);
noise_seg=noise(randstart:randstart+length(cleanspeech)-1);
scaled_noise=addnoise(cleanspeech,noise_seg,SNR);
NoisySignal=cleanspeech+scaled_noise;
%% process the data 
% the third input argument is the plot flag,
% when plot_flag=0,  results will not be given in figure;
% when plot_flag=1,  results will be given in figure;

% the fourth input argument is the pre-whitening flag,
% when prew_flag=0,  pre-whitening will be disabled, and 
% when prew_flag=1,  pre-whitening will be enabled

% the fifth input argument is the parameter set you can customize yourself,
plot_flag=1;
prew_flag=1;
tic
F0_result=BF0NLS(NoisySignal,samplingFreq,plot_flag,prew_flag);
toc