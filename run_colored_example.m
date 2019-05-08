clear
clc
close all;
%% read speech signals
[cleanspeech, samplingFreq] = audioread(['CleanSpeech.wav']);
%% generate noisy data
SNR=5;
load factoryNoise.mat
noise_s=resample(factoryNoise,samplingFreq,16000);
randstart=randi(length(noise_s)-length(cleanspeech)-100);
noise_seg=noise_s(randstart:randstart+length(cleanspeech)-1);
noise=addnoise_strict_snr(cleanspeech,noise_seg,SNR);
NoisySignal=cleanspeech+noise;
%% process the data 
% the third argument is the pre-whitening flag,
% when prew_flag=0,  pre-whitening will be disabled, and 
% when prew_flag=1,  pre-whitening will be enabled
plot_flag=1;
prew_flag=1;
F0_result=BF0NLS(NoisySignal,samplingFreq,plot_flag,prew_flag);