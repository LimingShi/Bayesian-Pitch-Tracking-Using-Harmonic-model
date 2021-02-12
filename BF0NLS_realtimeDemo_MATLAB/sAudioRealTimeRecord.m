% Based on https://se.mathworks.com/help/audio/ref/audiodevicereader-system-object.html
clear
clc
close all

% setup
samplingFreq = 8000; % Hz
recordingTime = 10; % s
segmentTime = 0.025; % s
segmentLength = round(segmentTime*samplingFreq);
deviceName = 'Default'; % use getAudioDevices(audioDeviceReader) to see available devices
micChannelNo = 1;

% initialise the recording object
% too see avai
AudioRecObj = audioDeviceReader(...
    'SampleRate', samplingFreq, ...
    'SamplesPerFrame', segmentLength, ...
    'NumChannels', length(micChannelNo), ...
    'Device', deviceName);


% 
% AudioRecObj = audioDeviceReader(...
%     'SampleRate', samplingFreq, ...
%     'SamplesPerFrame', segmentLength, ...
%     'NumChannels', length(micChannelNo), ...
%     'Device', deviceName);
% Call setup to reduce the computational load of initialization in an 
% audio stream loop.
AudioRecObj.setup();

% make a recording
nSegments = floor(recordingTime/segmentTime);
recordedData = nan(segmentLength, nSegments);
actFrameTimes = nan(1,nSegments);
for ii = 1:nSegments
%     tic
    [recordedSegment, numOverrun] = AudioRecObj();
%     toc
    % here we do some processing
    recordedData(:,ii) = recordedSegment(:);
    pause(0.01)
    % end of processing
    actFrameTimes(ii) = toc;
end
AudioRecObj.release();

%% plotting
figure(1)
plot((1:recordingTime*samplingFreq)'/samplingFreq, recordedData(:))
xlabel('time [s]')
ylabel('value [.]')
title('Recorded signal')

figure(2)
plot((1:nSegments)*segmentTime, cumsum(actFrameTimes),'.')
hold on
plot((1:nSegments)*segmentTime, (1:nSegments)*segmentTime, '--')
hold off
legend('Processed segment time', 'Real-time limit')
xlabel('sample recording time [s]')
ylabel('finished processing a segment [s]')
title('Processing time vs recorded time')