function PitchTraking
close all;
warning off;
%Syntax: PitchTraking
% Run PitchTraking by typing "PitchTraking" in your command line
%
% Liming Shi 3
% Generated: 20.02.2019
%



% Initialize and hide the GUI as it is being constructed.
aud = audiodevinfo;
if isempty(ver('Signal'))
	errordlg('Signal Processing Toolbox required','!! Error !!')
elseif isempty(aud.input)
    errordlg('No input device found for audio recording. After plugged in you have to restart Matlab.' ,'!! Error !!')
else

f=figure('Visible','off','Position',[0,0,900,600],'Units','normalized','Toolbar','figure','numbertitle','on','Visible','off');%,'PaperSize',[20,13]
set(f,'ToolBar', 'none');
set(f,'MenuBar', 'none');
% f = figure('numbertitle', 'off');

% t
hp1 = uipanel('Title','Length of record','FontSize',14,'BackgroundColor','white','Units','normalized','Position',[0.1 0.5 .6 .2],'FontUnits','normalized','visible','off');
hp2 = uipanel('Title','Display length','FontSize',14,'BackgroundColor','white','Units','normalized','Position',[0.1 0.25 .6 .2],'FontUnits','normalized','visible','off');
hp3 = uipanel('Title','Frequency window','FontSize',14,'BackgroundColor','white','Units','normalized','Position',[0.1 0.1 .25 .1],'FontUnits','normalized','visible','off');
TextH1 = text(0.2, 0.4, 'Authors: Liming Shi, Jesper Kjaer Nielsen,','FontSize',20,'Units','normalized');
TextH2 = text(0.1, 0.32, 'Jesper Rindom Jensen and Mads Graesboll Christensen','FontSize',20,'Units','normalized');
TextH3 = text(0.17, 0.24, 'Audio Analysis Lab, CREATE, Aalborg University','FontSize',20,'Units','normalized');
TextH4 = text(0.29, 0.16, '\{ls, jkn, jrj, mgc\}@create.aau.dk','FontSize',20,'Units','normalized');
TextH5 = text(0.32, 0.08, 'https://audio.create.aau.dk','FontSize',20,'Units','normalized');

TextH6 = text(0.75, 0.75, 'White noise','FontSize',13,'Units','normalized','Color','r');
TextH7 = text(0.75, 0.70, 'Infinite SNR','FontSize',13,'Units','normalized','Color','r');
TextH8 = text(0.75, 0.81, 'No','FontSize',13,'Units','normalized','Color','r');
TextH9 = text(0.25, 0.81, 'No','FontSize',13,'Units','normalized','Color','r');
axis off;
% set(TextH,'ToolBar', 'none');
% set(TextH,'MenuBar', 'none');
% 
% 
% % Construct the components.
% %Text
% htext=uicontrol('Parent',hp2,'Style','text','String','-','FontSize',10,'Units','normalized','Position',[.4,.1,.1,.6],'FontUnits','normalized','HorizontalAlignment','center');
% htextHz=uicontrol('Parent',hp2,'Style','text','String','Hz','FontSize',10,'Units','normalized','Position',[.9,.1,.1,.6],'FontUnits','normalized','HorizontalAlignment','center');

%Button
pre_button_real = uicontrol('String','Prewhitening','BackgroundColor','green','FontSize',14,'Units','normalized','Position',[.2,.75,.1,.05],'FontUnits','normalized','Callback', @buttionwhite2);
hbuttonStart = uicontrol('String','Start live recording','BackgroundColor','green','FontSize',14,'Units','normalized','Position',[.2,.55,.2,.2],'FontUnits','normalized','Callback', @buttonStart_Callback);
hbuttonStartUnlimited = uicontrol('String','Insert speech recording','BackgroundColor','green','FontSize',14,'Units','normalized','Position',[.6,.55,.2,.1],'FontUnits','normalized','visible','on','Callback', @buttonoffline);
pre_button=uicontrol('String','Prewhitening','BackgroundColor','green','FontSize',14,'Units','normalized','Position',[.6,.75,.10,.05],'FontUnits','normalized','visible','on','Callback', @buttionwhite);
noiseSNRCreation=uicontrol('String','SNR','BackgroundColor','green','FontSize',14,'Units','normalized','Position',[.6,.65,.10,.05],'FontUnits','normalized','visible','on','Callback', @buttionNOISENSR);
noiseCreation=uicontrol('String','Noise Type','BackgroundColor','green','FontSize',14,'Units','normalized','Position',[.6,.7,.10,.05],'FontUnits','normalized','visible','on','Callback', @buttionNOISE);
hbuttonStartAgain = uicontrol('String','back','FontSize',10,'Units','normalized','Position',[.85,.01,.125,.05],'FontUnits','normalized','visible','off','Callback', @buttonStartAgain_Callback);
hbuttonPlay1 = uicontrol('String','play again','FontSize',10,'Units','normalized','Position',[.70,.01,.125,.05],'FontUnits','normalized','visible','off','Callback', @buttonPlay_Callback);
% hbuttonStartAgain = uicontrol('String','back','FontSize',10,'Units','normalized','Position',[.85,.01,.125,.05],'FontUnits','normalized','visible','off','Callback', @buttonStartAgain_Callback);
hbuttonStop = uicontrol('String','stop recording','FontSize',10,'Units','normalized','Position',[.85,.95,.125,.05],'FontUnits','normalized','visible','off','Callback', @buttonStop_Callback);
hbuttonSave = uicontrol('String','save','FontSize',10,'Units','normalized','Position',[.65,.01,.1,.05],'FontUnits','normalized','visible','off','Callback', @buttonSave_Callback);
hbuttonSaveAs = uicontrol('String','save as','FontSize',10,'Units','normalized','Position',[.75,.01,.1,.05],'FontUnits','normalized','visible','off','Callback', @buttonSaveAs_Callback);
hbuttonPlay = uicontrol('String','play','FontSize',10,'Units','normalized','Position',[.05,.01,.1,.05],'FontUnits','normalized','visible','off','Callback', @buttonPlay_Callback);
hbuttonPlayAll = uicontrol('String','play all','FontSize',10,'Units','normalized','Position',[.15,.01,.1,.05],'FontUnits','normalized','visible','off','Callback', @buttonPlayAll_Callback);
hbuttonShowFigures = uicontrol('String','open figures','FontSize',10,'Units','normalized','Position',[.3,.01,.125,.05],'FontUnits','normalized','visible','off','Callback', @buttonShowFigures_Callback);


%Slider
hsliderRecordLength = uicontrol('Parent',hp1,'Style','slider','Min',1,'Max',30,'SliderStep',[1 1]./29,'Value',10,'Units','normalized','Position',[.1 .5 .8 .3],'FontUnits','normalized','Callback',@sliderRecordLength_Callback);  
hsliderShowLength = uicontrol('Parent',hp2,'Style','slider','Min',1,'Max',30,'SliderStep',[1 1]./29,'Value',5,'Units','normalized','Position',[.1 .5 .8 .3],'FontUnits','normalized','Callback',@sliderShowLength_Callback);  

%Edit fields
heditRecordLength=uicontrol('Parent',hp1,'Style','edit','String',22,'Units','normalized','Position',[.5 .2 .4 .3],'FontUnits','normalized','Callback',@editRecordLength_Callback);  
heditShowLength=uicontrol('Parent',hp2,'Style','edit','String',5,'Units','normalized','Position',[.5 .2 .4 .3],'FontUnits','normalized','Callback',@editShowLength_Callback);  
heditFrequencyWindow1=uicontrol('Parent',hp3,'Style','edit','String',0,'Units','normalized','Position',[.1 .1 .3 .8],'FontUnits','normalized','Callback',@editFrequencyWindow1_Callback);  
heditFrequencyWindow2=uicontrol('Parent',hp3,'Style','edit','String',3000,'Units','normalized','Position',[.5 .1 .4 .8],'FontUnits','normalized','Callback',@editFrequencyWindow2_Callback);  


ha1=axes('Units','Pixels','Position',[70,220,800,150],'Units','normalized','FontUnits','normalized','Layer','top','visible','off'); 
ha2=axes('Units','Pixels','Position',[70,50,800,150],'Units','normalized','FontUnits','normalized','Layer','top','visible','off'); 
ha3=axes('Units','Pixels','Position',[70,390,800,75],'Units','normalized','FontUnits','normalized','Layer','top','visible','off'); 
ha5=axes('Units','Pixels','Position',[70,485,800,75],'Units','normalized','FontUnits','normalized','Layer','top','visible','off'); 
ha4=axes('Units','Pixels','Position',[70,220,800,150],'Units','normalized','FontUnits','normalized','Layer','top','visible','off'); 

% Global variables
global RecordLength ShowLength FrequencyWindow1 FrequencyWindow2 myRecording fs mag Release Unlimited stoprec timeVector_final maxNoHarmonics used_noise used_snr noise used_prew used_prew2;
% global FullDictionarySpec1 numberSpeechVectors NumberNoiseCodebookEntries;

% speech_codebook=load('ls_codebook.mat');     %%32 speech codebook entries (from 4 
%%speakers of CMU arctic database, 54 minutes of speech material) 
% noise_codebook=load('ls_babble_codebook.mat');  %%14 noise codebook entries (trained on


% FullDictionarySpec1 = [speech_codebook.GeneralSpeechDictionary noise_codebook.Aurora8kHzNoiseDictionary];
% numberSpeechVectors = size(speech_codebook.GeneralSpeechDictionary,2); 
% NumberNoiseCodebookEntries = size(noise_codebook.Aurora8kHzNoiseDictionary,2); 
%     plotHandles = cell(1,1);

used_noise='White';
used_snr=inf;
used_prew=0;
used_prew2=0;

timerPeriod = .032; % s
samplingFreq=8000;
segmentTime = 0.032; %     seconds
segmentLength = round(segmentTime*samplingFreq); % samples
f0Bounds = [70, 400]/samplingFreq; % cycles/sample
maxNoHarmonics = 10;
segmentShift=timerPeriod;
f0Estimator = BayesianfastF0NLS(segmentLength, maxNoHarmonics, f0Bounds,5/samplingFreq,0.7);
%Initialisation
RecordLength = str2double(get(heditRecordLength,'String'));
ShowLength = str2double(get(heditShowLength,'String'));
FrequencyWindow1 = str2double(get(heditFrequencyWindow1,'String'));
FrequencyWindow2 = str2double(get(heditFrequencyWindow2,'String'));
Release = version('-release');
Release = str2double(Release(1:4));
Unlimited = 0;
stoprec = 0;

nShift = round(segmentShift*samplingFreq); ; % samples
nSegments = floor((RecordLength*samplingFreq-segmentLength)/nShift)+1;
timeVector_final = (0:nSegments-1)*segmentShift+segmentTime/2;

% Initialize the GUI.
set(f,'Name','Pitch tracking')  % Assign the GUI a name to appear in the window title.
movegui(f,'center')     % Move the GUI to the center of the screen.
set(f,'Visible','on');	% Make the GUI visible.

% axes(f)
% matlabImage = imread('logo.png');
% logo_handle=imshow(matlabImage);
% axis off
% axis image

end

%% BUTTONS
% Push button callbacks. Each callback plots current_data in the
% specified plot type.
function buttonStart_Callback(~,~)
    
    

    RecordLength = str2double(get(heditRecordLength,'String'));
    ShowLength = str2double(get(heditShowLength,'String'));
    FrequencyWindow1 = str2double(get(heditFrequencyWindow1,'String'));
    FrequencyWindow2 = str2double(get(heditFrequencyWindow2,'String'));
    
    set(hp1,'visible','off');    
    set(hp2,'visible','off');
    set(hp3,'visible','off');
    
    set(ha1,'visible','on');    
    set(ha2,'visible','on');  
    set(ha3,'visible','on');
    set(ha4,'visible','on');
    set(ha5,'visible','on');


    set(hbuttonStart,'visible','off') 
    set(hbuttonStartUnlimited,'visible','off')
    set(noiseCreation,'visible','off')     
    set(noiseSNRCreation,'visible','off') 
    set(hbuttonPlay,'visible','off') 
    set(hbuttonPlayAll,'visible','off') 
    set(hbuttonShowFigures,'visible','off') 
    set(hbuttonSave,'visible','off') 
    set(hbuttonSaveAs,'visible','off')   
    set(pre_button,'visible','off')   
    set(pre_button_real,'visible','off') 
    
    set(TextH1,'Visible','off')
    set(TextH2,'Visible','off')
    set(TextH3,'Visible','off')
    set(TextH4,'Visible','off')
    set(TextH5,'Visible','off')
    
    set(TextH6,'Visible','off')
    set(TextH7,'Visible','off')
       set(TextH8,'Visible','off')
    set(TextH9,'Visible','off')
    

    set(hbuttonStop,'visible','off')
    timerPeriod = 0.045; % s
    samplingFreq=8000;
    segmentTime = 0.0450; %     seconds
    segmentLength = round(segmentTime*samplingFreq); % samples
    f0Bounds = [70, 500]/samplingFreq; % cycles/sample
    maxNoHarmonics = 10;
    segmentShift=timerPeriod;
%     f0Estimator.unvoicing_scaled_alpha_buffer=[];
%     f0Estimator.scaled_alpha_buffer=nan;
%     f0Estimator.scaled_alpha_buffer2=nan;
    f0Estimator = BayesianfastF0NLS(segmentLength, maxNoHarmonics, f0Bounds,5/samplingFreq,0.7);
    %Initialisation
    RecordLength = str2double(get(heditRecordLength,'String'));
    ShowLength = str2double(get(heditShowLength,'String'));
    FrequencyWindow1 = str2double(get(heditFrequencyWindow1,'String'));
    FrequencyWindow2 = str2double(get(heditFrequencyWindow2,'String'));
    Release = version('-release');
    Release = str2double(Release(1:4));
    Unlimited = 0;
    stoprec = 0;
    
    nShift = round(segmentShift*samplingFreq); ; % samples
    nSegments = floor((RecordLength*samplingFreq-segmentLength)/nShift)+1;
    timeVector_final = (0:nSegments-1)*segmentShift+segmentTime/2;

    
    
    
    
    
    
    

    liverecording
%     myRecording=
    set(hbuttonStartAgain,'visible','on')
    set(hbuttonPlay1,'visible','on')
end

function buttonoffline(~,~)
    [filename,filepath]=uigetfile({'*.wav;*.WAV'},'Select an audio recording');
    if and(sum(filename==0),sum(filepath==0));
        return;
    end
    filename=strcat(filepath,filename);
    [speechSignal,samplingfrequency]=audioread(filename);
    if samplingfrequency>=8000
        speechSignal=resample(speechSignal,8000,samplingfrequency);
        samplingfrequency=8000;
    end
    
    fs=samplingfrequency;
    noise=get_noise_signals(used_noise);
    noise.noisetype=fields(noise);
    noisesig=noise.(noise.noisetype{1});
    noisesig=noisesig/max(noisesig);noisesig=noisesig(:);
    noisesig=resample(noisesig,samplingfrequency,19980);
    rand_num=randi(length(noisesig)-length(speechSignal)-100);
    noisesig=noisesig(rand_num:length(speechSignal)+rand_num-1);
    
    noisesig=addnoise_strict_snr(speechSignal,noisesig,used_snr);
    speechSignalClean=speechSignal;
    speechSignal=speechSignal+noisesig;
%     speechSignal=resample(speechSignal,fs,samplingfrequency);
    % normalization step
    sumE = sqrt(speechSignal'*speechSignal/length(speechSignal));
    scale = sqrt(3.1623e-5)/sumE; % scale to -45-dB loudness level
    speechSignal=speechSignal*scale;
    myRecording=speechSignal;
    
    nData = length(speechSignal);
    
    % set up
    segmentTime = 0.0375; %     seconds
    segmentLength = round(segmentTime*samplingfrequency); % samples
    segmentShift = 0.0375; % seconds
    nShift = round(segmentShift*samplingfrequency); % samples
    nSegments = floor((nData+segmentLength/2-segmentLength)/nShift)+1;
    f0Bounds = [70, 400]/samplingfrequency; % cycles/sample
    if used_prew
    maxNoHarmonics = 20;
    else
     maxNoHarmonics = 10;
    end
    f0Estimator = BayesianfastF0NLS(segmentLength, maxNoHarmonics, f0Bounds,5/samplingfrequency,.7);
    speechSignal=[zeros(round(segmentLength/2),1);speechSignal];
    speechSignalClean=[zeros(round(segmentLength/2),1);speechSignalClean];
    timeVector_final = (0:nSegments-1)*segmentShift+segmentTime/2;
    
    

    RecordLength = str2double(get(heditRecordLength,'String'));
    ShowLength = str2double(get(heditShowLength,'String'));
    FrequencyWindow1 = str2double(get(heditFrequencyWindow1,'String'));
    FrequencyWindow2 = str2double(get(heditFrequencyWindow2,'String'));
    
    set(hp1,'visible','off');    
    set(hp2,'visible','off');
    set(hp3,'visible','off');
    
    set(ha1,'visible','on');    
    set(ha2,'visible','on');  
    set(ha3,'visible','on');
    set(ha4,'visible','on');
    set(ha5,'visible','on');

    set(hbuttonStart,'visible','off') 
    set(hbuttonStartUnlimited,'visible','off')
    set(hbuttonPlay,'visible','off') 
    set(noiseCreation,'visible','off') 
    set(noiseSNRCreation,'visible','off') 
    set(hbuttonPlayAll,'visible','off') 
    set(hbuttonShowFigures,'visible','off') 
    set(hbuttonSave,'visible','off') 
    set(hbuttonSaveAs,'visible','off')   
    set(hbuttonStop,'visible','off')
    set(pre_button,'visible','off')
    set(pre_button_real,'visible','off')    
    
    set(TextH1,'Visible','off')
    set(TextH2,'Visible','off')
    set(TextH3,'Visible','off')
    set(TextH4,'Visible','off')
    set(TextH5,'Visible','off')
    set(TextH6,'Visible','off')
    set(TextH7,'Visible','off')
    set(TextH8,'Visible','off')
    set(TextH9,'Visible','off')
    
    

    offlinerecording(speechSignal,speechSignalClean,samplingfrequency,nSegments);

    set(hbuttonStartAgain,'visible','on')
    set(hbuttonPlay1,'visible','on')
end


function buttonStartAgain_Callback(~,~)

    set(ha1,'visible','off');    
    set(ha2,'visible','off');
    set(ha3,'visible','off');
    set(ha4,'visible','off');
    set(ha5,'visible','off');
    cla(ha1)
    cla(ha2)
    cla(ha3)
    cla(ha4)
    cla(ha5)
    
    set(hp1,'visible','off');    
    set(hp2,'visible','off');
    set(hp3,'visible','off');
    set(hbuttonPlay1,'visible','off')
    

    set(hbuttonStart,'visible','on') 
    set(noiseCreation,'visible','on') 
    set(noiseSNRCreation,'visible','on') 
    set(hbuttonStartUnlimited,'visible','on') 
    set(pre_button,'visible','on') 
    set(pre_button_real,'visible','on')     
    
    set(hbuttonPlay,'visible','off') 
    set(hbuttonPlayAll,'visible','off') 
    set(hbuttonShowFigures,'visible','off') 
    set(hbuttonSave,'visible','off') 
    set(hbuttonSaveAs,'visible','off')    
    set(hbuttonStartAgain,'visible','off')
    set(hbuttonStop,'visible','off')
    
        set(TextH1,'Visible','on')
    set(TextH2,'Visible','on')
    set(TextH3,'Visible','on')
    set(TextH4,'Visible','on')
    set(TextH5,'Visible','on')
    
        
    set(TextH6,'Visible','on')
    set(TextH7,'Visible','on')
    set(TextH8,'Visible','on')
    set(TextH9,'Visible','on')    
    

end

function buttonStop_Callback(~,~)
    stoprec = 1;
end

function buttonSave_Callback(~,~)   
    [y,m,d,h,min,sec]=datevec(now);
    audiowrite([num2str(y,'%04.0f') num2str(m,'%02.0f') num2str(d,'%02.0f') '-' num2str(h,'%02.0f') ''''  num2str(min,'%02.0f') '''''' num2str(floor(sec),'%02.0f') '.wav'],myRecording,fs);
end

function buttonSaveAs_Callback(~,~) 
    [y,m,d,h,min,sec]=datevec(now);    
    [file,path] = uiputfile([num2str(y,'%04.0f') num2str(m,'%02.0f') num2str(d,'%02.0f') '-' num2str(h,'%02.0f') ''''  num2str(min,'%02.0f') '''''' num2str(floor(sec),'%02.0f') '.wav'],'Save record');
    audiowrite([path file],myRecording,fs);
end

function buttonPlay_Callback(~,~)   
    xx = round(get(ha1,'Xlim')*fs);
    soundsc(myRecording(max(1,xx(1)):min(xx(2),size(myRecording,1))), fs);   
    audiowrite('myRecording.wav',myRecording,fs);
end

function buttonPlayAll_Callback(~,~)   
    sound(myRecording, fs);
end

function buttonShowFigures_Callback(~,~)   
    figure    
    plot((1:size(myRecording,1))./fs,myRecording)
    ylim([-1.2 1.2]*mag)
    xlim([0 max(size(myRecording,1)/fs,ShowLength)])
    box on
    grid on
    
    figure
    spectrogram(myRecording,2^9,2^7,2^12,fs)
    view(-90,90) 
    set(gca, 'ydir', 'reverse', 'box', 'on');        
    % spectrogram since R2015a with legend and in kHz       
    if Release<2015
        set(gca, 'xlim', [FrequencyWindow1 FrequencyWindow2], 'ylim' , [0 max(size(myRecording,1)/fs,ShowLength)]);
    else
        set(gca, 'xlim', [FrequencyWindow1 FrequencyWindow2]/1000, 'ylim' , [0 max(size(myRecording,1)/fs,ShowLength)]);
    end  
end

%% EDITFIELDS
function editRecordLength_Callback(~,~)
      RecordLength = str2double(get(heditRecordLength,'String'));
      set(hsliderRecordLength,'value',RecordLength);
end

function editShowLength_Callback(~,~)
      ShowLength = str2double(get(heditShowLength,'String'));
      set(hsliderShowLength,'value',ShowLength);
end

function editFrequencyWindow1_Callback(~,~)
      FrequencyWindow1 = str2double(get(heditFrequencyWindow1,'String'));
end
function editFrequencyWindow2_Callback(~,~)
      FrequencyWindow2 = str2double(get(heditFrequencyWindow2,'String'));
end

%% SLIDER
function sliderRecordLength_Callback(~,~)
      RecordLength = round(get(hsliderRecordLength,'value'));
      set(heditRecordLength,'string',num2str(RecordLength));
end

function sliderShowLength_Callback(~,~)
      ShowLength = round(get(hsliderShowLength,'value'));
      set(heditShowLength,'string',num2str(ShowLength));
end


%% GENERAL FUNCTIONS
function liverecording
%     RecordLength=RecordLength;
    
    fs = 8000;
    nBits = 16;
    mag = 1.05;
%     time_vec= [0:(RecordLength*fs-1)]/fs;
    plot(ha2,0,0);
    ylim(ha2,[-mag mag])
    xlim(ha2,[0 RecordLength-2])
    xlabel(ha2,'Time [s]')
    
    
    idx_last = 1;
    plot_freq=500;
    freq_ind=[0:512]/512*fs/2;
    [~,inx]= min(abs(freq_ind-plot_freq));
    imagesc(ha1,timeVector_final, freq_ind(1:inx),nan(inx,length(timeVector_final)));
    ylim(ha1,[0,plot_freq])
    xlim(ha1,[0 RecordLength-2])
     axis(ha1,'xy');
     hold on;
     plot(ha4,timeVector_final,nan(length(timeVector_final),1),'r-', 'linewidth',2);
      ylim(ha4,[0,plot_freq])
    xlim(ha4,[0 RecordLength-2])
     axis(ha4,'off');
    
    
%     xlabel(ha3,'Time [s]')
%     idx_last = 1;
    
    plot(ha3,0,0);
    ylim(ha3,[0 1])
    xlim(ha3,[0 RecordLength-2])
    
     plot(ha5,0,0);
    ylim(ha5,[0 maxNoHarmonics])
    xlim(ha5,[0 RecordLength-2])

    
    
    
    AudioRecObj = audioDeviceReader(...
    'SampleRate', fs, ...
    'SamplesPerFrame', segmentLength, ...
    'NumChannels', 1, ...
    'Device', 'Default');

  
% 
% if used_prew2==1
% %     recObj.method=@processData_MMSE;
%     
% %    recObj.set('TimerPeriod', timerPeriod, ...
% %     'TimerFcn', {@processData_MMSE});
% else
% %     recObj.method=@processData_orig;
% %     recObj.set('TimerPeriod', timerPeriod, ...
% %     'TimerFcn', {@processData_orig});
% end

    segmentLength=segmentTime*samplingFreq;
    recObj.UserData{1}=1;
    recObj.UserData{2}=1:segmentLength;
    recObj.UserData{3}=segmentTime*samplingFreq;
    recObj.UserData{4}=@f0Estimator.estimate;
    recObj.UserData{5}=samplingFreq;
    recObj.UserData{6}=[];
    recObj.UserData{8}=ha1;
    recObj.UserData{9}=ha2;
    recObj.UserData{10}=ha3;
    recObj.UserData{15}=ha4;
    recObj.UserData{17}=ha5;
    recObj.UserData{11}=RecordLength;
    recObj.UserData{20}=timeVector_final;
    recObj.UserData{21}=freq_ind(1:inx);
    recObj.UserData{22}=inx;
    recObj.UserData{23}=plot_freq;
    
    recObj.UserData{24}=maxNoHarmonics;
    recObj.UserData{25}=[];
    recObj.UserData{26}=[];
    
    AudioRecObj.setup();
    myRecording1=nan(segmentLength,nSegments);
    
    for ii = 1:nSegments
        [recordedSegment] = AudioRecObj();
        recObj.frame=recordedSegment;
        % here we do some processing
        if used_prew2==1
        recObj=processData_MMSE(recObj);
        else
            recObj=processData_orig(recObj);
        end
        pause(eps);
        % end of processing
        myRecording1(:,ii)=recordedSegment;
%         actFrameTimes(ii) = toc;
    end
    myRecording=myRecording1(:);
    AudioRecObj.release();
    
    
    
    
    
    
%     recordblocking(recObj, RecordLength);
%     
%     myRecording=getaudiodata(recObj);
%     if stoprec==1
%         stop(recObj);
%     end

    endrecording
end



function offlinerecording(resampled,clean,sampfreq,nSegments)
%     profile on
    RecordLength_offline=length(resampled)/sampfreq;
    

    nBits = 16;
    mag = 1.05;
%     time_vec= [0:(RecordLength*fs-1)]/fs;
    plot(ha2,0,0);
    ylim(ha2,[-mag mag])
    xlim(ha2,[0 RecordLength_offline])
    xlabel(ha2,'Time [s]')
    
    
    idx_last = 1;
    plot_freq=500;
    freq_ind=[0:128]/128*sampfreq/2;
    [~,inx]= min(abs(freq_ind-plot_freq));
    imagesc(ha1,timeVector_final, freq_ind(1:inx),nan(inx,length(timeVector_final)));
    ylim(ha1,[0,plot_freq])
    xlim(ha1,[0 RecordLength_offline])
     axis(ha1,'xy');
     hold on;
     plot(ha4,timeVector_final,nan(length(timeVector_final),1),'r-', 'linewidth',2);
      ylim(ha4,[0,plot_freq])
    xlim(ha4,[0 RecordLength_offline])
     axis(ha4,'off');
    
    
%     xlabel(ha3,'Time [s]')
%     idx_last = 1;
    
    plot(ha3,0,0);
    ylim(ha3,[0 1])
    xlim(ha3,[0 RecordLength_offline])
    
     plot(ha5,0,0);
    ylim(ha5,[0 maxNoHarmonics])
    xlim(ha5,[0 RecordLength_offline])
%     xlabel(ha3,'Time [s]')
%     idx_last = 1;
    timerPeriod=segmentTime;
%     recObj = audiorecorder(fs,nBits,1,0);

%     recObj.set('TimerPeriod', timerPeriod, ...
%     'TimerFcn', {@processData}, ...
%     'StartFcn', 'disp("Start recording ...")', ...
%     'StopFcn', 'disp("... finished recording. ")');


%     recObj.set('TimerPeriod', timerPeriod, ...
%     'TimerFcn', {@offlineprocessData});

    segmentLength=segmentTime*sampfreq;
    recObj.UserData{1}=1;
    recObj.UserData{2}=1:segmentLength;
    recObj.UserData{3}=segmentTime*sampfreq;
    recObj.UserData{4}=@f0Estimator.estimate;
    recObj.UserData{5}=sampfreq;
    recObj.UserData{6}=[];
    recObj.UserData{8}=ha1;
    recObj.UserData{9}=ha2;
    recObj.UserData{10}=ha3;
    recObj.UserData{15}=ha4;
    recObj.UserData{17}=ha5;
    recObj.UserData{11}=RecordLength_offline;
    recObj.UserData{20}=timeVector_final;
    recObj.UserData{21}=freq_ind(1:inx);
    recObj.UserData{22}=inx;
    recObj.UserData{23}=plot_freq;
    
    recObj.UserData{24}=maxNoHarmonics;
    recObj.UserData{25}=resampled;
    recObj.UserData{26}=clean;
%      recObj.UserData{12}=plotHandles;
    

%     pause(2)
    soundsc(resampled,fs);
if used_prew==1
    for ii=1:nSegments
%         recObj=offlineprocessData(recObj,FullDictionarySpec1,numberSpeechVectors,NumberNoiseCodebookEntries);
        tic
        recObj=offlineprocessData_MMSE(recObj);
        t1=toc;
        if t1<segmentTime
%             pause(3*(segmentTime-t1)/4);
        end
    end
else
    for ii=1:nSegments
        tic
        recObj=offlineprocessData_orig(recObj);
        t1=toc;
        if t1<segmentTime
            pause(eps);
        end
        
    end
end
    
%     recordblocking(recObj, RecordLength_offline);
    
% 
%     if stoprec==1
%         stop(recObj);
%     end

    endrecording
    
%     profile viewer
end



function endrecording
    set(hbuttonPlay,'visible','off')
    set(hbuttonPlayAll,'visible','off') 
    set(hbuttonShowFigures,'visible','off') 
    set(hbuttonSave,'visible','off')
    set(hbuttonSaveAs,'visible','off')
    stoprec = 0;
    Unlimited = 0;
end



function buttionNOISE(~,~)
 list = {'White','Babble','Factory',...                   
'Pink','F16','Leopard','M109','Gust wind','Howl wind'};
[indx,tf] = listdlg('PromptString','Select a noise type:','ListString',list,'SelectionMode','single','ListSize',[100,95]);
%  list2 = {'-5 dB','0 dB','10 dB',...                   
% '20 dB','30 dB','Infinite'};
% list_temp={-5,0,10,20,30,inf};
% [indx2,tf2] = listdlg('PromptString','Select the SNR:','ListString',list2,'SelectionMode','single','ListSize',[100,95]);
if tf~=0;
set(TextH6,'String',strcat(list{indx},' noise'));
used_noise=list{indx};
end
% if tf2~=0
% set(TextH7,'String',strcat(list2{indx2},' SNR'));
% used_snr=list_temp{indx2};
% end


end

function buttionNOISENSR(~,~)
%  list = {'White','Babble','Factory',...                   
% 'Pink','F16','Leopard','M109'};
% [indx,tf] = listdlg('PromptString','Select a noise type:','ListString',list,'SelectionMode','single','ListSize',[100,95]);
 list2 = {'-5 dB','0 dB','5 dB',...                   
'10 dB','20 dB','Infinite'};
list_temp={-5,0,5,10,20,inf};
[indx2,tf2] = listdlg('PromptString','Select the SNR:','ListString',list2,'SelectionMode','single','ListSize',[100,95]);
% if tf~=0;
% set(TextH6,'String',strcat(list{indx},' noise'));
% used_noise=list{indx};
% end
if tf2~=0
set(TextH7,'String',strcat(list2{indx2},' SNR'));
used_snr=list_temp{indx2};
end


end

function buttionwhite(~,~)
%  list = {'White','Babble','Factory',...                   
% 'Pink','F16','Leopard','M109'};
% [indx,tf] = listdlg('PromptString','Select a noise type:','ListString',list,'SelectionMode','single','ListSize',[100,95]);
 list2 = {'Yes','No'};
list_temp={1,0};
[indx2,tf2] = listdlg('PromptString','Using prewhitening:','ListString',list2,'SelectionMode','single','ListSize',[100,95]);
% if tf~=0;
% set(TextH6,'String',strcat(list{indx},' noise'));
% used_noise=list{indx};
% end
if tf2~=0
set(TextH8,'String',strcat(list2{indx2}));
used_prew=list_temp{indx2};
end


end



function buttionwhite2(~,~)
%  list = {'White','Babble','Factory',...                   
% 'Pink','F16','Leopard','M109'};
% [indx,tf] = listdlg('PromptString','Select a noise type:','ListString',list,'SelectionMode','single','ListSize',[100,95]);
 list2 = {'Yes','No'};
list_temp={1,0};
[indx2,tf2] = listdlg('PromptString','Using prewhitening:','ListString',list2,'SelectionMode','single','ListSize',[100,95]);
% if tf~=0;
% set(TextH6,'String',strcat(list{indx},' noise'));
% used_noise=list{indx};
% end
if tf2~=0
set(TextH9,'String',strcat(list2{indx2}));
used_prew2=list_temp{indx2};
end


end

function noise=get_noise_signals(noiseType)
%     noise=1;
    switch noiseType
        case 'White'
            noise=load('white.mat');
        case 'Babble'
             noise=load('babble.mat');
        case 'Factory'
             
             noise=load('factory2.mat');
        case 'Pink'
             noise=load('pink.mat');
        case 'F16'
            noise=load('f16.mat');
        case 'Leopard'
            noise=load('leopard.mat');
        case 'M109'
            noise=load('m109.mat');
        case 'Gust wind'
            noise=load('gustWind.mat');
        case 'Howl wind'
            noise=load('howlwind.mat');
    end
end

function output_noise=addnoise_strict_snr(sig,input_noise,snr)
input_noise=input_noise(1:length(sig));
noise_std_var=sqrt(10^(-snr/10)*(sig(:)'*sig(:))/(input_noise(:)'*input_noise(:)));
output_noise=noise_std_var*input_noise;
end



end