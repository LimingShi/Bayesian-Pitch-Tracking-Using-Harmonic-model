function Obj=offlineprocessData_orig(Obj, event)
% tic
    % this is a so-called callback function. See more here
    % https://se.mathworks.com/help/matlab/matlab_prog/timer-callback-functions.html
    x= Obj.UserData{25};    
    x_clean= Obj.UserData{26};  
    time_indx=Obj.UserData{2};
    frame=x(time_indx);
    frame_clean=x_clean(time_indx);
    [f0Estimates,order,vad_prob]=Obj.UserData{4}(frame,1);
    if vad_prob<=.5
        f0Estimates=nan;
    end
    Obj.UserData{6}(Obj.UserData{1})=f0Estimates*Obj.UserData{5};
    win_input=frame_clean.*gausswin(length(frame));
    pow_spec=periodgram_shi(win_input,256);
%     Obj.UserData{7}(Obj.UserData{1},Obj.UserData{6}(Obj.UserData{1}),order,pow_spec);
%     plot(Obj.UserData{1},f0Estimates*Obj.UserData{5},'.');hold on;

%     Obj.UserData{8}{1}.XData(time_indx)=time_indx/Obj.UserData{5};
%     Obj.UserData{8}{1}.YData(time_indx)=x(time_indx);
    if Obj.UserData{1}==1
        sample_vec=[1:Obj.UserData{11}*Obj.UserData{5}]/Obj.UserData{5};
        Obj.UserData{12}{1}=plot(Obj.UserData{9},sample_vec,nan(length(sample_vec),1));  
        xlabel(Obj.UserData{9},'Time [s]')
        ylabel(Obj.UserData{9},'Time series');
        
        Obj.UserData{12}{2}(1)=imagesc(Obj.UserData{8},Obj.UserData{20}, Obj.UserData{21},nan(Obj.UserData{22},length(Obj.UserData{20})));
        ylim(Obj.UserData{8},[0,Obj.UserData{23}])
        xlim(Obj.UserData{8},[0 Obj.UserData{11}])
        axis(Obj.UserData{8},'xy');
        hold on;
        Obj.UserData{12}{2}(2)=plot(Obj.UserData{15},Obj.UserData{20},nan(length(Obj.UserData{20}),1),'r-', 'linewidth',2);
          ylim(Obj.UserData{15},[0,Obj.UserData{23}])
        xlim(Obj.UserData{15},[0 Obj.UserData{11}])
        axis(Obj.UserData{15},'off');
        ylabel(Obj.UserData{8},'Frequency')
        
        
         Obj.UserData{12}{3}=plot(Obj.UserData{10},Obj.UserData{20},nan(length(Obj.UserData{20}),1),'r-', 'linewidth',2);
        ylim(Obj.UserData{10},[0,1])
        xlim(Obj.UserData{10},[0 Obj.UserData{11}])
         ylabel(Obj.UserData{10},'Voicing Probability');    
         
                 
         Obj.UserData{12}{4}=plot(Obj.UserData{17},Obj.UserData{20},nan(length(Obj.UserData{20}),1),'r.', 'linewidth',2);
        ylim(Obj.UserData{17},[0,Obj.UserData{24}])
        xlim(Obj.UserData{17},[0 Obj.UserData{11}])
         ylabel(Obj.UserData{17},'Order');  
        
    end
%     soundsc(x(time_indx))
     Obj.UserData{12}{1}.YData(time_indx)=x(time_indx);
     xlim(Obj.UserData{9},[0,Obj.UserData{11}]);
     Obj.UserData{12}{2}(1).CData(:,Obj.UserData{1})=10*log10(abs(pow_spec(1:Obj.UserData{22})));
     Obj.UserData{12}{2}(2).YData(Obj.UserData{1})=[ Obj.UserData{6}(Obj.UserData{1})];
     xlim(Obj.UserData{8},[0,Obj.UserData{11}]);
     
     Obj.UserData{12}{3}.YData(Obj.UserData{1})=vad_prob;
     xlim(Obj.UserData{10},[0 Obj.UserData{11}]);
     
     Obj.UserData{12}{4}.YData(Obj.UserData{1})=order;
      xlim(Obj.UserData{17},[0 Obj.UserData{11}]);
    
     Obj.UserData{1}=Obj.UserData{1}+1;
    Obj.UserData{2}=Obj.UserData{2}+Obj.UserData{3};
    
    

% toc
end
function pow=periodgram_shi(x,nfft,fs)
if nargin<3
    fs=1;
end
x=x(:);
data_len=length(x);
% matr=fft_shi(data_len,nfft);
%normalized rectangular window consistent with energy
f=fft(x,nfft)/sqrt(data_len);
%normalized rectangular window consistent with matlab
% f=matr*x/sqrt(data_len);
pow=abs(f).^2;
end