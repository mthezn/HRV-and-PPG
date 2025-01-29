clear all
close all 
clc 

meanTaco = []; %column 1: HRV reference, column 2: HRV evaluated from ppg
stdTaco = []; %column 1: HRV reference, column 2: HRV evaluated from ppg
meanDelay = []; % mean of the delays between R peaks and Systolic peaks
stdDelay = []; %std of the delays

for j = 5:5
    load("./DATA//DATA/S0" + j +".mat")    
    fs = dati.fs;
    ecg = dati.ECG;  
    t = dati.Time; 
    N = length(ecg); 

    d = 1
    z = 1 %indexes for the previous arrays
    %% ecg spectra to find noise 
    % this section is about the analysis of the fourier transform of every
    % ecg in order to find if there is some noise superimposed
    
    Xecg = fft_agile(ecg) ;
    f = 0:fs/N:fs/2; 
    
    
    %Plot of the non filtered signal spectrum
    % figure; 
    % plot(f,Xecg) 
    % xlim([0,100] )  
    % xlabel("Frequency(f)")
    % ylabel("|X(f)|") 
    % title("Subject " + j + " not filtered ecg spectrum")
    %% filtering  
    %this section is about filtering our signal inside the typical
    %frequencies of an ecg (0.05-150Hz) 
 
    %High pass filter 
    filter_order = 4;
    fcut = 0.05;       % cut-off frequency [Hz]
    wc = fcut/(fs/2); % NORMALIZED cut-off frequency
    [b,a] = butter(filter_order,wc,'high'); 
    ecg_filt = filter(b,a,ecg) ;  
    
    
    %low pass filter
    fcut = 150;  
    wc = fcut/(fs/2); % NORMALIZED cut-off frequency
    [b,a] = butter(filter_order,wc,'low'); 
    ecg_filt = filter(b,a,ecg_filt) ;
    

    %veryfing the effective filtering of the signal 
   
    Xecg_filt = fft_agile(ecg_filt) ;
    %Plot of the filtered signal  
    % figure;
    % plot(f,Xecg_filt);
    % xlim([0,200]) 
    % xlabel("Frequency(f)")
    % ylabel("|X(f)|") 
    % title("Subject " + j + " filtered ecg spectrum" )
    % 
    

    %take a window of the signals to better visualize it
    step = fs * 60 * 5; %we will take a five minutes window 
    
    Nstart =  1 ;
   
    Nstop = Nstart + step;  
    
  
    
    while (Nstop <= length(dati.ECG)) 
    
    ecg = ecg_filt(Nstart:Nstop);
    %ecg = ecg - mean(ecg); 
    ppg = dati.PPG(Nstart:Nstop);
    %ppg = ppg - mean(ppg); 
    %ppg = normalize(ppg); %normalizing the signal will be usufeul to use the same algorithm parameters for every signal
    t = dati.Time(Nstart:Nstop); 
    N = length(ecg);
    
    Minstart = round((Nstart*(1/fs))/60) 
    Minstop = round((Nstop*(1/fs))/60)
    max_press = 200 
    min_press = 20 

    %ppg_rescaled = rescale(ppg,min_press,max_press); %we rescale the ppg signal in order to make it inside a range of values which is similar to the ones 
                                                     % of arterial blood
                                                     % pressure so the
                                                     % functions can work 
    
    %A plot just to visualize teh two signals one over the other
    % figure;
    % plot(t,ecg - mean(ecg)) 
    % hold on 
    % plot(t,ppg -mean(ppg)) 
    % xlabel("Time(s)")
    % hold off  
    % title("Subject " + j + " ecg and ppg from minute " + (Nstart*(1/fs))/60 + " to " + (Nstop*(1/fs))/60 )
    % 
    
    %% estimate HRV from PPG (findpeaks) 
    
 [pks,pos] = findpeaks(ppg, "MinPeakProminence",0.8,'MinPeakDistance',700) ;  %freqMax assumed hypotetically 200bpm = 3,3 beats/second = 1 beat 
                                                                                %every 0,333 s -> 680 samples

    posTime = pos*1/fs  %trasformo in tempo gli indici di posizione 
  
    RR_ppg = diff(posTime)*1000;    

    %t_RR_ecg = t(Rpeaks);  

    %taco_time(1) = [];  
    %taco_time(end) = []
    t_RR_ppg = t(pos); 
    t_RR_ppg(1) = []
    %find onsets   

    %resample blood pressure to 125Hz (in order to use wabp)
   
   
   
    
    %% estimate HRV from ecg 
    
    [tachogram,locs]=pan_tompkins(ecg,fs,t,0)  ; 
    
    
    taco_time = tachogram(:,1);
    taco_distance = tachogram(:,2); 
    
    %in order to synchronize the two signals and start from the same beat
    if(t_RR_ppg(1) < taco_time(1)) 
        t_RR_ppg = t_RR_ppg(2:end);  
        RR_ppg = RR_ppg(2:end); 

    end

     
     %statisticalValues 
   %the column is the number of the subject
   % for every couple of values the first one is the reference one 
   %and the secondo is the one from the hrv evaluated 
   t_RR_ppg_noOut = t_RR_ppg
   [RR_ppg_noOut,posRem] = rmoutliers(RR_ppg,"mean") 
   
   
    meanTaco(d,j) = mean(taco_distance); 
    t_RR_ppg_noOut(posRem) = []

    d = d + 1
    meanTaco(d,j) = mean(RR_ppg);
    stdTaco(z,j) = std(taco_distance);
    z = z + 1
    
    stdTaco(z,j) = std(RR_ppg);
    
    
    %% plotting 
     figure; 
    % plot(t,ecg)
    hold on 
     plot(t,ppg)
      plot(t(pos),ppg(pos),"Marker","o")
    % plot(t(pos),ppg(pos),"Marker","*") 
    title("Subject " + j + " check peaks from minute " + Minstart + " to " + Minstop  )
    hold off 
    xlabel("t(s)")
    ylabel("a.u")
    % 


    figure
   
    plot(taco_time,taco_distance)
    hold on 
    plot(t_RR_ppg,RR_ppg);
    plot(t_RR_ppg_noOut,RR_ppg_noOut)
    
    xlabel('time [s]')
    ylabel('\Delta time [ms]') 
    hold off 
    title("Subject " + j + " Tachogram from minute " + Minstart + " to minute " + Minstop )
    legend("Reference tachogram", "Tachogram from ppg","Tachogram without outliers") 

    ylim([300 1500])
    %% PSD of the tachograms 
    %resampling to 3 Hz in order to have a better explanation of the low
    %frequency components
    
    f_rs2 = 3; 
    
    RR_ppg_r_noOut = interp1(t_RR_ppg_noOut,RR_ppg_noOut,(t_RR_ppg_noOut(1):1/f_rs2:t_RR_ppg_noOut(end))); 
    t_RR_ppg_r_noOut = t_RR_ppg_noOut(1): 1/f_rs2 : t_RR_ppg_noOut(end); 
    
    RR_ppg_r = interp1(t_RR_ppg,RR_ppg,(t_RR_ppg(1):1/f_rs2:t_RR_ppg(end))); 
    t_RR_ppg_r = t_RR_ppg(1): 1/f_rs2 : t_RR_ppg(end); 
    
    taco_distance_r = interp1(taco_time,taco_distance,(taco_time(1):1/f_rs2:taco_time(end))); 
    taco_time_r = t_RR_ppg(1): 1/f_rs2 : t_RR_ppg(end); 
    
    RR_ppg_d_noOut = detrend(RR_ppg_r_noOut)
    RR_ppg_d = detrend(RR_ppg_r);
    taco_distance_d = detrend(taco_distance_r); 
    
    %compute the periodogram of tacogram from ppg with welch 
    %4 windows 
    k = 4 ;
    M = ceil(length(RR_ppg_d)/k + 1) ;
    noverlap = ceil(M/2);
    [PSD_ppg,f_ppg] = pwelch(RR_ppg_d/1000,hamming(M),noverlap,[],f_rs2) ;
    
    %compute the periodogram of tacogram from ecg with welch 
    %4 windows 
    k = 4 ;
    M = ceil(length(taco_distance_d)/k) ;
    noverlap = ceil(M/2);
    [PSD_ecg,f_ecg] = pwelch(taco_distance_d/1000,hamming(M),noverlap,[],f_rs2) ;
     
    k = 4 ;
    M = ceil(length(RR_ppg_noOut)/k + 1) ;
    noverlap = ceil(M/2);
    [PSD_ppg_noOut,f_ppg_noOut] = pwelch(RR_ppg_d_noOut/1000,hamming(M),noverlap,[],f_rs2) ;
    


    figure; 
    plot(f_ppg,PSD_ppg)
    hold on 
    plot(f_ecg,PSD_ecg) 
    plot(f_ppg_noOut,PSD_ppg_noOut)
    legend("PSD of tachogram from PPG", "PSD of tachogram from ECG","PSD of tachogram from PPG no outliers")
    title("HRV PSD subject " + j + " from minute " + Minstart +  " to " + Minstop)
    xlabel('Frequency [Hz]')
    ylabel('s^{2}/Hz')
    hold off 

    %% studying the delay between R peaks and systolic peaks
    % pos_picchi_ecg = locs/fs;
    % pos_picchi_ppg = r(:,1) /f_rs1;
    % 
    % %synchronizing the two signals to start from the same beat
    % if(pos_picchi_ppg(1) < pos_picchi_ecg(1)) 
    %     pos_picchi_ppg = pos_picchi_ppg(2:end);
    % end
    % index = min(length(pos_picchi_ecg),length(pos_picchi_ppg));
    % delay = ( pos_picchi_ppg(1:index) - pos_picchi_ecg(1:index) ) ;
    % 
    % figure 
    % plot(delay)
    % hold on 
    % yline(mean(delay))
    % title("Delay between R peak and systolic peak for subject " + j + " from minute " + Minstart + " to " + Minstop) 
    % ylabel("t(s)")
    % xlabel("#beat") 
    % hold off 
    % 
    % %statistical values 
    % meanDelay(j) = mean(delay);
    % stdDelay(j) = std(delay);


    %% Estimate respiratory signal from PPG
   
    Resp = dati.RESP(Nstart:Nstop);
    Resp = detrend(Resp);
    t = (1:length(Resp))./fs;
    
    fcut = 1.7;
    fs_new = 4;
    
    [PSD_Resp, f_resp] = psd(Resp,fs,fcut,fs_new);

    figure
    plot(f_resp,PSD_Resp)

    xlim([0 fs_new/2])

    
    
    % define the filter 
    n_order = 4;      
    fcut = 0.4;       
    
    wn = fcut/(fs/2); % NORMALIZED cut-off frequency
    
    % plot butterworth filter (IIR HP) characteristics 
    
    [b,a] = butter(n_order,wn,'low');
    %figure; freqz(b,a,1024,fs);
    
    % use the filter
    
    Resp_from_ppg = filtfilt(b,a,dati.PPG-mean(dati.PPG)); 
    Resp_from_ppg= detrend(Resp_from_ppg);

    %delete the delay
    Resp_orig = dati.RESP-mean(dati.RESP);
    [c,lags] = xcorr(Resp_from_ppg(Nstart:Nstop),Resp_orig(Nstart:Nstop)); 
    c = c/max(c);
    [max_c,pos_max_c] = max(c);
    t_c = lags(pos_max_c);
    Resp_from_ppg1 = Resp_from_ppg(Nstart+abs(t_c):Nstop);
    Resp_from_ppg1 = detrend(Resp_from_ppg1);

    %methos 2: envelope
    [env_upper, env_lover] = envelope(dati.PPG(Nstart:Nstop-abs(t_c))-mean(dati.PPG(Nstart:Nstop-abs(t_c))),fs,'peaks');
    Resp_env = env_lover - mean(env_lover);
    Resp_env=detrend(Resp_env);
    

    %plot PSD of Resp from ppg vs PSD of Resp from envelope vs PSD of Resp
    %as reference
    
    fcut = 1.7;
    fs_new = 4;
    
    % PSD resp from PPG
    [PSD_resp_ppg,f_resp_ppg]=psd(Resp_from_ppg1,fs,fcut,fs_new);

    figure
    plot(f_resp_ppg,PSD_resp_ppg)
    hold on

    %PSD resp from envelope
    [PSD_resp_env,f_resp_env]=psd(Resp_env,fs,fcut,fs_new);
    plot(f_resp_env,PSD_resp_env)
    

    %PSD resp as reference
    [PSD_Resp, f_resp] = psd(Resp_orig(Nstart:Nstop-abs(t_c))-mean(Resp_orig(Nstart:Nstop-abs(t_c))),fs,fcut,fs_new);
    plot(f_resp,PSD_Resp)
    xlim([0 fs_new/2])
    title("RESP PSD for subject " + j + " from minute " + Minstart + " to minute " + Minstop )
    xlabel('Frequency [Hz]'); ylabel('a.u.^{2}/Hz')
    legend('PSD resp from PPG','PSD resp from envelope', 'PSD resp as reference')
    hold off

    %Original respiratory signal vs respiratory signal extracted from PPG:
    %window of 60 s
    figure
    plot(dati.Time(Nstart:Nstop-abs(t_c)),Resp_from_ppg1*5)
    hold on
    plot(dati.Time(Nstart:Nstop-abs(t_c)),Resp_orig(Nstart:Nstop-abs(t_c))-mean(Resp_orig(Nstart:Nstop-abs(t_c))))
    plot(dati.Time(Nstart:Nstop-abs(t_c)),Resp_env)
    axis tight
    xlabel('t(s)');
    ylabel('a.u.');
    title("RESP for subject " + j + " from minute " + Minstart + " to minute " + Minstop ); 
    legend('RESP from filtered PPG','RESP as reference','RESP from envelope')
    hold off

    %% Consider the whole ecg signal
    
    % 
    % ecg = ecg_filt;
    % ecg = ecg - mean(ecg); 
    % t = dati.Time; 
    % 
    % % estimate HRV from ecg 
    % 
    % [tachogram,locs]=pan_tompkins(ecg,fs,t,0)  ; 
    % 
    % %resampling to 4 Hz
    % f_rs = 4; 
    % 
    % taco_time = tachogram(:,1);
    % taco_distance = tachogram(:,2); 
    % taco_distance_r = interp1(taco_time,taco_distance,(taco_time(1):1/f_rs:taco_time(end))); 
    % taco_time_r = t_RR_ppg(1): 1/f_rs : t_RR_ppg(end); 
    % 
    % RR_ppg_d = detrend(RR_ppg_r);
    % taco_distance_d = detrend(taco_distance_r);  

    %compute the periodogram of tacogram from ecg with welch 
    %4 windows 
    k = 4 ;
    M = ceil(length(taco_distance_d)/k) ;
    noverlap = ceil(M/2);
    [PSD_ecg,f_ecg] = pwelch(taco_distance_d/1000,hamming(M),noverlap,[],f_rs2) ;
    
    %compute the PSD for the whole respiratory signal estimated from PPG
    fcut = 1.7;
    fs_new = 4;
    [PSD_resp_ppg,f_resp_ppg]=psd(Resp_from_ppg,fs,fcut,fs_new);
   


    % single plot HRV and RESP PSDs 

    figure
    subplot(211)
    plot(f_ecg,PSD_ecg)
    title("HRV PSD subject " + j + " from minute " + Minstart+ " to " + Minstop) ;
    xlim([0 0.5])
    subplot(212)
    plot(f_resp_ppg,PSD_resp_ppg)
    title("RESP PSD " + j + " from minute " + Minstart + " to " + Minstop); 
    xlim([0 0.5])
    
    Nstart = Nstop
   
    Nstop = Nstart + step 
    pause
    end 

end 
    


    
%% funzioni 
function  X = fft_agile(x)

X = fft(x);
N = length(x);
X = abs(X); %modulo dei coefficenti(ampiezza dello spettro),normalizzata per la
             %ampiezza dello spettro

X = 2*X(1:N/2 + 1); %spettro simmetrico quindi ne considero solo metà 
end  


function  [PSD_x,f] = psd(x,fs,fcut,fs_new)
% downsample RESP to avoid aliasing
%low-pass before downsampling to avoid aliasing perchè il segnale potrebbe
%avere delle frequenze più alte di 2fs
b = fir1(120,fcut/(fs/2),'low');
%figure,freqz(b,1);

x_lp = filtfilt(b,1,x);

%now downsample
x_ds = downsample(x_lp,fs/fs_new); %basically, take one sample every fs/fs_new
t = (1:length(x))./fs;
t_ds = downsample(t,fs/fs_new);

% %check it is ok
% figure
% plot(t,x)
% hold on
% plot(t_ds,x_ds)
% xlabel('time [s]'), ylabel('a.u.')
% hold off

% with Welch
%4 windows - Welch
k = 4;
L = ceil((2*length(x_ds))/(k+1));
noverlap = L/2;
[PSD_x,f] = pwelch(x_ds,hamming(L),noverlap,[],fs_new); %4 minutes

end



