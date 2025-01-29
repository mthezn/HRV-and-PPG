function [tachogram,locs]=pan_tompkins(ECG,fs_ecg,time,varargin)

    %translated and modified by Giulio Steyde 07/08/2023
    %revisitation of the Pan-Tompkins algorithm written at the PheelLab.
    %Credits to Phd Pierluigi Reali, PhD Giulia Tacchino and Professor Maria Bianchi

    % INPUT:
    % ECG: either row or column vector 
    % fs_ecg: integer
    % time: either row or column vector

    % OUTPUT:
    % tacogramma: two column matrix: the first is the time istant in which the R peak was detected, the second is the time distance from the previous peak (in ms)
    % locs: R peaks position (in samples)

    %--------PARAMETERS---------%
    Q_factor = 5;   
    moving_average_samples = round(0.150 * fs_ecg); %150 ms average duration of the QRS complex in adult humans

    %determine if negative peaks can be allowed
    if nargin == 3
        %default is no
        allow_neg = 0; %allow negative r peaks (the electrode configuration was upside-down or something)
    else
        allow_neg = varargin{1};
    end

    % - check input is ok
    if isrow(ECG), ECG = ECG'; end
    if iscolumn(time), time = time'; end
    
    % - band pass filter
    n = 20;             % filter coefficients
    fc = 17;            % central frequency
    Q = Q_factor;       % Q factor
    BW = fc/Q;          % BandWidth
    banda = [(fc-BW/2)*2/fs_ecg (fc+BW/2)*2/fs_ecg];    % normalize wrt Nyquist: fs/2
    b = fir1(n,banda);  % coeff of band pass filter
    y = filtfilt(b,1,ECG); % filter the ECG

    
    % - derivative filter
    H = [1/5 1/10 0 -1/10 -1/5];  % filter coefficients
    y = filter(H,1,y);            % Ora derivo il segnale che prima avevo filtrato passa-banda:
                                  % evidenzio così il QRS che è una variazione veloce del segnale
    MDeriv=length(H)-1;           % PierMOD: Calcoliamo il ritardo imposto dal filtro FIR
    ritDeriv = round(MDeriv/2);   %
    y(1:MDeriv)=0;                % PierMOD: Transitorio di avvio del filtro derivatore:
                                  %          escludiamo i primi M (ordine modello) campioni dai passaggi
                                  %          successivi.
    
    % NB: il filtro derivatore qui usato è di tipo NON RICORSIVO (è un filtro "all-zeros" infatti
    %     a denominatore della FdT abbiamo messo 1 --> FIR) e SIMMETRICO (con
    %     simmetria dispari) e, quindi, a fase lineare. H sono i
    %     coefficienti (in questo caso sono tutti zeri) del filtro.
    % NB2: Nell'help di Matlab si sconsiglia l'uso di "filtfilt" con i
    %      filtri derivatori--> le differenze in termini di ampiezza
    %      relativa tra i picchi si sentono molto rispetto ad utilizzare il
    %      semplice "filter": ho provato a re-implementare il
    %      filtro derivatore calcolando il ritardo generato e riallineando
    %      il segnale ottenuto, anziché utilizzare "filtfilt". La
    %      differenza si vede soprattutto con gli ECG rumorosi (proprio
    %      quelli che fanno funzionare male la funzione originale).

    
    % - squaring + moving average
    y = y.^2;   
    N = moving_average_samples;
    % il filtraggio con il filtro scritto sotto equivale ad una media mobile
    % perchè nel tempo ho la convoluzione tra il segnale e il filtro; quindi
    % moltiplico 32 campioni per 1/N, equivale a fare la media su finestre di
    % 32 campioni che scorrono sul segnale
    H(1:N) = 1/N;   % creo un vettore che ha N termini tutti uguali a 1/N
    y = filtfilt(H,1,y);    % filtro il segnale
    y = y./max(y); 

    
    % - find QRS complex
	fCardiacaMax = 300; %bpm (maximum possible cardiac frequency)
	fCardiacaMaxHz = fCardiacaMax/60; %Hz
    ncamp_min = round(fs_ecg/fCardiacaMaxHz);
    [~,locs] = findpeaks(y,'MinPeakProminence',.1,'MinPeakDistance',ncamp_min); 
    
    
    % - correct delay induced by filter
    locs = locs-ritDeriv; 
    
    threshold = round(ncamp_min/2);  %identify threshold to look for r peaks in the baseline ECG
    for z=1:length(locs)
        
        % lower bound
        minc=locs(z)-threshold;
        if minc<=0, minc=1; end   
        
        % upper bound
        maxc=locs(z)+threshold;      
        if maxc>length(ECG), maxc=length(ECG); end 
        
        % find highest peak in the window
        if allow_neg
            [~,pos] = max( abs( ECG(minc:maxc)-mean(ECG(minc:maxc)) ) );%can be negative
        else
            [~,pos] = max(ECG(minc:maxc)-mean(ECG(minc:maxc)));%cannot be negative
        end
        locs(z) = pos+minc-1;
        
    end

    % - TACHOGRAM
    tachogram = zeros(length(locs)-1,1); %inizializzo il vettore
    for i = 2:length(locs)
        tachogram(i-1) = (locs(i)-locs(i-1))/fs_ecg;
    end
    tachogram = [time(locs(2:end))' tachogram];

    %to milliseconds
    tachogram(:,2) = 1000*tachogram(:,2); 

    

end

