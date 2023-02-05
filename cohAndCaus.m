function [cau1D,cau2D,Cauf] = cohAndCaus(lfp1,lfp2,t_point,t_dura,Fs,corc)
%% Computes coherence and causality between two different lfp signals at various timepoints
%
% INPUTS
% lfp1, first lfp trace (hippocampus field potential, FP08)
% lfp2, second lfp trace (PF field potential, FP07)
% t_point = timepoint
% t_dura = peri-timepoint length, write as vector (e.g. [1.5,2.5])
% Fs, sampling frequency (Hz)
% corc = 1 for coherence output, 2 for causality output, 3 for both
%
% OUTPUTS
% cau1D = Causality in the direction of lfp1 on lfp2
% cau2D = Causality in the direction of lfp2 on lfp1
% Cauf = Frequencies (x-values) for causality graph


%%% lfp1 = FP08;
%%% lfp2 = FP07;
%%% t_point = 115;
%%% t_duration = [1.5,2.5];
%%% Fs = 1000;    
%%% corc = 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THIS WORKS:
% load ('myData2.mat')
%%% [caus1,caus2,cf] = cohAndCaus(FP08,FP07,115,[1.5,2.5],1000,1);
%%% [caus1,caus2,cf] = cohAndCaus(FP08,FP07,115,[1.5,2.5],1000,2);
%%% [caus1,caus2,cf] = cohAndCaus(FP08,FP07,115,[1.5,2.5],1000,3);

%% Filter the signals
% fNQ=Fs/2;
% Wn=[4 12]/fNQ;
% filtD1=fir1(100,Wn); %creates the filter
% 
% flfp3 = filtfilt(filtD1,1,transpose(FP09));
% flfp3 = transpose(flfp3);
% 
% lfp1 = lfp1+FP09;
% lfp2(25:343216) = lfp2(25:343216) + FP09(1:343192);
%% Parameter settings for coherence and causality

% Parameters setting for coherence
ParaCoh.tapers =[3 5];
ParaCoh.pad = 0;
ParaCoh.Fs = Fs;
ParaCoh.fpass = [0 80];

% Parameters setting for causality
ParaCau.regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
ParaCau.icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

ParaCau.morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
ParaCau.momax     = 30;     % maximum model order for model order estimation

ParaCau.acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

ParaCau.alpha     = 0.05;   % significance level for significance test
ParaCau.mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

ParaCau.Fs        = Fs;    % sample rate (Hz)
ParaCau.fres      = 1000;     % frequency resolution (empty for automatic calculation)
ParaCau.fpass = [0 80];   % the range of frequency for causality

%% Calculate coherence and causality

% Snip trace to the desired timepoint and duration
ctp = t_point * Fs;
ctd = t_dura .* Fs;
slfp1 = lfp1(ctp-ctd(1):ctp+ctd(2));
slfp2 = lfp2(ctp-ctd(1):ctp+ctd(2));

if corc == 1 || corc == 3
    if corc == 3
        figure
    end
% Calculate coherence
    [C,~,~,~,~,f]=coherencyc(slfp1,slfp2,ParaCoh);
    CohValue = C;
    Cohf = f;


    plot(Cohf, CohValue)
    xlabel('Frequency (Hz)')
    ylabel('Coherence')
end

%% calculate the causality
if corc == 2 || corc == 3
    if corc == 3
        figure
    end
    data(1,:) = slfp1;
    data(2,:) = slfp2;
    
    [Results] = Cal_FreGC_Cui(data,ParaCau);
    CauValue = Results.FreGC;
    Cauf = Results.FreGC_Fres;

    cau1D = squeeze(CauValue(2,1,:));
    cau2D = squeeze(CauValue(1,2,:));


    plot(Cauf,cau1D)
    hold on
    plot(Cauf,cau2D)
    xlabel('Frequency (Hz)')
    ylabel('Granger Causality')
    hold off
end

end