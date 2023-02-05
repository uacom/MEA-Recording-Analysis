clear all
clc
close all
% NOTE, when running this script, do not save and overwrite the myData2.mat
% data file!!!!!!!
load ('C:\Users\Shenfeng Qiu\Documents\MATLAB\myData2.mat')
% load ('C:\Users\Shenfeng Qiu\Documents\MATLAB\myData.mat')


%% select data
TimePoints = [5.55,15.32,32.14,46.55,60.25,76.45]; % event points
Ta = 2.5; 
Tb = 3.5;
Fs = 1000;

FP07_sel = LFP_Select_Qiu(FP07,TimePoints,Ta,Tb,Fs);
FP08_sel = LFP_Select_Qiu(FP08,TimePoints,Ta,Tb,Fs);
SPK02b = SPK_Select_Qiu(SPK02b,TimePoints,Ta,Tb,Fs);
SPK04a = SPK_Select_Qiu(SPK04a,TimePoints,Ta,Tb,Fs);
SPK09a = SPK_Select_Qiu(SPK09a,TimePoints,Ta,Tb,Fs);


%% LFP signal estimation
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

% main estimation
[CohValue,CauValue,Cohf,Cauf] = LFP_Estimation_Qiu_MultiTrials(FP07_sel,FP08_sel,ParaCoh,ParaCau);

% plot figures

% coherence, figure(1)
Coh_mean = mean(CohValue,1);
Coh_sem = std(CohValue,0,1)/sqrt(size(CohValue,1));
ShadedErrorBar_Qiu(Coh_mean,Coh_sem,Cohf,'Coherence')

% Casuality: Xij means j casual i, figure(2)
Cau_mean(1,:) = mean(squeeze(CauValue(:,2,1,:)),1);
Cau_mean(2,:) = mean(squeeze(CauValue(:,1,2,:)),1);
Cau_sem(1,:) = std(squeeze(CauValue(:,2,1,:)),0,1)/sqrt(size(squeeze(CauValue(:,2,1,:)),1));
Cau_sem(2,:) = std(squeeze(CauValue(:,1,2,:)),0,1)/sqrt(size(squeeze(CauValue(:,1,2,:)),1));
ShadedErrorBar_Qiu(Cau_mean,Cau_sem,Cauf,'Casuality')


%% For the whole time series, figures(3, 4,5,6,7,8...) 
% Phase locking between LFP and Spike
% 
fpass = [1 13];
Fs = 1000;
phs_res = 16;
for i = 1 : 1 : length(FP07_sel)
    figure
    [Result{i}] = LFP_SPK_Phase(FP07_sel{i},SPK02b{i},Fs,fpass,phs_res);
end

%% spike and LFP lagged PLV calculation. Comment this block if not needed.
clear all
clc

load ('C:\Users\Shenfeng Qiu\Documents\MATLAB\myData2.mat')

TimePoints = [5.55,15.32,32.14,46.55,60.25,76.45]; % event points
Ta = 1.5; 
Tb = 2.5;
Fs = 1000;
fpass = [1 13];
TimeRang = [-150 150]; % input as miliseconds
fp = FP07;
spk = SPK02b;

[Lagged_PLV] = LFP_SPK_LaggedPhase(fp,spk,Fs,fpass,TimePoints,Ta,Tb,TimeRang); % return Lagged_PLV as a n*m matrix

%% sort in descending order for the largest PLV value, and plot, figure(9)
[~,sortIndex] = sort(max(Lagged_PLV'), 'descend');
Lagged_PLV_sort = Lagged_PLV(sortIndex,:);
figure
imagesc(Lagged_PLV_sort); colormap jet;

autoArrangeFigures(0, 0, 2)




