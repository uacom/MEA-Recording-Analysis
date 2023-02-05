clear
clc

load spkMatrix.mat

spkMat = s_rmShortISI(spkMat, 4, 1000); % see s_rmShortISI.m, to remove short ISIs of <4 ms

%% parameters setting according to the original s_periEventAnalysis.m

fs = 1000; % sampling frequency
tEvent = 1; % time (s) for baseline
binWidth = 50; % window
Time_Interval = [1 1.5]; % the time for estimating whether the unit is excited or inhibited
Thred_Ex = 1.2; % 
Thred_In = -1.2; %

%% estimated the firing rate with the binWidth for the whole time
for i=1:size(spkMat,1)
    for j=1:size(spkMat,2)/binWidth
        count1 = ((j-1)*binWidth)+1;
        count2 = j*binWidth;
        FiringRate(i,j) = sum(spkMat(i,count1:count2))/(binWidth/fs);  % the firing rate for each unit in each window
    end
end

%% calculate the range of firing rate in the baseline for each unit
N_Base = tEvent*fs/binWidth;
for i = 1 : 1 : size(spkMat,1)
    FiringRate_Base_Mean(i) = mean(FiringRate(i,1:N_Base));
    FiringRate_Base_Std(i) = std(FiringRate(i,1:N_Base));
end

% z-scores
for i = 1 : 1 : size(spkMat,1)
    for j = 1 : 1 : size(spkMat,2)/binWidth
        zScores(i,j) = (FiringRate(i,j) - FiringRate_Base_Mean(i))/FiringRate_Base_Std(i);
    end
end

%% deciding the excited unit or inhibited unit in specific time interval
N_Time_Interval = Time_Interval*fs/binWidth;
for i = 1 : 1 : size(spkMat,1)
    FiringRate_Interval(i) = sum(spkMat(i,Time_Interval(1)*fs+1:Time_Interval(2)*fs))/(Time_Interval(2)-Time_Interval(1));
    Thred(i) = (FiringRate_Interval(i)-FiringRate_Base_Mean(i))/FiringRate_Base_Std(i);
end
Units_Ex = find(Thred>Thred_Ex); % numbers of excited units
Units_In = find(Thred<Thred_In); % numbers of inhibited units
