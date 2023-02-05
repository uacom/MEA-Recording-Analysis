function [] = LFP_analysis(LFP_PFC, LFP_HPC, SPK, Ta, Tb, timePoints)
% raw/filtered LFP data with raster plot and phase histogram (Tony Nehme 09112018)

%%% load('C:\Users\Shenfeng Qiu\Documents\MATLAB\File2_053118.mat');
%%% LFP_analysis(FP_PFC(:,1), FP_HPC(:,1), unit11a, 2.5, 2.5,[5.55,15.32,32.14,46.55,60.25,76.45]);
    %%%f_range is freq range for bandpass filter
    %%%Ta is time before timestamp
    %%%Tb is time after timestamp
    %%% for freq range definition: delta 0.5-4 Hz, theta 4-12 Hz, beta
    %%% 12-30; gamma 30-100 Hz. 
%     vars for debugging
%     
%     LFP_PFC = FP_PFC(:,1);
%     LFP_HPC = FP_HPC(:,1);
%     SPK = unit11a;
%     Ta = 2.5;
%     Tb = 2.5;
%     timePoints = [5.55,15.32,32.14,46.55,60.25,76.45];

% % % % % % % % % % % % % % % % % % % % % % 

%%select data
Fs = 1000;
selLFP_PFC = LFP_Select_Qiu(LFP_PFC,timePoints,Ta,Tb,Fs);
selLFP_HPC = LFP_Select_Qiu(LFP_HPC,timePoints,Ta,Tb,Fs);
selSPK = SPK_Select_Qiu(SPK,timePoints,Ta,Tb,Fs);

%% filtering LFP data

% theta filter
dt = 1/Fs;
fNQ=Fs/2;
Wn=[4 12]/fNQ;
filtD1=fir1(100,Wn); %creates the filter

LFP_PFC_tfilt=filtfilt(filtD1,1,transpose(LFP_PFC));
LFP_PFC_tfilt=transpose(LFP_PFC_tfilt);

LFP_HPC_tfilt=filtfilt(filtD1,1,transpose(LFP_HPC));
LFP_HPC_tfilt=transpose(LFP_HPC_tfilt);

for i = 1:size(selLFP_PFC, 2) %applies filter to each time point
    selLFP_PFC_tfilt{1, i}=filtfilt(filtD1,1,transpose(selLFP_PFC{1,i}));
    selLFP_PFC_tfilt{1, i}=transpose(selLFP_PFC_tfilt{1, i});
end
for i = 1:size(selLFP_HPC, 2)
    selLFP_HPC_tfilt{1, i}=filtfilt(filtD1,1,transpose(selLFP_HPC{1,i}));
    selLFP_HPC_tfilt{1, i}=transpose(selLFP_HPC_tfilt{1, i});
end

% delta filter
Wn=[0.5 4]/fNQ;
filtD1=fir1(100,Wn); %creates the filter

LFP_PFC_dfilt=filtfilt(filtD1,1,transpose(LFP_PFC));
LFP_PFC_dfilt=transpose(LFP_PFC_dfilt);

LFP_HPC_dfilt=filtfilt(filtD1,1,transpose(LFP_HPC));
LFP_HPC_dfilt=transpose(LFP_HPC_dfilt);

for i = 1:size(selLFP_PFC, 2) %applies filter to each time point
    selLFP_PFC_dfilt{1, i}=filtfilt(filtD1,1,transpose(selLFP_PFC{1,i}));
    selLFP_PFC_dfilt{1, i}=transpose(selLFP_PFC_dfilt{1, i});
end
for i = 1:size(selLFP_HPC, 2)
    selLFP_HPC_dfilt{1, i}=filtfilt(filtD1,1,transpose(selLFP_HPC{1,i}));
    selLFP_HPC_dfilt{1, i}=transpose(selLFP_HPC_dfilt{1, i});
end

% 100 filter
Wn=[0.5 100]/fNQ;
filtD1=fir1(100,Wn); %creates the filter
LFP_PFC_100filt=filtfilt(filtD1,1,transpose(LFP_PFC));
LFP_PFC_100filt=transpose(LFP_PFC_100filt);

LFP_HPC_100filt=filtfilt(filtD1,1,transpose(LFP_HPC));
LFP_HPC_100filt=transpose(LFP_HPC_100filt);

for i = 1:size(selLFP_PFC, 2) %applies filter to each time point
    selLFP_PFC_100filt{1, i}=filtfilt(filtD1,1,transpose(selLFP_PFC{1,i}));
    selLFP_PFC_100filt{1, i}=transpose(selLFP_PFC_100filt{1, i});
end
for i = 1:size(selLFP_HPC, 2)
    selLFP_HPC_100filt{1, i}=filtfilt(filtD1,1,transpose(selLFP_HPC{1,i}));
    selLFP_HPC_100filt{1, i}=transpose(selLFP_HPC_100filt{1, i});
end
%% graphing raw LFP with filtered LFP and spike raster plot

%% theta
% figure(1)
tz = dt:dt:length(LFP_HPC)*dt;
figure
for i = 1:length(timePoints)
    t = (timePoints(i)-Ta+dt:dt:timePoints(i)+Tb);
    
    fphit=angle(hilbert(LFP_PFC_tfilt));
    sel_phase = find(fphit >= -pi/4); % SQ comment, define the start phase to be shaded.
    % sel_phase = find(fphit >= -pi/4);
    sel_phase = sel_phase(find(fphit(sel_phase) <= pi/4)); % SQ comment, define the end phase to be shaded.
    % sel_phase = sel_phase(find(fphit(sel_phase) <= pi/4)); 
    
    subplot(2,length(timePoints),i)
        selsel_phase = sel_phase(find(tz(sel_phase) >= timePoints(i)-Ta));
        selsel_phase = selsel_phase(find(tz(selsel_phase) <= timePoints(i)+Tb));
        ii = 1;
        m = 1;
            for k=2:length(selsel_phase)
                if(selsel_phase(k)-selsel_phase(k-1) == 1)
                    sel3_phase{ii}(m) = selsel_phase(k-1);
                    sel3_phase{ii}(m+1) = selsel_phase(k);
                    m = m+1;
                else
                    ii = ii+1;
                    m = 1;
                end
            end 
            
        for n=1:length(sel3_phase)
            if isempty(sel3_phase{n}) == 0
                patch([tz(min(sel3_phase{n})) tz(min(sel3_phase{n})), tz(max(sel3_phase{n})) tz(max(sel3_phase{n}))],[-1 1 1 -1],[0.8 0.8 0.8],'EdgeColor','none','HandleVisibility','off');
            end
        end
        hold on
        plot(t,selLFP_PFC{i}, 'k','LineWidth',0.4)
        plot(t,selLFP_PFC_tfilt{i}, 'r', 'LineWidth', 1)
        for j=1:length(selSPK{i})
            line([selSPK{i}(j)+timePoints(i)-Ta selSPK{i}(j)+timePoints(i)-Ta], [0.8 1])
        end     
        hold off
        title('PFC')
        xlabel('Time (s)')
        ylabel('Voltage (mV)')
        axis([timePoints(i)-Ta timePoints(i)+Tb -1 1])
        legend('Raw data','Theta-filtered data','Location','southwest')
        
        clear fphit
        clear sel_phase
        clear selsel_phase
        clear sel3_phase
    
    fphit=angle(hilbert(LFP_HPC_tfilt));
    sel_phase = find(fphit >= -pi/4);
    sel_phase = sel_phase(find(fphit(sel_phase) <= pi/4));    
    
    subplot(2,length(timePoints),i+length(timePoints))
        selsel_phase = sel_phase(find(tz(sel_phase) >= timePoints(i)-Ta));
        selsel_phase = selsel_phase(find(tz(selsel_phase) <= timePoints(i)+Tb));
        ii = 1;
        m = 1;
            for k=2:length(selsel_phase)
                if(selsel_phase(k)-selsel_phase(k-1) == 1)
                    sel3_phase{ii}(m) = selsel_phase(k-1);
                    sel3_phase{ii}(m+1) = selsel_phase(k);
                    m = m+1;
                else
                    ii = ii+1;
                    m = 1;
                end
            end 
            
        for n=1:length(sel3_phase)
            if isempty(sel3_phase{n}) == 0
                patch([tz(min(sel3_phase{n})) tz(min(sel3_phase{n})), tz(max(sel3_phase{n})) tz(max(sel3_phase{n}))],[-1 1 1 -1],[0.8 0.8 0.8],'EdgeColor','none','HandleVisibility','off');
            end
        end
        hold on
        plot(t,selLFP_HPC{i}, 'k', 'LineWidth',0.4)
        plot(t,selLFP_HPC_tfilt{i}, 'r', 'LineWidth', 1)
        for j=1:length(selSPK{i})
            line([selSPK{i}(j)+timePoints(i)-Ta selSPK{i}(j)+timePoints(i)-Ta], [0.8 1])
        end
        hold off
        title('HPC')
        xlabel('Time (s)')
        ylabel('Voltage (mV)')
        axis([timePoints(i)-Ta timePoints(i)+Tb -1 1])
        legend('Raw data','Theta-filtered data','Location','southwest')
        
        clear fphit
        clear sel_phase
        clear selsel_phase
        clear sel3_phase
end

% % delta
% figure
% for i = 1:length(timePoints)
%     t = (timePoints(i)-Ta+dt:dt:timePoints(i)+Tb);
%     subplot(2,length(timePoints),i)
%         plot(t,selLFP_PFC{i}, '.','MarkerSize',4)
%         hold on
%         plot(t,selLFP_PFC_dfilt{i}, 'r', 'LineWidth', 2)
%         for ii=1:length(selSPK{i})
%             line([selSPK{i}(ii)+timePoints(i)-Ta selSPK{i}(ii)+timePoints(i)-Ta], [0.8 1])
%         end
%         hold off
%         title('PFC')
%         xlabel('Time (s)')
%         ylabel('Voltage (mV)')
%         axis([timePoints(i)-Ta timePoints(i)+Tb -1 1])
%         legend('Raw data','Delta-filtered data','Location','southwest')
%     subplot(2,length(timePoints),i+length(timePoints))
%         plot(t,selLFP_HPC{i}, '.','MarkerSize',4)
%         hold on
%         plot(t,selLFP_HPC_dfilt{i}, 'r', 'LineWidth', 2)
%         for ii=1:length(selSPK{i})
%             line([selSPK{i}(ii)+timePoints(i)-Ta selSPK{i}(ii)+timePoints(i)-Ta], [0.8 1])
%         end
%         hold off
%         title('HPC')
%         xlabel('Time (s)')
%         ylabel('Voltage (mV)')
%         axis([timePoints(i)-Ta timePoints(i)+Tb -1 1])
%         legend('Raw data','Delta-filtered data','Location','southwest')
% end

%% calculating phase and plotting for each timepoint (HPC)

% theta, HPC, figure(2)
figure
for i=1:length(timePoints) %to divide subplots into 2 rows
    subplot(2,round(length(timePoints)/2),i)
    t = (timePoints(i)-Ta+dt:dt:timePoints(i)+Tb); %create t index
    selfphit{1,i}=angle(hilbert(selLFP_HPC_tfilt{1,i})); %hilbert transform for phase
    plot(t, selfphit{1,i}, 'k:');
    hold on
        for j=1:length(selSPK{i})
            seltSpk{i}(j) = find((abs(t-(selSPK{i}(j)+timePoints(i)-Ta)) < 0.0008),1); %grabs t indices of spikes
        end
    selSpkPhase{i} = selfphit{i}(seltSpk{i}); %grabs phase of each spike
    plot(t(seltSpk{i}),selSpkPhase{i},'o') %plots spikes
    hold off
    xlim([timePoints(i)-Ta timePoints(i)+Tb])
    chr = int2str(i); %for graph titles
    title(strcat('HPC Timepoint ',chr));
    xlabel('Time (s)');
    ylabel('Theta Phase (radians)');
end

% %delta
% figure
% for i=1:length(timePoints) %to divide subplots into 2 rows
%     subplot(2,round(length(timePoints)/2),i)
%     t = (timePoints(i)-Ta+dt:dt:timePoints(i)+Tb); %create t index
%     selfphid{1,i}=angle(hilbert(selLFP_HPC_dfilt{1,i})); %hilbert transform for phase
%     plot(t, selfphid{1,i}, 'k:');
%     hold on
%         for ii=1:length(selSPK{i})
%             seltSpk{i}(ii) = find((abs(t-(selSPK{i}(ii)+timePoints(i)-Ta)) < 0.0008),1); %grabs t indices of spikes
%         end
%     selSpkPhase{i} = selfphid{i}(seltSpk{i}); %grabs phase of each spike
%     plot(t(seltSpk{i}),selSpkPhase{i},'o') %plots spikes
%     hold off
%     xlim([timePoints(i)-Ta timePoints(i)+Tb])
%     chr = int2str(i); %for graph titles
%     title(strcat('Timepoint',chr));
%     xlabel('Time (s)');
%     ylabel('HPC Delta Phase (radians)');
% end

%% plotting phase and histogram for whole trace

% theta. figure(3)
figure
subplot(2,1,1)
t = dt:dt:(dt*length(LFP_HPC));
fphit=angle(hilbert(LFP_HPC_tfilt));
plot(t,fphit,'k:')
hold on
for i=1:length(SPK)
    tSpk(i) = find((abs(t-SPK(i)) < 0.0008),1); %grabs t indices of spikes
end
spkPhase = fphit(tSpk);
plot(t(tSpk),spkPhase,'o', 'MarkerEdgeColor','r')
hold off
title('HPC Phase Plot')
xlabel('Time (s)')
ylabel('Theta Phase (radians)')

% % histogram
% subplot(2,1,2)
% degSpkPhase = spkPhase .* (180/pi); %convert phase to degrees
% degSpkPhase = degSpkPhase + 180; %start at 0
% degSpkPhase2 = degSpkPhase + 360;
% degSpkPhase = cat(1,degSpkPhase,degSpkPhase2); %0 to 720
% histogram(degSpkPhase, 36, 'Normalization', 'probability')
% xlim([0 720])
% xticks([0:180:720])
% yticklabels(yticks*100)
% xlabel('Hippocampal Theta Phase (degrees)')
% ylabel('Percentage of Spikes')

%histogram of only the spikes from the selected timestamps
subplot(2,1,2)
spkPhase = selSpkPhase{:}; %combine all of the timepoint spike data
degSpkPhase = spkPhase .* (180/pi); %convert phase to degrees
degSpkPhase = degSpkPhase + 180; %start at 0
degSpkPhase2 = degSpkPhase + 360;
degSpkPhase = cat(1,degSpkPhase,degSpkPhase2); %0 to 720
histogram(degSpkPhase, 36, 'Normalization', 'probability')
xlim([0 720])
xticks([0:180:720])
yticklabels(yticks*100)
title('HPC Phase-distribution Histogram for spikes')
xlabel('Theta Phase (degrees)')
ylabel('Percentage of Spikes')
% 
% % delta
% figure
% subplot(2,1,1)
% t = dt:dt:(dt*length(LFP_HPC));
% fphid=angle(hilbert(LFP_HPC_dfilt));
% plot(t,fphid,'k:')
% hold on
% for i=1:length(SPK)
%     tSpk(i) = find((abs(t-SPK(i)) < 0.0008),1); %grabs t indices of spikes
% end
% spkPhase = fphid(tSpk);
% plot(t(tSpk),spkPhase,'o', 'MarkerEdgeColor','r')
% hold off
% title('Whole Trace')
% xlabel('Time (s)')
% ylabel('Hippocampal Delta Phase (radians)')
% 
% 
% %histogram of only the spikes from the selected timestamps
% subplot(2,1,2)
% spkPhase = selSpkPhase{:}; %combine all of the timepoint spike data
% degSpkPhase = spkPhase .* (180/pi); %convert phase to degrees
% degSpkPhase = degSpkPhase + 180; %start at 0
% degSpkPhase2 = degSpkPhase + 360;
% degSpkPhase = cat(1,degSpkPhase,degSpkPhase2); %0 to 720
% histogram(degSpkPhase, 36, 'Normalization', 'probability')
% xlim([0 720])
% xticks([0:180:720])
% yticklabels(yticks*100)
% xlabel('Hippocampal Delta Phase (degrees)')
% ylabel('Percentage of Spikes')
% title('Selected Timepoints')
% xlabel('Hippocampal Delta Phase (degrees)')
% ylabel('Percentage of Spikes')

%% calculating phase and plotting for each timepoint (PFC)

% theta, PFC, figure(4)
figure
for i=1:length(timePoints) %to divide subplots into 2 rows
    subplot(2,round(length(timePoints)/2),i)
    t = (timePoints(i)-Ta+dt:dt:timePoints(i)+Tb); %create t index
    selfphit{1,i}=angle(hilbert(selLFP_PFC_tfilt{1,i})); %hilbert transform for phase
    plot(t, selfphit{1,i}, 'k:');
    hold on
        for j=1:length(selSPK{i})
            seltSpk{i}(j) = find((abs(t-(selSPK{i}(j)+timePoints(i)-Ta)) < 0.0008),1); %grabs t indices of spikes
        end
    selSpkPhase{i} = selfphit{i}(seltSpk{i}); %grabs phase of each spike
    plot(t(seltSpk{i}),selSpkPhase{i},'o') %plots spikes
    hold off
    xlim([timePoints(i)-Ta timePoints(i)+Tb])
    chr = int2str(i); %for graph titles
    title(strcat('PFC Timepoint ',chr));
    xlabel('Time (s)');
    ylabel('Theta Phase (radians)');
end


% theta, PFC, figure(5)
figure
subplot(2,1,1)
t = dt:dt:(dt*length(LFP_PFC));
fphit=angle(hilbert(LFP_PFC_tfilt));
plot(t,fphit,'k:')
hold on
for i=1:length(SPK)
    tSpk(i) = find((abs(t-SPK(i)) < 0.0008),1); %grabs t indices of spikes
end
spkPhase = fphit(tSpk);
plot(t(tSpk),spkPhase,'o', 'MarkerEdgeColor','r')
hold off
title('PFC Phase Plot')
xlabel('Time (s)')
ylabel('Theta Phase (radians)')

%histogram of only the spikes from the selected timestamps
subplot(2,1,2)
spkPhase = selSpkPhase{:}; %combine all of the timepoint spike data
degSpkPhase = spkPhase .* (180/pi); %convert phase to degrees
degSpkPhase = degSpkPhase + 180; %start at 0
degSpkPhase2 = degSpkPhase + 360;
degSpkPhase = cat(1,degSpkPhase,degSpkPhase2); %0 to 720
histogram(degSpkPhase, 36, 'Normalization', 'probability')
xlim([0 720])
xticks([0:180:720])
yticklabels(yticks*100)
title('PFC Phase-distribution Histogram for spikes')
xlabel('Theta Phase (degrees)')
ylabel('Percentage of Spikes')


%% coherence

% Parameters setting for coherence
ParaCoh.tapers = [3 5];
ParaCoh.pad = 0;
ParaCoh.Fs = Fs;
ParaCoh.fpass = [0 100];

% main estimation
for i=1:length(selLFP_HPC)
    [C,~,~,~,~,f] = coherencyc(selLFP_PFC_100filt{i},selLFP_HPC_100filt{i},ParaCoh);
    coh(i,:) = C;
end

windowWidth = 15;
kernel = ones(windowWidth,1) / windowWidth;
meanCoh = mean(coh);
filtmeanCoh = filter(kernel,1,meanCoh);
SEMcoh = (std(coh)/sqrt(length(coh(:,1))));

%% figure(6)
figure
f2 = [f,fliplr(f)];
inBetweenPos = [filtmeanCoh, fliplr(filtmeanCoh+SEMcoh)];
fill(f2,inBetweenPos,[1 0.7 0.7]);
hold on;
inBetweenNeg = [filtmeanCoh-SEMcoh, fliplr(filtmeanCoh)];
fill(f2,inBetweenNeg,[1 0.7 0.7]);
hold on

plot(f,filtmeanCoh+SEMcoh,'color',[1 0.7 0.7])
hold on
plot(f,filtmeanCoh-SEMcoh,'color',[1 0.7 0.7])
hold on

plot(f,filtmeanCoh,'r','LineWidth',1.5)
title('Coherence between PFC and HPC')
xlabel('Frequency (Hz)')
ylabel('Coherence')
hold off

%% figure(7)
figure

cohDelta = find(f>=0.5);
cohDelta = cohDelta(find(f(cohDelta)<=4));
cohDelta = meanCoh(cohDelta);

cohWaves(1,:) = mean(cohDelta);
semWaves(1,:) = std(cohDelta)./sqrt(size(cohDelta,2));

cohTheta = find(f>=4);
cohTheta = cohTheta(find(f(cohTheta)<=12));
cohTheta = meanCoh(cohTheta);

cohWaves(2,:) = mean(cohTheta);
semWaves(2,:) = std(cohTheta)./sqrt(size(cohTheta,2));

cohBeta = find(f>=12);
cohBeta = cohBeta(find(f(cohBeta)<=30));
cohBeta = meanCoh(cohBeta);

cohWaves(3,:) = mean(cohBeta);
semWaves(3,:) = std(cohBeta)./sqrt(size(cohBeta,2));


cohGamma = find(f>=30);
cohGamma = cohGamma(find(f(cohGamma)<=100));
cohGamma = meanCoh(cohGamma);

cohWaves(4,:) = mean(cohGamma);
semWaves(4,:) = std(cohGamma)./sqrt(size(cohGamma,2));

bar(cohWaves)
hold on
for i=1:4
    errorbar(i,cohWaves(i),semWaves(i),'k','LineWidth',1.5)
end
ylabel('Coherence')
xticklabels({'Delta','Theta','Beta','Gamma'})